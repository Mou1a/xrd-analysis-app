import streamlit as st
import io
import re
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

import warnings
warnings.filterwarnings("ignore")

# ── STREAMLIT PAGE CONFIG ──────────────────────────────────────────────────────
st.set_page_config(page_title="ICP-OES Data Processor", page_icon="🧪", layout="wide")

# ── constants ──────────────────────────────────────────────────────────────────
MC_ELEMENTS = {"Al", "Ca", "Fe", "K", "Mg", "Na", "S", "Ti"}

# Colours
C_NAVY = "1F3864"; C_BLUE = "2F5496"; C_GREEN = "1F5C1F"; C_LGRN = "E2EFDA"; C_DGRN = "1A5C1A"; C_MGRN = "2D7A2D"
C_PURP = "6D2C6D"; C_LPURP = "F3E0F7"; C_MPURP = "8E3D8E"; C_BRWN = "7B3000"; C_LBRN = "EAD1DC"; C_MBRN = "A04000"
C_YEL = "FFF2CC"; C_WHITE = "FFFFFF"; C_GOOD = "C6EFCE"; C_WARN = "FFEB9C"; C_BAD = "FFC7CE"; C_LGREY = "F5F5F5"

# ═══════════════════════════════════════════════════════════════════════════════
# 1. HEURISTIC PARSING & UTILS
# ═══════════════════════════════════════════════════════════════════════════════

def clean_number(x):
    if pd.isna(x): return None
    if isinstance(x, (int, float)): 
        return float(x) if not np.isnan(x) and not np.isinf(x) else None
        
    s = str(x).strip()
    if s in ("--", "####", "", "Uncal"): return None
    if re.search(r"\bu\b|u\s*$", s, re.IGNORECASE) or s.startswith('<'): return 0.0
    if re.search(r"\bo\b|o\s*$", s, re.IGNORECASE) or s.startswith('>'): return None
        
    try: return float(re.sub(r"[a-zA-Z\s]+$", "", s).strip())
    except ValueError: return None

def clean_scalar(x):
    if isinstance(x, pd.Series): x = x.iloc[0] if len(x) > 0 else None
    return clean_number(x)

def normalize(col): 
    return str(col).replace(" ppm", "").replace(" c/s", "").strip()

def clean_col_name(c):
    return re.sub(r'\s*(ppm|mg/l|ppb|c/s|cps)\s*$', '', str(c), flags=re.I).strip().upper()

def extract_element(line_name):
    m = re.match(r"^([A-Z][a-z]?)\s+", str(line_name).strip())
    return m.group(1) if m else str(line_name).split()[0]

def extract_all_data_blocks(raw):
    blocks, header_rows = [], []
    for i in range(len(raw)):
        row_strs = raw.iloc[i].astype(str)
        nm_count = row_strs.str.contains(r"nm", flags=re.IGNORECASE).sum()
        if nm_count >= 5:
            if not header_rows or i - header_rows[-1] > 2:
                header_rows.append(i)
                
    for idx, h_row in enumerate(header_rows):
        end_idx = header_rows[idx+1] if idx + 1 < len(header_rows) else len(raw)
        df = raw.iloc[h_row:end_idx].copy()
        
        headers = [str(h).strip() if pd.notna(h) and str(h).strip() != "" else f"Unnamed_{i}" for i, h in enumerate(df.iloc[0].tolist())]
        df.columns = headers
        df = df.iloc[1:].reset_index(drop=True)
        df = df.loc[:, ~df.columns.duplicated()].copy()
        
        df = df.dropna(how="all")
        if "Solution Label" in df.columns:
            df = df[df["Solution Label"].astype(str).str.strip() != "Solution Label"]
        valid_rows = df.astype(str).replace(r'^\s*$', np.nan, regex=True).notna().sum(axis=1) > 2
        df = df[valid_rows].reset_index(drop=True)
        if not df.empty: blocks.append(df)
            
    return blocks

def get_block_numeric_mean(df):
    nm_cols = [c for c in df.columns if "nm" in str(c).lower()]
    total_sum, count = 0, 0
    for col in nm_cols:
        s = df[col].apply(clean_number).dropna()
        total_sum += s.sum()
        count += len(s)
    return total_sum / count if count > 0 else 0

def identify_blocks(blocks):
    intens_df, conc_df, unadj_df = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    if not blocks: return conc_df, intens_df, unadj_df
    
    stats = [{"df": df, "mean": get_block_numeric_mean(df)} for df in blocks]
    stats.sort(key=lambda x: x["mean"], reverse=True)
    
    if len(stats) >= 1:
        if len(stats) == 1:
            conc_df = stats[0]["df"]
        elif stats[0]["mean"] > 5000 or (len(stats) > 1 and stats[0]["mean"] > 10 * stats[1]["mean"]):
            intens_df = stats[0]["df"]
            remaining = stats[1:]
            if len(remaining) == 1:
                conc_df = remaining[0]["df"]
            elif len(remaining) >= 2:
                conc_df = remaining[0]["df"]
                unadj_df = remaining[1]["df"]
        else:
            conc_df = stats[0]["df"]
            unadj_df = stats[1]["df"]
            
    return conc_df, intens_df, unadj_df

def get_calibrated_lines(conc_df):
    calibrated, discarded = set(), set()
    for col in conc_df.columns[2:]:
        col_name = normalize(col)
        uncal_count = conc_df[col].astype(str).str.contains("Uncal", case=False, na=False).sum()
        if uncal_count >= 3: discarded.add(col_name)
        else: calibrated.add(col_name)
    return calibrated, discarded

def parse_qc_label(label):
    s = str(label).strip()
    if "QC" not in s.upper(): return None, None
    # Smarter REGEX: Matches "QC MC 10" OR "MC QC 10"
    m = re.search(r"(?:QC\s*([A-Za-z]+)|([A-Za-z]+)\s*QC)\s*([0-9.]+)", s, re.I)
    if not m: return "Unknown", None
    
    # Extract the Type (MC or Si) regardless of which side of "QC" it was on
    qc_type = m.group(1) if m.group(1) else m.group(2)
    return qc_type, float(m.group(3))

def closest_target(value, targets):
    v = clean_scalar(value)
    if v is None or not targets: return None
    return min(targets, key=lambda t: abs(float(v) - float(t)))

def get_base_name(x):
    s = str(x).strip()
    has_f = bool(re.search(r"[\s_]+f$", s, re.IGNORECASE))
    s = re.sub(r"[\s_]+f$", "", s, flags=re.IGNORECASE)
    s = re.sub(r'(?<=\d)[a-eA-E](?=\s|_|$)', '', s)
    s = re.sub(r'[\s\-_]+[a-eA-E](?=\s|_|$)', '', s)
    s = re.sub(r'(?<=[a-zA-Z])[A-E](?=\s|_|$)', '', s)
    s = re.sub(r'[\s\-_]*rep\d*(?=\s|_|$)', '', s, flags=re.IGNORECASE)
    s = re.sub(r'\s+', ' ', s).strip()
    if has_f: s += "_f"
    return s

def format_val(val):
    if not isinstance(val, (int, float)) or pd.isna(val) or np.isinf(val): return "N/A"
    tgt = abs(val)
    if tgt < 10: return f"{val:,.2f}"
    elif tgt < 100: return f"{val:,.1f}"
    else: return f"{val:,.0f}"

# ═══════════════════════════════════════════════════════════════════════════════
# 2. DATA PROCESSING STEPS
# ═══════════════════════════════════════════════════════════════════════════════

def step1_contamination(conc_df, intens_df, cal_lines):
    if intens_df.empty: return pd.DataFrame()
    
    lc = "Solution Label"
    results = []
    blank_mask = intens_df[lc].astype(str).str.contains(r"\bBlank\b|^B$", case=False, na=False)
    if not blank_mask.any(): return pd.DataFrame()
    blank_row = intens_df[blank_mask].iloc[0]

    conc_stds = conc_df[conc_df[lc].astype(str).str.contains(r"Standard|STD", case=False, na=False, regex=True)].copy()
    int_stds  = intens_df[intens_df[lc].astype(str).str.contains(r"Standard|STD", case=False, na=False, regex=True)].copy()

    conc_cols = {clean_col_name(c): c for c in conc_df.columns[2:]}
    int_cols  = {clean_col_name(c): c for c in intens_df.columns[2:]}

    for base_key, int_col in int_cols.items():
        if base_key not in conc_cols: continue
        conc_col = conc_cols[base_key]; line_key = normalize(conc_col)
        if line_key not in cal_lines: continue
        
        element = extract_element(line_key)
        sel_label, sel_conc, sel_intens = None, None, None
        
        for _, row in conc_stds.iterrows():
            exp = clean_number(row[conc_col])
            if exp is not None:
                sel_label = str(row[lc]).strip(); sel_conc = exp
                mi = int_stds[int_stds[lc].astype(str).str.strip() == sel_label]
                if not mi.empty: sel_intens = clean_number(mi.iloc[0][int_col])
                break
        
        blank_i = clean_number(blank_row[int_col])
        neg_inv = False
        
        if blank_i is None or sel_intens is None or sel_intens == 0:
            ratio = None; flag = "Check"; neg_inv = True
        elif sel_intens < 0:
            ratio = None; flag = "YES - negative std"; neg_inv = True
        else:
            ratio = blank_i / sel_intens
            flag = "YES" if ratio > 0.10 else "NO"
            
        est_b = (sel_conc * ratio) if (ratio is not None and sel_conc is not None) else None
        
        results.append({
            "Element": element, "Line": line_key, "Selected standard": sel_label,
            "Std expected conc (ppm)": sel_conc, "Blank intensity": blank_i, "Std intensity": sel_intens,
            "Blank/Std intensity ratio": ratio, "Contamination flag >10%": flag,
            "Est. blank conc (ppm)": est_b, "Detection limit = 2x est. blank (ppm)": (2 * est_b) if est_b else None, "_rank_penalty": 1 if neg_inv else 0
        })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values(["Element", "_rank_penalty", "Blank/Std intensity ratio"], na_position='last')
        df = df.drop(columns=["_rank_penalty"]).reset_index(drop=True)
    return df

def step2_qc_statistics(conc_df, unadj_df, cal_lines):
    lc = "Solution Label"
    results = []
    main = conc_df.reset_index(drop=True)
    
    unadj_occ = {}
    if not unadj_df.empty:
        for idx, r in unadj_df.iterrows():
            l = str(r.get(lc, "")).strip()
            unadj_occ.setdefault(l, []).append(idx)
            
    conc_occ_seen = {}
    
    for i in range(len(main)):
        row = main.iloc[i]
        label_raw = row.get(lc, "")
        label = str(label_raw).strip()
        
        conc_occ_seen[label] = conc_occ_seen.get(label, 0) + 1
        occ_idx = conc_occ_seen[label] - 1
        
        if "QC" not in str(label_raw).upper(): continue
            
        qc_type, target = parse_qc_label(label_raw)
        qc_type_upper = str(qc_type).upper() if qc_type else "UNKNOWN"
        
        prev_row = main.iloc[i-1] if i > 0 else None
        prev_label = prev_row[lc] if prev_row is not None else None
        prev_is_blank = prev_label is not None and bool(re.search(r"\bBLANK\b|^B$", str(prev_label).strip(), re.IGNORECASE))

        u_idx = None
        if label in unadj_occ:
            if occ_idx < len(unadj_occ[label]):
                u_idx = unadj_occ[label][occ_idx]
            else:
                u_idx = unadj_occ[label][-1]

        for col in main.columns[2:]:
            line_key = normalize(col); base_col_key = clean_col_name(line_key); elem = extract_element(line_key)
            if line_key not in cal_lines: continue
            if qc_type_upper == "MC" and elem not in MC_ELEMENTS: continue
            if qc_type_upper == "SI" and elem != "Si": continue
            if qc_type_upper in ["CL", "C"] and elem != "Cl": continue
            
            m_unadj = None
            if u_idx is not None and not unadj_df.empty:
                unadj_col = next((uc for uc in unadj_df.columns if clean_col_name(uc) == base_col_key), None)
                if unadj_col:
                    m_unadj = clean_scalar(unadj_df.iloc[u_idx][unadj_col])
            
            if m_unadj is None: m_unadj = clean_number(row[col])
            bc = clean_number(prev_row[col]) if prev_is_blank else None
            err = ((m_unadj - target) / target * 100) if (m_unadj is not None and target) else None
            bratio = (bc / target * 100) if (bc is not None and target) else None
            
            results.append({
                "QC label": label_raw, "QC type": qc_type, "Target concentration": target, 
                "Element": elem, "Line": line_key, "Unadjusted Measured": m_unadj, 
                "Recovery %": (m_unadj / target * 100) if (m_unadj and target) else None,
                "Error %": err, "QC flag >10% error": "YES" if (err and abs(err) > 10) else "NO",
                "Previous row label": prev_label, "Blank before QC conc": bc, 
                "Blank before QC / target %": bratio,
                "Blank before QC flag >10%": "YES" if (bratio and abs(bratio) > 10) else "NO" if prev_is_blank else "No blank immediately before QC"
            })
    return pd.DataFrame(results)

def step3_qc_summary(qc_df):
    if qc_df.empty: return pd.DataFrame(), pd.DataFrame()
    d = qc_df.copy()
    d["Abs_error_pct"] = pd.to_numeric(d["Error %"], errors="coerce").abs()
    
    summ = d.groupby(["Element", "Target concentration", "Line"], as_index=False).agg(
        Avg_error_pct=("Error %", "mean"), Avg_abs_error_pct=("Abs_error_pct", "mean"),
        Std_error_pct=("Error %", "std"), N=("Error %", "count")
    )
    summ["_sort_err"] = summ["Avg_abs_error_pct"].fillna(999999)
    summ["Rank"] = summ.groupby(["Element", "Target concentration"])["_sort_err"].rank(method="first").astype(int)
    
    best = summ[summ["Rank"] == 1].copy()
    return summ, best

def step4_data_overview(conc_df, unadj_df, best_lines_df, contam_df, summary_df):
    lc = "Solution Label"
    if summary_df.empty: return {}, [], pd.DataFrame(), pd.DataFrame()
    
    uni_map = {}
    for e, g in summary_df.groupby("Element"):
        means = g.groupby("Line")["Avg_abs_error_pct"].mean().dropna()
        uni_map[e] = means.idxmin() if not means.empty else g["Line"].iloc[0]

    elements = sorted(uni_map.keys())
    targets = sorted(pd.to_numeric(best_lines_df["Target concentration"], errors='coerce').dropna().unique())
    base_dl_map = {r["Line"]: r["Detection limit = 2x est. blank (ppm)"] for _, r in contam_df.iterrows()} if not contam_df.empty else {}
    
    pub_rows = []
    main = conc_df.reset_index(drop=True)
    
    valid_label_mask = main[lc].astype(str).str.strip() != ""
    not_qc_mask = ~main[lc].astype(str).str.upper().str.contains(r"^(QC|STANDARD|STD)|\bBLANK\b|^B$", case=False, na=False)
    s_idx = main.index[valid_label_mask & not_qc_mask].tolist()
    
    unadj_occ = {}
    if not unadj_df.empty:
        for idx, r in unadj_df.iterrows():
            l = str(r.get(lc, "")).strip()
            unadj_occ.setdefault(l, []).append(idx)
            
    conc_occ_seen = {}
    
    for i in range(len(main)):
        row = main.iloc[i]
        label = str(row.get(lc, "")).strip()
        
        conc_occ_seen[label] = conc_occ_seen.get(label, 0) + 1
        occ_idx = conc_occ_seen[label] - 1
        
        if i not in s_idx: continue
            
        u_idx = None
        if label in unadj_occ:
            if occ_idx < len(unadj_occ[label]): u_idx = unadj_occ[label][occ_idx]
            else: u_idx = unadj_occ[label][-1]

        df_ratios = []
        if u_idx is not None and not unadj_df.empty:
            for col in main.columns[2:]:
                base_col_key = clean_col_name(normalize(col))
                unadj_col = next((uc for uc in unadj_df.columns if clean_col_name(uc) == base_col_key), None)
                if unadj_col:
                    adj_v = clean_number(main.iloc[i][col])
                    unadj_v = clean_number(unadj_df.iloc[u_idx][unadj_col])
                    if adj_v is not None and unadj_v is not None and adj_v > 0.0001 and unadj_v > 0.0001:
                        df_ratios.append(adj_v / unadj_v)
        
        DF = np.median(df_ratios) if df_ratios else None
        
        if DF is None or pd.isna(DF):
            m_dil = re.search(r'1:(\d+)', label)
            DF = float(m_dil.group(1)) if m_dil else 1.0
            
        if DF < 0.9: DF = 1.0
        if abs(DF - round(DF)) / DF < 0.05: DF = round(DF)
        else: DF = round(DF, 2)

        out = {"Sample": label, "_DF": DF}
        for elem in elements:
            uv = None
            if u_idx is not None and not unadj_df.empty:
                for uc in unadj_df.columns:
                    if extract_element(normalize(uc)) == elem:
                        uv = clean_scalar(unadj_df.iloc[u_idx][uc])
                        if uv is not None: break
            
            if uv is None:
                for ac in main.columns:
                    if extract_element(normalize(ac)) == elem:
                        av = clean_scalar(main.iloc[i][ac])
                        if av is not None:
                            uv = av / DF
                            break
            
            t = closest_target(uv, targets)
            
            ranked_lines = []
            if t is not None:
                sub = summary_df[(summary_df["Element"] == elem) & (summary_df["Target concentration"] == t)].sort_values("Rank")
                ranked_lines = sub["Line"].apply(normalize).tolist()
            
            uni_sub = summary_df[summary_df["Element"] == elem].groupby("Line")["Avg_abs_error_pct"].mean().sort_values()
            uni_ranked_lines = uni_sub.index.map(normalize).tolist()
            
            lines_to_try, seen = [], set()
            for l in ranked_lines + uni_ranked_lines:
                if l not in seen:
                    lines_to_try.append(l); seen.add(l)
            
            best_ln = lines_to_try[0] if lines_to_try else None
            val = None
            
            for ln in lines_to_try:
                adj_col = next((c for c in main.columns if normalize(c) == ln), None)
                v = clean_scalar(main.iloc[i][adj_col]) if adj_col else None
                if v is not None:
                    best_ln = ln; val = v
                    break
            
            err_pct = None
            if t is not None:
                err_row = summary_df[(summary_df["Element"] == elem) & (summary_df["Target concentration"] == t) & (summary_df["Line"] == best_ln)]
                if not err_row.empty: err_pct = err_row.iloc[0]["Avg_error_pct"]
            if err_pct is None or pd.isna(err_pct):
                err_row = summary_df[(summary_df["Element"] == elem) & (summary_df["Line"] == best_ln)]
                if not err_row.empty: err_pct = err_row["Avg_error_pct"].mean()
            
            base_dl = base_dl_map.get(best_ln)
            adj_dl = (base_dl * DF) if base_dl is not None else None
            
            out[f"{elem} (ppm)"] = val
            out[f"{elem}_line"] = best_ln
            out[f"{elem}_DL"] = adj_dl
            out[f"{elem}_Err"] = err_pct
            
        pub_rows.append(out)
    
    pub_df = pd.DataFrame(pub_rows)
    rep_df = pd.DataFrame()
    if not pub_df.empty:
        pdf = pub_df.copy(); pdf["_base"] = pdf["Sample"].apply(get_base_name)
        r_rows = []
        for base, grp in pdf.groupby("_base", sort=False):
            if not base or len(grp) < 2: continue
            r = {"Sample (base)": base, "N": len(grp)}
            for e in elements:
                v = pd.to_numeric(grp[f"{e} (ppm)"], errors='coerce').dropna()
                mean_val = v.mean() if not v.empty else None
                r[f"{e} Mean"] = mean_val
                r[f"{e} SD"] = v.std() if len(v) >= 2 else None
                if len(v) >= 2 and mean_val and mean_val != 0: r[f"{e} RSD%"] = (v.std() / mean_val) * 100
                else: r[f"{e} RSD%"] = None
                
                r[f"{e} DL"] = grp[f"{e}_DL"].iloc[0]
                r[f"{e} Err"] = grp[f"{e}_Err"].iloc[0]
            r_rows.append(r)
        rep_df = pd.DataFrame(r_rows)
    return uni_map, elements, pub_df, rep_df

# ═══════════════════════════════════════════════════════════════════════════════
# 3. EXCEL WRITER (SINGLE SHEET, MULTI-BLOCK LAYOUT)
# ═══════════════════════════════════════════════════════════════════════════════

def _cell(ws, r, c, val="", bg="FFFFFF", fc="000000", bold=False, ha="center"):
    fmt = None
    if isinstance(val, (int, float)) and not pd.isna(val) and not np.isinf(val):
        tgt = abs(val)
        if tgt < 10: fmt = "#,##0.00"
        elif tgt < 100: fmt = "#,##0.0"
        else: fmt = "#,##0"
        
    if val is None or pd.isna(val) or (isinstance(val, float) and (np.isnan(val) or np.isinf(val))): 
        val = ""
        
    cell = ws.cell(row=r, column=c, value=val)
    cell.font, cell.fill = Font(bold=bold, color=fc, size=10, name="Arial"), PatternFill("solid", start_color=bg)
    cell.alignment = Alignment(horizontal=ha, vertical="center", wrap_text=True)
    cell.border = Border(left=Side(style='thin', color="BBBBBB"), right=Side(style='thin', color="BBBBBB"), top=Side(style='thin', color="BBBBBB"), bottom=Side(style='thin', color="BBBBBB"))
    
    if fmt: cell.number_format = fmt
    return cell

def _sec(ws, r, c, text, ncols, bg):
    ws.merge_cells(start_row=r, start_column=c, end_row=r, end_column=c + ncols - 1)
    cell = ws.cell(row=r, column=c, value=text)
    cell.font, cell.fill, cell.alignment = Font(bold=True, color="FFFFFF", size=11, name="Arial"), PatternFill("solid", start_color=bg), Alignment(horizontal="center")

def _sub_sec(ws, r, c, text, ncols, bg):
    ws.merge_cells(start_row=r, start_column=c, end_row=r, end_column=c + ncols - 1)
    cell = ws.cell(row=r, column=c, value=text)
    cell.font, cell.fill, cell.alignment = Font(bold=True, color="FFFFFF", size=10, name="Arial"), PatternFill("solid", start_color=bg), Alignment(horizontal="center")

def write_sheet(wb, s_name, elements, pub_df, rep_df, summ_df, qc_df, contam_df, uni_map):
    ws = wb.create_sheet(title=s_name[:31]); ws.sheet_view.showGridLines = False
    SR, HR, ca = 2, 3, 2; ne = len(elements); no = 1 + ne
    
    _sec(ws, SR, ca, "DATA OVERVIEW (Best-Line Conc. ppm)", no, C_GREEN)
    _cell(ws, HR, ca, "Sample", bg=C_GREEN, fc="FFFFFF", bold=True)
    ws.column_dimensions[get_column_letter(ca)].width = 24
    for j, e in enumerate(elements): 
        _cell(ws, HR, ca+1+j, e, bg=C_MGRN, fc="FFFFFF", bold=True)
        ws.column_dimensions[get_column_letter(ca+1+j)].width = 13

    for i, (_, row) in enumerate(pub_df.iterrows()):
        r, bg = HR+1+i, (C_LGRN if i%2==0 else C_WHITE)
        _cell(ws, r, ca, row["Sample"], bg=bg, ha="left", bold=True)
        
        for j, e in enumerate(elements): 
            val = row[f"{e} (ppm)"]
            dl_val = row.get(f'{e}_DL')
            err_val = row.get(f'{e}_Err')
            
            cell_bg = C_BAD if isinstance(err_val, (int, float)) and pd.notna(err_val) and abs(err_val) > 10 else bg
            
            display_val = val
            if isinstance(val, (int, float)) and pd.notna(val):
                if isinstance(dl_val, (int, float)) and pd.notna(dl_val) and dl_val > 0:
                    if val <= dl_val or val == 0.0:
                        display_val = f"< {format_val(dl_val)}"
                elif val == 0.0:
                    display_val = "< DL"
                    
            _cell(ws, r, ca+1+j, display_val, bg=cell_bg)

    lr = HR + len(pub_df) + 2
    _sub_sec(ws, lr, ca, "Emission lines, DLs, & QC Error", no, "3A6A3A")
    for i, (_, row) in enumerate(pub_df.iterrows()):
        r, bg = lr+1+i, (C_LGREY if i%2==0 else C_WHITE)
        _cell(ws, r, ca, row["Sample"], bg=bg, ha="left")
        for j, e in enumerate(elements): 
            val_str = row.get(f'{e}_line')
            val_str = str(val_str) if val_str else ""
            
            dl_val = row.get(f'{e}_DL'); err_val = row.get(f'{e}_Err')
            dl_str = format_val(dl_val)
            err_str = f"{err_val:.1f}%" if isinstance(err_val, (int, float)) and not pd.isna(err_val) else "N/A"
            
            cell_bg = C_BAD if isinstance(err_val, (int, float)) and pd.notna(err_val) and abs(err_val) > 10 else bg
            
            val = f"{val_str}\n(DL: {dl_str}, Err: {err_str})" if val_str else ""
            _cell(ws, r, ca+1+j, val, bg=cell_bg).font = Font(size=8)

    if not rep_df.empty:
        curr = lr + len(pub_df) + 3
        for b_type in ["Mean", "SD", "RSD%"]:
            _sec(ws, curr, ca, f"REPLICATE {b_type.upper()}", 2+ne, C_DGRN); curr += 1
            _cell(ws, curr, ca, "Sample (base)", bg=C_MGRN, fc="FFFFFF"); _cell(ws, curr, ca+1, "N", bg=C_MGRN, fc="FFFFFF")
            for j, e in enumerate(elements): _cell(ws, curr, ca+2+j, e, bg=C_MGRN, fc="FFFFFF")
            for i, (_, row) in enumerate(rep_df.iterrows()):
                r, bg = curr+1+i, (C_LGRN if i%2==0 else C_WHITE)
                _cell(ws, r, ca, row["Sample (base)"], bg=bg, ha="left"); _cell(ws, r, ca+1, row["N"], bg=bg)
                for j, e in enumerate(elements): 
                    val = row[f"{e} {b_type}"]
                    dl_val = row.get(f'{e} DL')
                    err_val = row.get(f'{e} Err')
                    mean_val = row.get(f"{e} Mean")
                    
                    cell_bg = C_BAD if isinstance(err_val, (int, float)) and pd.notna(err_val) and abs(err_val) > 10 else bg
                    
                    display_val = val
                    is_below_dl = False
                    if isinstance(mean_val, (int, float)) and pd.notna(mean_val):
                        if isinstance(dl_val, (int, float)) and pd.notna(dl_val) and dl_val > 0:
                            if mean_val <= dl_val or mean_val == 0.0: is_below_dl = True
                        elif mean_val == 0.0: is_below_dl = True
                            
                    if is_below_dl:
                        if b_type == "Mean":
                            display_val = f"< {format_val(dl_val)}" if isinstance(dl_val, (int, float)) and pd.notna(dl_val) else "< DL"
                        else:
                            display_val = "" 
                            
                    _cell(ws, r, ca+2+j, display_val, bg=cell_bg)
            curr += len(rep_df) + 2
            
        _sub_sec(ws, curr, ca, "Replicate Emission lines, DLs, & QC Error", no, "3A6A3A")
        pdf_copy = pub_df.copy(); pdf_copy["_base"] = pdf_copy["Sample"].apply(get_base_name)
        for i, (_, row) in enumerate(rep_df.iterrows()):
            r, bg = curr+1+i, (C_LGREY if i%2==0 else C_WHITE)
            base_name = row["Sample (base)"]
            _cell(ws, r, ca, base_name, bg=bg, ha="left")
            
            match = pdf_copy[pdf_copy["_base"] == base_name]
            match_row = match.iloc[0] if not match.empty else pd.Series()
            
            for j, e in enumerate(elements): 
                val_str = match_row.get(f'{e}_line')
                val_str = str(val_str) if pd.notna(val_str) and val_str else ""
                dl_val = match_row.get(f'{e}_DL'); err_val = match_row.get(f'{e}_Err')
                
                dl_str = format_val(dl_val)
                err_str = f"{err_val:.1f}%" if isinstance(err_val, (int, float)) and not pd.isna(err_val) else "N/A"
                cell_bg = C_BAD if isinstance(err_val, (int, float)) and pd.notna(err_val) and abs(err_val) > 10 else bg
                
                val = f"{val_str}\n(DL: {dl_str}, Err: {err_str})" if val_str else ""
                _cell(ws, r, ca+1+j, val, bg=cell_bg).font = Font(size=8)
        curr += len(rep_df) + 2

        _sec(ws, curr, ca, "REPLICATE SUMMARY (Mean ± SD)", 2+ne, C_DGRN)
        _cell(ws, curr+1, ca, "Sample (base)", bg=C_MGRN, fc="FFFFFF")
        _cell(ws, curr+1, ca+1, "N", bg=C_MGRN, fc="FFFFFF")
        for j, e in enumerate(elements): 
            _cell(ws, curr+1, ca+2+j, e, bg=C_MGRN, fc="FFFFFF")
            
        for i, (_, row) in enumerate(rep_df.iterrows()):
            r, bg = curr+2+i, (C_LGRN if i%2==0 else C_WHITE)
            _cell(ws, r, ca, row["Sample (base)"], bg=bg, ha="left")
            _cell(ws, r, ca+1, row["N"], bg=bg)
            
            for j, e in enumerate(elements): 
                mean_val = row.get(f"{e} Mean")
                sd_val = row.get(f"{e} SD")
                dl_val = row.get(f"{e} DL")
                err_val = row.get(f"{e} Err")
                
                cell_bg = C_BAD if isinstance(err_val, (int, float)) and pd.notna(err_val) and abs(err_val) > 10 else bg
                
                display_val = ""
                is_below_dl = False
                if isinstance(mean_val, (int, float)) and pd.notna(mean_val):
                    if isinstance(dl_val, (int, float)) and pd.notna(dl_val) and dl_val > 0:
                        if mean_val <= dl_val or mean_val == 0.0: is_below_dl = True
                    elif mean_val == 0.0: is_below_dl = True
                        
                if is_below_dl:
                    display_val = f"< {format_val(dl_val)}" if isinstance(dl_val, (int, float)) and pd.notna(dl_val) else "< DL"
                else:
                    if isinstance(mean_val, (int, float)) and pd.notna(mean_val):
                        m_str = format_val(mean_val)
                        if isinstance(sd_val, (int, float)) and pd.notna(sd_val):
                            s_str = format_val(sd_val)
                            display_val = f"{m_str} ± {s_str}"
                        else:
                            display_val = m_str
                            
                _cell(ws, r, ca+2+j, display_val, bg=cell_bg)
        curr += len(rep_df) + 3

    cb = ca + no + 2; cols_b = ["Element", "Target concentration", "Line", "Avg_error_pct", "Avg_abs_error_pct", "Std_error_pct", "N", "Rank"]
    _sec(ws, SR, cb, "QC SUMMARY (Line Ranking)", len(cols_b), C_NAVY)
    for j, h in enumerate(cols_b): 
        _cell(ws, HR, cb+j, h, bg=C_BLUE, fc="FFFFFF", bold=True)
        ws.column_dimensions[get_column_letter(cb+j)].width = 16
    for i, (_, row) in enumerate(summ_df.iterrows()):
        r, bg = HR+1+i, (C_YEL if i%2==0 else C_WHITE)
        for j, h in enumerate(cols_b): _cell(ws, r, cb+j, row.get(h, ""), bg=bg)

    cc = cb + len(cols_b) + 2; cols_c = ["QC label", "Element", "Line", "Unadjusted Measured", "Target concentration", "Recovery %", "Error %", "QC flag >10% error", "Blank before QC / target %"]
    _sec(ws, SR, cc, "QC ANALYSIS (Data-Driven Link)", len(cols_c), C_PURP)
    for j, h in enumerate(cols_c): 
        _cell(ws, HR, cc+j, h, bg=C_MPURP, fc="FFFFFF", bold=True)
        ws.column_dimensions[get_column_letter(cc+j)].width = 18
    for i, (_, row) in enumerate(qc_df.iterrows()):
        r, bg = HR+1+i, (C_LPURP if i%2==0 else C_WHITE)
        for j, h in enumerate(cols_c): _cell(ws, r, cc+j, row.get(h, ""), bg=bg)

    cd = cc + len(cols_c) + 2; cols_d = ["Element", "Line", "Blank intensity", "Std intensity", "Blank/Std intensity ratio", "Detection limit = 2x est. blank (ppm)"]
    _sec(ws, SR, cd, "CONTAMINATION CHECK", len(cols_d), C_BRWN)
    for j, h in enumerate(cols_d): 
        _cell(ws, HR, cd+j, h, bg=C_MBRN, fc="FFFFFF", bold=True)
        ws.column_dimensions[get_column_letter(cd+j)].width = 16
    for i, (_, row) in enumerate(contam_df.iterrows()):
        r, bg = HR+1+i, (C_LBRN if i%2==0 else C_WHITE)
        for j, h in enumerate(cols_d): _cell(ws, r, cd+j, row.get(h, ""), bg=bg)

# ═══════════════════════════════════════════════════════════════════════════════
# 4. STREAMLIT APP RUNNER
# ═══════════════════════════════════════════════════════════════════════════════

st.title("🧪 ICP-OES Auto-Processor")
st.markdown("Upload your raw ICP Excel file to automatically identify blocks, clean data, compute robust replicate statistics, and rank emission lines.")

with st.expander("📖 Guide: How to format your input file", expanded=False):
    st.markdown("""
    **1. File & Sheet Naming**
    * **File Type:** Excel workbook (`.xlsx`).
    * **Sheet Names:** The app automatically scans all sheets. If it finds valid ICP data blocks on a sheet, it will process it regardless of the sheet's name.

    **2. Data Structure**
    * The sheet must contain a column named exactly **`Solution Label`**.
    * Element column headers must include the wavelength and **`nm`** (e.g., `Al 237.312 nm ppm`).

    **3. Naming Conventions (`Solution Label`)**
    * **QC Standards:** Must contain **`QC`**, followed by the type and target concentration (e.g., `QC MC 10`, `MC QC 10`, or `Si QC 5`).
    * **Blanks:** Must contain the word **`Blank`** (e.g., `Cal Blank`) or be exactly **`B`**. 
      *(Note: Experimental samples containing "DI", like `DI Water`, are treated as normal experimental samples and will be kept in your final overview).*
    * **Replicates:** To calculate Mean and SD automatically, use a consistent base name followed by a trailing letter or number (e.g., `Sample a`, `Sample b` or `Sample 1`, `Sample 2`).

    **4. Limits & Flagging**
    * The app detects `u` or `<` flags and treats them as `0.0`.
    * Mean calculations that fall below the dynamic method detection limit will be reported as `< [DL]`.
    """)

uploaded_file = st.file_uploader("Upload Excel File (.xlsx)", type=["xlsx"])

if uploaded_file is not None:
    if st.button("Process Data", type="primary"):
        with st.spinner("Processing ICP Data..."):
            wb = Workbook()
            sheets_count = 0
            log_messages = []
            
            try:
                xl = pd.ExcelFile(uploaded_file)
                
                target_sheets = [s for s in xl.sheet_names if "ICP" in s.upper()]
                if not target_sheets:
                    target_sheets = xl.sheet_names
                    
                for s_name in target_sheets:
                    log_messages.append(f"**Processing:** `{s_name}`...")
                    try:
                        raw = pd.read_excel(xl, sheet_name=s_name, header=None, dtype=str).fillna("")
                        
                        blocks = extract_all_data_blocks(raw)
                        conc_df, intens_df, unadj_df = identify_blocks(blocks)
                        
                        if conc_df.empty:
                            log_messages.append(f"  - ⏭️ *Skipped:* No valid Concentration block found in `{s_name}`.")
                            continue
                            
                        if intens_df.empty:
                            log_messages.append(f"  - ⚠️ *Warning:* No Intensity block detected. Skipping Contamination check.")
                        if unadj_df.empty:
                            log_messages.append(f"  - ⚠️ *Warning:* No Unadjusted block detected. Unadjusted columns will be blank.")
                        
                        cal_lines, _ = get_calibrated_lines(conc_df)
                        contam_df = step1_contamination(conc_df, intens_df, cal_lines)
                        qc_df = step2_qc_statistics(conc_df, unadj_df, cal_lines)
                        summ_df, best_df = step3_qc_summary(qc_df)
                        uni_map, elems, pub_df, rep_df = step4_data_overview(conc_df, unadj_df, best_df, contam_df, summ_df)
                        
                        if not pub_df.empty: 
                            write_sheet(wb, s_name, elems, pub_df, rep_df, summ_df, qc_df, contam_df, uni_map)
                            sheets_count += 1
                            log_messages.append(f"  - ✅ Successfully generated report for `{s_name}`.")
                    except Exception as e: 
                        log_messages.append(f"  - 💥 *CRITICAL ERROR* on `{s_name}`: {e}")
                        
                for msg in log_messages:
                    st.write(msg)
                    
                if sheets_count > 0:
                    if "Sheet" in wb.sheetnames: del wb["Sheet"]
                    output = io.BytesIO()
                    wb.save(output)
                    output.seek(0)
                    
                    st.success("🎉 Processing Complete! Download your formatted report below.")
                    st.download_button(
                        label="📥 Download Processed Excel Report",
                        data=output,
                        file_name=f"{uploaded_file.name.replace('.xlsx', '')}_Processed.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
                else:
                    st.error("No valid data found to process. Please check the Excel file formatting.")
                    
            except Exception as e:
                st.error(f"Failed to read the Excel file. Error: {e}")
