import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import io
import re
import textwrap

st.set_page_config(page_title="XRD Analysis App", layout="wide")

st.title("XRD Plotter & Quantitative Analysis")
st.markdown("Upload your raw XRD text file and the exported CSV database file to generate a combined spectra and donut chart.")

col1, col2 = st.columns(2)
with col1:
    txt_file = st.file_uploader("1. Upload Raw XRD Data (.txt)", type=['txt'])
with col2:
    csv_file = st.file_uploader("2. Upload Database Export (.csv)", type=['csv'])

def simplify_text(text):
    return re.sub(r'[^a-zA-Z0-9]', '', str(text).lower())

if txt_file and csv_file:
    st.success("Files loaded successfully.")
    
    # ---------------------------------------------------------
    # 1. PARSE RAW TXT
    # ---------------------------------------------------------
    raw_2theta = []
    raw_intensity = []
    lines = txt_file.getvalue().decode("utf-8", errors="ignore").splitlines()
    data_started = False
    for line in lines:
        if line.strip() == '[Data]':
            data_started = True
            continue
        if data_started:
            parts = line.strip().split(',')
            if len(parts) >= 2:
                try:
                    raw_2theta.append(float(parts[0]))
                    raw_intensity.append(float(parts[1]))
                except ValueError:
                    pass
    
    # ---------------------------------------------------------
    # 2. PARSE CSV (Peaks + Robust Quant Extraction)
    # ---------------------------------------------------------
    df = pd.read_csv(csv_file, header=None)
    
    header_row_idx = -1
    for idx, row in df.iterrows():
        row_strs = [str(x).strip().lower() for x in row.values]
        if "2theta (°)" in row_strs or "2theta" in "".join(row_strs):
            header_row_idx = idx
            break

    # Extract Peak Database
    peaks = {}
    if header_row_idx != -1:
        name_row = df.iloc[header_row_idx - 1]
        header_row = df.iloc[header_row_idx]
        
        current_comp = None
        for col_idx in range(len(header_row)):
            if pd.notna(name_row[col_idx]) and str(name_row[col_idx]).strip() != "":
                current_comp = str(name_row[col_idx]).strip()
            
            if "2theta" in str(header_row[col_idx]).lower():
                if col_idx + 1 < len(header_row) and "i fix" in str(header_row[col_idx + 1]).lower():
                    twotheta_vals = []
                    ifix_vals = []
                    for val_idx in range(header_row_idx + 1, len(df)):
                        tt = df.iloc[val_idx, col_idx]
                        ifix = df.iloc[val_idx, col_idx + 1]
                        if pd.isna(tt) or str(tt).strip() == "":
                            break
                        try:
                            twotheta_vals.append(float(tt))
                            ifix_vals.append(float(ifix))
                        except ValueError:
                            pass
                    if current_comp and len(twotheta_vals) > 0:
                        if current_comp not in peaks:
                            peaks[current_comp] = {'2theta': twotheta_vals, 'ifix': ifix_vals}
                        else:
                            peaks[current_comp]['2theta'].extend(twotheta_vals)
                            peaks[current_comp]['ifix'].extend(ifix_vals)

    # Dynamic Extraction for Quantitative Data
    pre_fill_data = []
    cryst_val = ""
    
    comp_col, pdf_col, sq_col, cryst_col = 1, 3, 4, 6
    found_header = False
    
    for idx in range(max(1, header_row_idx if header_row_idx != -1 else len(df))):
        row_strs = [str(x).strip().lower() for x in df.iloc[idx].values]
        if "compound name" in row_strs:
            found_header = True
            comp_col = row_strs.index("compound name")
            if "pdf name" in row_strs: pdf_col = row_strs.index("pdf name")
            if "s-q %" in row_strs: sq_col = row_strs.index("s-q %")
            if "crystallinity" in row_strs: cryst_col = row_strs.index("crystallinity")
            
            for j in range(idx + 1, header_row_idx if header_row_idx != -1 else len(df)):
                data_row = df.iloc[j]
                
                if "2theta" in "".join([str(x).lower() for x in data_row.values]):
                    break
                    
                comp = str(data_row[comp_col]).strip() if comp_col < len(data_row) else ""
                if comp and comp.lower() not in ["nan", "none", ""]:
                    pdf = str(data_row[pdf_col]).strip() if pdf_col < len(data_row) else ""
                    if pdf.lower() == "nan": pdf = ""
                    
                    sq_str = str(data_row[sq_col]).strip() if sq_col < len(data_row) else ""
                    try:
                        sq_float = float(sq_str)
                        sq_val = sq_float * 100 if sq_float <= 1.0 else sq_float
                    except ValueError:
                        sq_val = 0.0
                        
                    pre_fill_data.append({
                        "Compound Name": comp,
                        "PDF Name": pdf,
                        "S-Q %": sq_val
                    })
                    
                c_str = str(data_row[cryst_col]).strip() if cryst_col < len(data_row) else ""
                if c_str and c_str.lower() not in ["nan", "none", ""]:
                    try:
                        c_float = float(c_str)
                        cryst_val = str(round(c_float * 100 if c_float <= 1.0 else c_float, 2))
                    except ValueError:
                        pass
            break

    if not pre_fill_data:
        pre_fill_data = [{"Compound Name": "", "PDF Name": "", "S-Q %": 0.0}]

    df_quant = pd.DataFrame(pre_fill_data)

    # ---------------------------------------------------------
    # 3. INTERACTIVE UI: DATA EDITOR
    # ---------------------------------------------------------
    st.markdown("---")
    st.subheader("Quantitative Data Editor")
    st.markdown("Verify the extracted data below. **You can edit cells, delete rows, or add new rows directly in this table.**")
    
    edited_df = st.data_editor(
        df_quant, 
        num_rows="dynamic", 
        use_container_width=True,
        column_config={
            "S-Q %": st.column_config.NumberColumn("S-Q %", min_value=0.0, max_value=100.0, format="%.2f")
        }
    )
    
    col_c1, col_c2 = st.columns([1, 3])
    with col_c1:
        include_cryst = st.checkbox("Show Crystallinity in Plot", value=(cryst_val != ""))
    with col_c2:
        if include_cryst:
            cryst_input = st.text_input("Crystallinity (%)", value=cryst_val)
        else:
            cryst_input = ""

    st.sidebar.header("Plot Settings")
    stick_scale = st.sidebar.slider("Reference Stick Height", min_value=0.1, max_value=1.0, value=0.4, step=0.05)
    base_font = st.sidebar.slider("Text Size", min_value=10, max_value=30, value=18, step=1)
    
    title_font = base_font + 4
    cryst_font = base_font + 2

    compounds_info = []
    for idx, row in edited_df.iterrows():
        name = str(row["Compound Name"]).strip()
        if name and name.lower() != "nan" and row["S-Q %"] > 0:
            compounds_info.append({
                'name': name,
                'pdf': str(row["PDF Name"]).strip() if pd.notna(row["PDF Name"]) and str(row["PDF Name"]).lower() not in ["nan", "none"] else "",
                'sq': float(row["S-Q %"])
            })

    if st.button("Generate Figure", type="primary"):
        plt.rcParams['font.family'] = 'Arial'
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), gridspec_kw={'width_ratios': [2.5, 1]})
        colors = ['#d62728', '#1f77b4', '#2ca02c', '#9467bd', '#17becf', '#ff7f0e']
        
        if len(raw_intensity) > 0:
            max_raw = max(raw_intensity)
            min_raw = min(raw_intensity)
        else:
            max_raw, min_raw = 100, 0
            
        max_ifix = max([max(data['ifix']) for data in peaks.values() if len(data['ifix']) > 0] + [1])
        stick_max_height = (max_raw - min_raw) * stick_scale if max_raw > min_raw else 50
        offset = stick_max_height * 1.1 - min_raw if max_raw > min_raw else 5

        ax1.plot(raw_2theta, [val + offset for val in raw_intensity], color='black', label='Raw Data', linewidth=1.5)

        for i, (comp, data) in enumerate(peaks.items()):
            c = colors[i % len(colors)]
            x = data['2theta']
            y = [val / max_ifix * stick_max_height for val in data['ifix']]
            
            pdf_str = ""
            comp_simple = simplify_text(comp)
            for ci in compounds_info:
                ci_name_simple = simplify_text(ci['name'])
                if len(ci_name_simple) > 2 and (ci_name_simple in comp_simple or comp_simple in ci_name_simple):
                    if ci['pdf']:
                        pdf_str = f"\n{ci['pdf']}"
                    break
            
            # --- FIXED LEGEND TEXT ---
            # Wrap long compound names (stops at 30 characters and drops to next line)
            wrapped_comp = "\n".join(textwrap.wrap(comp, width=30))
            label_str = f"{wrapped_comp}{pdf_str}"
            
            ax1.plot([], [], color=c, label=label_str, linewidth=3)
            ax1.vlines(x, ymin=0, ymax=y, color=c, linewidth=2.0)

        ax1.set_xlabel('2$\\theta$ Angle (degrees)', fontsize=title_font, fontweight='bold')
        ax1.set_ylabel('Normalized Intensity', fontsize=title_font, fontweight='bold')
        ax1.set_xlim([min(raw_2theta) if raw_2theta else 5, max(raw_2theta) if raw_2theta else 80])
        
        # --- FIXED HEADROOM ---
        # Add 25% empty space at the top of the plot to prevent legend overlap
        highest_plot_point = max_raw + offset
        ax1.set_ylim(bottom=0, top=highest_plot_point * 1.25)
        
        ax1.set_yticks([]) 
        ax1.tick_params(axis='x', labelsize=base_font)
        
        # --- FIXED LEGEND PADDING ---
        legend_font = max(10, base_font - 2) 
        # labelspacing=1.5 adds vertical breathing room between multiline items
        ax1.legend(fontsize=legend_font, frameon=False, loc='upper right', labelspacing=1.5)
        
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        if len(compounds_info) > 0:
            labels = []
            sizes = []
            for comp in compounds_info:
                # Keep pie chart labels reasonably wrapped too
                name_part = comp['name'].replace(" (", "\n(")
                name_part = "\n".join(textwrap.wrap(name_part, width=25))
                labels.append(name_part)
                sizes.append(comp['sq'])
                
            colors_donut = colors[:len(labels)]
            
            wedges, texts, autotexts = ax2.pie(
                sizes, labels=labels, autopct='%1.1f%%', 
                startangle=90, colors=colors_donut, 
                textprops={'fontsize': base_font, 'fontweight': 'bold'},
                labeldistance=1.1,
                pctdistance=0.75
            )

            centre_circle = plt.Circle((0,0), 0.50, fc='white')
            ax2.add_artist(centre_circle)

            if include_cryst and cryst_input.strip() != "":
                c_text = cryst_input.strip()
                if not c_text.endswith("%"):
                    c_text += "%"
                ax2.text(0, 0, f"Crystallinity:\n{c_text}", ha='center', va='center', fontsize=cryst_font, fontweight='bold')

            ax2.axis('equal')
        else:
            ax2.text(0.5, 0.5, "No Quantitative Data Provided", ha='center', va='center', fontsize=base_font)
            ax2.axis('off')

        plt.tight_layout()
        st.pyplot(fig)
        
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=300, bbox_inches='tight')
        st.download_button(
            label="Download High-Res PNG",
            data=buf.getvalue(),
            file_name="XRD_Combined_Analysis.png",
            mime="image/png"
        )
