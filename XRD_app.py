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

    # Extract Quantitative Data
    df_quant_raw = []
    cryst_val = ""
    comp_col, pdf_col, sq_col, cryst_col = 1, 3, 4, 6
    
    for idx in range(max(1, header_row_idx if header_row_idx != -1 else len(df))):
        row_strs = [str(x).strip().lower() for x in df.iloc[idx].values]
        if "compound name" in row_strs:
            if "compound name" in row_strs: comp_col = row_strs.index("compound name")
            if "pdf name" in row_strs: pdf_col = row_strs.index("pdf name")
            if "s-q %" in row_strs: sq_col = row_strs.index("s-q %")
            if "crystallinity" in row_strs: cryst_col = row_strs.index("crystallinity")
            
            for j in range(idx + 1, header_row_idx if header_row_idx != -1 else len(df)):
                data_row = df.iloc[j]
                if "2theta" in "".join([str(x).lower() for x in data_row.values]): break
                    
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
                        
                    df_quant_raw.append({"name": comp, "pdf": pdf, "sq": sq_val})
                    
                c_str = str(data_row[cryst_col]).strip() if cryst_col < len(data_row) else ""
                if c_str and c_str.lower() not in ["nan", "none", ""]:
                    try:
                        c_float = float(c_str)
                        cryst_val = str(round(c_float * 100 if c_float <= 1.0 else c_float, 2))
                    except ValueError:
                        pass
            break

    # Link the Peaks directly to the Data Table
    pre_fill_data = []
    assigned_q = set()
    
    for p_key in peaks.keys():
        best_name = p_key
        best_pdf = ""
        best_sq = 0.0
        p_key_simp = simplify_text(p_key)
        
        for q in df_quant_raw:
            q_name_simp = simplify_text(q["name"])
            if len(q_name_simp) > 2 and (q_name_simp in p_key_simp or p_key_simp in q_name_simp):
                best_name = q["name"]
                best_pdf = q["pdf"]
                best_sq = q["sq"]
                assigned_q.add(q["name"])
                break
                
        pre_fill_data.append({
            "Peak Mapping (DO NOT EDIT)": p_key,
            "Compound Name": best_name,
            "PDF Name": best_pdf,
            "S-Q %": best_sq
        })

    # Add any extra quantitative data (like Amorphous) that has no peak sticks
    for q in df_quant_raw:
        if q["name"] not in assigned_q:
            pre_fill_data.append({
                "Peak Mapping (DO NOT EDIT)": "No Peak Data",
                "Compound Name": q["name"],
                "PDF Name": q["pdf"],
                "S-Q %": q["sq"]
            })

    if not pre_fill_data:
        pre_fill_data = [{"Peak Mapping (DO NOT EDIT)": "", "Compound Name": "", "PDF Name": "", "S-Q %": 0.0}]

    df_quant = pd.DataFrame(pre_fill_data)

    # ---------------------------------------------------------
    # 3. INTERACTIVE UI: DATA EDITOR
    # ---------------------------------------------------------
    st.markdown("---")
    st.subheader("Quantitative Data Editor")
    st.markdown("Verify the extracted data below. **The plot reads directly from this table.**")
    
    edited_df = st.data_editor(
        df_quant, 
        num_rows="dynamic", 
        use_container_width=True,
        column_config={
            "Peak Mapping (DO NOT EDIT)": st.column_config.TextColumn(disabled=True),
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
    fig_height = st.sidebar.slider("Figure Height (Stretch Y-Axis)", min_value=6, max_value=16, value=10, step=1)
    stick_scale = st.sidebar.slider("Reference Stick Height", min_value=0.1, max_value=1.0, value=0.4, step=0.05)
    raw_lw = st.sidebar.slider("Raw Data Line Thickness", min_value=0.2, max_value=2.0, value=0.8, step=0.1)
    base_font = st.sidebar.slider("Text Size", min_value=10, max_value=30, value=18, step=1)
    
    title_font = base_font + 4
    cryst_font = base_font + 2

    # ---------------------------------------------------------
    # 5. GENERATE PLOT
    # ---------------------------------------------------------
    if st.button("Generate Figure", type="primary"):
        plt.rcParams['font.family'] = 'Arial'
        
        # INCREASED HEIGHT to stretch the aspect ratio 
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, fig_height), gridspec_kw={'width_ratios': [2.5, 1]})
        colors = ['#d62728', '#1f77b4', '#2ca02c', '#9467bd', '#17becf', '#ff7f0e', '#8c564b', '#e377c2']
        
        if len(raw_intensity) > 0:
            max_raw = max(raw_intensity)
            min_raw = min(raw_intensity)
        else:
            max_raw, min_raw = 100, 0
            
        max_ifix = max([max(data['ifix']) for data in peaks.values() if len(data['ifix']) > 0] + [1])
        stick_max_height = (max_raw - min_raw) * stick_scale if max_raw > min_raw else 50
        offset = stick_max_height * 1.1 - min_raw if max_raw > min_raw else 5

        # Plot Raw Data with adjustable line width
        ax1.plot(raw_2theta, [val + offset for val in raw_intensity], color='black', label='Raw Data', linewidth=raw_lw)

        legend_handles = []
        legend_labels = []
        raw_handle = plt.Line2D([0], [0], color='black', linewidth=raw_lw)
        legend_handles.append(raw_handle)
        legend_labels.append('Raw Data')

        donut_labels = []
        donut_sizes = []
        donut_colors = []

        color_idx = 0

        # Loop purely over the Data Editor Table to build the plot
        for idx, row in edited_df.iterrows():
            comp_name = str(row["Compound Name"]).strip()
            sq = float(row["S-Q %"])
            if not comp_name or comp_name.lower() in ["nan", "none"]: continue
            
            c = colors[color_idx % len(colors)]
            pdf = str(row["PDF Name"]).strip()
            
            # One line string for the legend
            pdf_str = f" - {pdf}" if pdf and pdf.lower() not in ["nan", "none"] else ""
            label_str = f"{comp_name}{pdf_str}"
            
            p_id = row["Peak Mapping (DO NOT EDIT)"]
            
            # If the row has peak data, draw the sticks and add to legend
            if p_id in peaks:
                data = peaks[p_id]
                x = data['2theta']
                y = [val / max_ifix * stick_max_height for val in data['ifix']]
                
                proxy = plt.Line2D([0], [0], color=c, linewidth=3)
                legend_handles.append(proxy)
                legend_labels.append(label_str)
                
                ax1.vlines(x, ymin=0, ymax=y, color=c, linewidth=2.0)
                
                # Only iterate the color if we actually drew peak data
                color_idx += 1
            else:
                # E.g. Amorphous content. Use grey for the donut.
                c = "#aaaaaa"

            # Add valid data to Donut
            if sq > 0:
                # Wrap Pie chart names so they don't break the circle
                wrapped_name = "\n".join(textwrap.wrap(comp_name, width=20))
                donut_labels.append(wrapped_name)
                donut_sizes.append(sq)
                donut_colors.append(c)

        ax1.set_xlabel('2$\\theta$ Angle (degrees)', fontsize=title_font, fontweight='bold')
        ax1.set_ylabel('Normalized Intensity', fontsize=title_font, fontweight='bold')
        ax1.set_xlim([min(raw_2theta) if raw_2theta else 5, max(raw_2theta) if raw_2theta else 80])
        
        highest_plot_point = max_raw + offset
        ax1.set_ylim(bottom=0, top=highest_plot_point * 1.15)
        
        ax1.set_yticks([]) 
        ax1.tick_params(axis='x', labelsize=base_font)
        
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        # ALIGN LEGEND TO BOTTOM X-AXIS, UNDERNEATH DONUT CHART
        legend_font = max(10, base_font - 2) 
        ax1.legend(
            legend_handles, legend_labels, 
            fontsize=legend_font, 
            frameon=False, 
            loc='lower left', 
            bbox_to_anchor=(1.05, 0.0) # 1.05 puts it right of ax1. 0.0 aligns it perfectly with bottom X-axis.
        )

        # AXIS 2: DONUT CHART
        if len(donut_sizes) > 0:
            wedges, texts, autotexts = ax2.pie(
                donut_sizes, labels=donut_labels, autopct='%1.1f%%', 
                startangle=90, colors=donut_colors, 
                textprops={'fontsize': base_font, 'fontweight': 'bold'},
                labeldistance=1.1,
                pctdistance=0.75
            )

            centre_circle = plt.Circle((0,0), 0.50, fc='white')
            ax2.add_artist(centre_circle)

            if include_cryst and cryst_input.strip() != "":
                c_text = cryst_input.strip()
                if not c_text.endswith("%"): c_text += "%"
                ax2.text(0, 0, f"Crystallinity:\n{c_text}", ha='center', va='center', fontsize=cryst_font, fontweight='bold')

            ax2.axis('equal')
        else:
            ax2.axis('off')

        plt.subplots_adjust(right=0.9)
        
        st.pyplot(fig)
        
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=300, bbox_inches='tight')
        st.download_button(
            label="Download High-Res PNG",
            data=buf.getvalue(),
            file_name="XRD_Combined_Analysis.png",
            mime="image/png"
        )
