import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import io

st.set_page_config(page_title="XRD Analysis App", layout="wide")

st.title("XRD Plotter & Quantitative Analysis")
st.markdown("Upload your raw XRD text file and the exported CSV database file to generate a combined spectra and donut chart.")

col1, col2 = st.columns(2)
with col1:
    txt_file = st.file_uploader("1. Upload Raw XRD Data (.txt)", type=['txt'])
with col2:
    csv_file = st.file_uploader("2. Upload Database Export (.csv)", type=['csv'])

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
    # 2. PARSE CSV (Peaks + Attempted Quant Extraction)
    # ---------------------------------------------------------
    df = pd.read_csv(csv_file, header=None)
    
    # Find where the Peak data starts
    header_row_idx = -1
    for idx, row in df.iterrows():
        if "2Theta (°)" in str(row.values) or "2Theta" in str(row.values):
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
            
            if "2Theta" in str(header_row[col_idx]):
                if col_idx + 1 < len(header_row) and "I fix" in str(header_row[col_idx + 1]):
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

    # Attempt to pre-fill Quantitative Data
    pre_fill_data = []
    cryst_val = ""
    
    if header_row_idx > 0:
        for idx in range(header_row_idx):
            row = df.iloc[idx]
            
            # Look for crystallinity
            if idx == 1 and len(row) > 6:
                c = str(row[6]).strip()
                if c.lower() not in ["nan", "crystallinity", "none", ""]:
                    try:
                        c_float = float(c)
                        # Convert fraction to percentage if necessary
                        cryst_val = str(round(c_float * 100 if c_float <= 1.0 else c_float, 2))
                    except ValueError:
                        pass
                        
            # Look for Compound Data
            if len(row) > 4:
                comp = str(row[1]).strip()
                if comp and comp.lower() not in ["nan", "compound name", "none"]:
                    sample = str(row[0]).strip() if pd.notna(row[0]) else ""
                    formula = str(row[2]).strip() if pd.notna(row[2]) else ""
                    pdf = str(row[3]).strip() if pd.notna(row[3]) else ""
                    if pdf.lower() == "nan": pdf = ""
                    
                    try:
                        sq_float = float(row[4])
                        if not np.isnan(sq_float):
                            # Ensure percentage is scaled to 100 (not 1.0)
                            sq_val = sq_float * 100 if sq_float <= 1.0 else sq_float
                            pre_fill_data.append({
                                "Sample": sample if sample.lower() != "nan" else "",
                                "Compound Name": comp,
                                "Formula": formula if formula.lower() != "nan" else "",
                                "PDF Name": pdf,
                                "S-Q %": sq_val
                            })
                    except (ValueError, TypeError):
                        pass

    # If parsing failed completely, provide an empty structure for the user
    if not pre_fill_data:
        pre_fill_data = [{"Sample": "", "Compound Name": "", "Formula": "", "PDF Name": "", "S-Q %": 0.0}]

    df_quant = pd.DataFrame(pre_fill_data)

    # ---------------------------------------------------------
    # 3. INTERACTIVE UI: DATA EDITOR
    # ---------------------------------------------------------
    st.markdown("---")
    st.subheader("Quantitative Data Editor")
    st.markdown("Verify the extracted data below. **You can edit cells, delete rows, or add new rows directly in this table.**")
    
    # The editable table
    edited_df = st.data_editor(
        df_quant, 
        num_rows="dynamic", 
        use_container_width=True,
        column_config={
            "S-Q %": st.column_config.NumberColumn("S-Q %", min_value=0.0, max_value=100.0, format="%.2f")
        }
    )
    
    # Standalone Crystallinity Box
    cryst_input = st.text_input("Crystallinity (%)", value=cryst_val, help="Leave blank if not applicable.")

    # ---------------------------------------------------------
    # 4. PLOT SETTINGS SLIDERS
    # ---------------------------------------------------------
    st.sidebar.header("Plot Settings")
    stick_scale = st.sidebar.slider("Reference Stick Height", min_value=0.1, max_value=1.0, value=0.4, step=0.05)
    base_font = st.sidebar.slider("Text Size", min_value=10, max_value=30, value=18, step=1)
    
    title_font = base_font + 4
    cryst_font = base_font + 2

    # Parse edited DataFrame back into our plotting list
    compounds_info = []
    for idx, row in edited_df.iterrows():
        name = str(row["Compound Name"]).strip()
        if name and name.lower() != "nan" and row["S-Q %"] > 0:
            compounds_info.append({
                'name': name,
                'pdf': str(row["PDF Name"]).strip() if pd.notna(row["PDF Name"]) and str(row["PDF Name"]).lower() != "nan" else "",
                'sq': float(row["S-Q %"])
            })

    # ---------------------------------------------------------
    # 5. GENERATE PLOT
    # ---------------------------------------------------------
    if st.button("Generate Figure", type="primary"):
        plt.rcParams['font.family'] = 'Arial'
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), gridspec_kw={'width_ratios': [2.5, 1]})
        colors = ['#d62728', '#1f77b4', '#2ca02c', '#9467bd', '#17becf', '#ff7f0e']
        
        # AXIS 1: MAIN SPECTRA
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
            for ci in compounds_info:
                if ci['name'].lower() == comp.lower() and ci['pdf']:
                    pdf_str = f" - {ci['pdf']}"
                    break
            
            label_str = f"{comp}{pdf_str}"
            ax1.plot([], [], color=c, label=label_str, linewidth=3)
            ax1.vlines(x, ymin=0, ymax=y, color=c, linewidth=2.0)

        ax1.set_xlabel('2$\\theta$ Angle (degrees)', fontsize=title_font, fontweight='bold')
        ax1.set_ylabel('Normalized Intensity', fontsize=title_font, fontweight='bold')
        ax1.set_xlim([min(raw_2theta) if raw_2theta else 5, max(raw_2theta) if raw_2theta else 80])
        ax1.set_ylim(bottom=0)
        ax1.set_yticks([]) 
        ax1.tick_params(axis='x', labelsize=base_font)
        ax1.legend(fontsize=base_font, frameon=False, loc='upper right')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        # AXIS 2: DONUT CHART
        if len(compounds_info) > 0:
            labels = []
            sizes = []
            for comp in compounds_info:
                name_part = comp['name'].replace(" (", "\n(")
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

            if cryst_input.strip() != "":
                # Add % sign if the user just typed a number
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
