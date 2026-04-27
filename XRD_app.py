import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import io

# Setup the browser tab title and page width
st.set_page_config(page_title="XRD Analysis App", layout="wide")

st.title("XRD Plotter & Quantitative Analysis")
st.markdown("Upload your raw XRD text file and the exported CSV database file to generate a combined spectra and donut chart.")

# Create two columns for the drag-and-drop file uploaders
col1, col2 = st.columns(2)

with col1:
    txt_file = st.file_uploader("1. Upload Raw XRD Data (.txt)", type=['txt'])

with col2:
    csv_file = st.file_uploader("2. Upload Database Export (.csv)", type=['csv'])

# Only run the code once both files have been uploaded by the user
if txt_file and csv_file:
    st.success("Files loaded successfully. Generating plot...")
    
    # 1. Parse TXT (reading from the uploaded file stream)
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
                    
    # 2. Parse CSV (reading from the uploaded file stream)
    df = pd.read_csv(csv_file, header=None)
    
    compounds_info = []
    crystallinity = None
    
    for idx, row in df.iterrows():
        if str(row[0]).strip() == "Sample":
            continue
        if pd.notna(row[1]) and "Compound Name" not in str(row[1]):
            comp = str(row[1]).strip()
            pdf = str(row[3]).strip() if pd.notna(row[3]) else ""
            sq = row[4]
            if pd.notna(row[6]) and idx == 1:
                crystallinity = row[6]
            if pd.notna(sq):
                try:
                    compounds_info.append({'name': comp, 'pdf': pdf, 'sq': float(sq)})
                except:
                    pass
        if idx > 15:
            break
            
    header_row_idx = -1
    for idx, row in df.iterrows():
        if "2Theta (°)" in str(row.values) or "2Theta" in str(row.values):
            header_row_idx = idx
            break
            
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
                        except:
                            pass
                    if current_comp and len(twotheta_vals) > 0:
                        if current_comp not in peaks:
                            peaks[current_comp] = {'2theta': twotheta_vals, 'ifix': ifix_vals}
                        else:
                            peaks[current_comp]['2theta'].extend(twotheta_vals)
                            peaks[current_comp]['ifix'].extend(ifix_vals)

    # Adding a Sidebar Slider so the user can tweak the stick height dynamically
    st.sidebar.header("Plot Settings")
    stick_scale = st.sidebar.slider("Reference Stick Height (Relative to raw data)", min_value=0.1, max_value=1.0, value=0.4, step=0.05)
    
    # 3. Create Figure
    plt.rcParams['font.family'] = 'Arial'
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), gridspec_kw={'width_ratios': [2.5, 1]})
    colors = ['#d62728', '#1f77b4', '#2ca02c', '#9467bd', '#17becf', '#ff7f0e']
    
    # ---------------- AXIS 1: MAIN XRD SPECTRA ----------------
    if len(raw_intensity) > 0:
        max_raw = max(raw_intensity)
        min_raw = min(raw_intensity)
    else:
        max_raw, min_raw = 100, 0
        
    max_ifix = max([max(data['ifix']) for data in peaks.values() if len(data['ifix']) > 0] + [1])
    
    # Use the variable from the Streamlit slider to adjust the stick height
    stick_max_height = (max_raw - min_raw) * stick_scale if max_raw > min_raw else 50
    offset = stick_max_height * 1.1 - min_raw if max_raw > min_raw else 5

    ax1.plot(raw_2theta, [val + offset for val in raw_intensity], color='black', label='Raw Data', linewidth=1.5)

    for i, (comp, data) in enumerate(peaks.items()):
        c = colors[i % len(colors)]
        x = data['2theta']
        y = [val / max_ifix * stick_max_height for val in data['ifix']]
        
        pdf_str = ""
        for ci in compounds_info:
            if ci['name'] == comp and ci['pdf']:
                pdf_str = f" - {ci['pdf']}"
                break
        
        label_str = f"{comp}{pdf_str}"
        ax1.plot([], [], color=c, label=label_str, linewidth=3)
        ax1.vlines(x, ymin=0, ymax=y, color=c, linewidth=2.0)

    ax1.set_xlabel('2$\\theta$ Angle (degrees)', fontsize=22, fontweight='bold')
    ax1.set_ylabel('Normalized Intensity', fontsize=22, fontweight='bold')
    ax1.set_xlim([min(raw_2theta) if raw_2theta else 5, max(raw_2theta) if raw_2theta else 80])
    ax1.set_ylim(bottom=0)
    ax1.set_yticks([]) 
    ax1.tick_params(axis='x', labelsize=18)
    ax1.legend(fontsize=18, frameon=False, loc='upper right')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    # ---------------- AXIS 2: QUANTITATIVE DONUT CHART ----------------
    if len(compounds_info) > 0:
        labels = []
        for comp in compounds_info:
            name_part = comp['name'].replace(" (", "\n(")
            labels.append(name_part)
            
        sizes = [comp['sq'] * 100 if comp['sq'] <= 1.0 else comp['sq'] for comp in compounds_info] 
        colors_donut = colors[:len(labels)]
        
        wedges, texts, autotexts = ax2.pie(
            sizes, labels=labels, autopct='%1.1f%%', 
            startangle=90, colors=colors_donut, 
            textprops={'fontsize': 18, 'fontweight': 'bold'},
            labeldistance=1.1,
            pctdistance=0.75
        )

        centre_circle = plt.Circle((0,0), 0.50, fc='white')
        ax2.add_artist(centre_circle)

        if crystallinity is not None:
            try:
                cryst_val = float(crystallinity) * 100 if float(crystallinity) <= 1.0 else float(crystallinity)
                cryst_text = f"Crystallinity:\n{cryst_val:.1f}%"
            except:
                cryst_text = f"Crystallinity:\n{crystallinity}"
            ax2.text(0, 0, cryst_text, ha='center', va='center', fontsize=20, fontweight='bold')

        ax2.axis('equal')
    else:
        ax2.axis('off')

    plt.tight_layout()
    
    # Render the matplotlib figure directly inside the Streamlit web browser
    st.pyplot(fig)
    
    # Render a download button underneath the plot to save the high-res PNG
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=300, bbox_inches='tight')
    st.download_button(
        label="Download High-Res PNG",
        data=buf.getvalue(),
        file_name="XRD_Combined_Analysis.png",
        mime="image/png"
    )