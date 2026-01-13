import sys
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output
import matplotlib.colors as mcolors

def darken_color(hex_color, factor=0.6):
    """Darken a hex color by the given factor (not used in final version)."""
    rgb = np.array(mcolors.to_rgb(hex_color))
    dark_rgb = np.clip(rgb * factor, 0, 1)
    return mcolors.to_hex(dark_rgb)

if len(sys.argv) != 2:
    print("Usage: python mtdna_viewer.py <tsv_file>")
    sys.exit(1)

tsv_file = sys.argv[1]
print(f"Loading {tsv_file}...")

# Load TSV data (expected format: chrM, Position, Ref, Alt, sample1, sample2, ...)
data = pd.read_csv(tsv_file, sep='\t')
MT_LENGTH = 16569  # Human mtDNA reference length

print(f"Columns ({len(data.columns)}): {list(data.columns)}")
print(f"Shape: {data.shape}")

# Column indices: chrM(0), Position(1), Ref(2), Alt(3), samples(4+)
pos_idx, ref_idx, alt_idx, sample_start = 1, 2, 3, 4

positions = data.iloc[:, pos_idx].astype(int).values
refs = data.iloc[:, ref_idx].fillna('N/A').astype(str).values
alts = data.iloc[:, alt_idx].fillna('N/A').astype(str).values

samples = data.columns[sample_start:]
NUM_SAMPLES = len(samples)
print(f"Samples ({NUM_SAMPLES}): {samples[0]} ... {samples[-1]}")

def safe_het_convert(x):
    """Safely convert heteroplasmy values to fractions (0.0-1.0)."""
    if pd.isna(x): return 0.0
    x_str = str(x).strip()
    if '%' in x_str:
        try: return float(x_str.replace('%', '')) / 100
        except: return 0.0
    try:
        val = float(x_str)
        return val/100 if val > 1 else val
    except: return 0.0

# Heteroplasmy matrix and angular positions
het_matrix = data.iloc[:, sample_start:].map(safe_het_convert).fillna(0).values
thetas = positions / MT_LENGTH * 2 * np.pi

# Human mtDNA genes with coordinates (standard rCRS annotation)
genes = [
    ("D-loop", 0, 576), ("D-loop", 16024, 16569),
    ("MT-RNR1", 648, 1601), ("MT-RNR2", 1671, 3229),
    ("MT-ND1", 3307, 4262), ("MT-ND2", 4470, 5511),
    ("MT-CO1", 5904, 7445), ("MT-CO2", 7586, 8269),
    ("MT-ATP8", 8366, 8572), ("MT-ATP6", 8527, 9207),
    ("MT-CO3", 9207, 9990), ("MT-ND3", 10059, 10404),
    ("MT-ND4L", 10470, 10766), ("MT-ND4", 10760, 12137),
    ("MT-ND5", 12337, 14148), ("MT-ND6", 14149, 14673),
    ("MT-CYB", 14747, 15887),
]

def position_to_gene(pos):
    """Map mtDNA position to gene name."""
    for g, s, e in genes:
        if s <= e and s <= pos <= e: return g
        if s > e and (pos >= s or pos <= e): return g
    return "Intergenic"

# Distinct colors for each gene
gene_colors = {
    "D-loop":"#FF4444", "MT-RNR1":"#44FF44", "MT-RNR2":"#4444FF", "MT-ND1":"#FFAA00",
    "MT-ND2":"#AA00FF", "MT-CO1":"#00FFFF", "MT-CO2":"#FF44AA", "MT-ATP8":"#44AAFF",
    "MT-ATP6":"#AAFF44", "MT-CO3":"#FF8844", "MT-ND3":"#44FF88", "MT-ND4L":"#8844FF",
    "MT-ND4":"#FF4488", "MT-ND5":"#4488FF", "MT-ND6":"#88FF44", "MT-CYB":"#44FFAA",
    "Intergenic":"#AAAAAA"
}

# Variant statistics (>1% threshold)
threshold = 0.01
variant_count = (het_matrix > threshold).sum(axis=1)
max_het = np.max(het_matrix, axis=1)
r_scale = np.clip(max_het / np.max(max_het), 0.0, 1.0)
total_sites = int((variant_count > 0).sum())
print(f"Variant sites (>1%): {total_sites}")

# 0% reference radius (markers start exactly here)
ZERO_RADIUS = 0.20

def create_figure(selected_gene=None):
    """Create 3D circular mtDNA heteroplasmy visualization."""
    fig = go.Figure()

    # Filter data by selected gene
    if selected_gene and selected_gene != "":
        mask = np.array([position_to_gene(p) == selected_gene for p in positions])
    else:
        mask = np.ones(len(positions), dtype=bool)

    disp_idx = np.where(mask)[0]
    disp_pos = positions[disp_idx]
    disp_theta = thetas[disp_idx]
    disp_z = variant_count[disp_idx]
    disp_r = r_scale[disp_idx]
    disp_gene = np.array([position_to_gene(p) for p in disp_pos])
    disp_ref = refs[disp_idx]
    disp_alt = alts[disp_idx]

    # Plot variant markers (gene-color coded, radial position = heteroplasmy)
    MARKER_SIZE = 8
    for i, idx in enumerate(disp_idx):
        theta = disp_theta[i]
        z_peak = disp_z[i]
        r_sc = disp_r[i]
        gene = disp_gene[i]
        
        # Radial position: 0% circle → 100% (linear scale)
        r_offset = ZERO_RADIUS
        r_pos = r_offset + (r_sc * (1.10 - r_offset))
        
        marker_color = gene_colors.get(gene, "#CCCCCC")
        
        fig.add_trace(go.Scatter3d(
            x=[r_pos * np.cos(theta)],
            y=[r_pos * np.sin(theta)],
            z=[z_peak],
            mode='markers',
            marker=dict(size=MARKER_SIZE, color=marker_color),
            hovertemplate=(
                f'<b>{gene}</b><br>'
                f'Position: <b>{int(disp_pos[i])}</b><br>'
                f'REF: {disp_ref[i]} → ALT: {disp_alt[i]}<br>'
                f'<b>Max heteroplasmy: {max_het[idx]:.1%}</b><br>'
                f'<b>Samples: {int(z_peak)}/{NUM_SAMPLES}</b><extra></extra>'
            ),
            showlegend=False,
            name=gene
        ))

    # Reference circles
    n_ring = 200
    full_thetas = np.linspace(0, 2 * np.pi, n_ring)
    
    # 0% reference circle (black, bold)
    fig.add_trace(go.Scatter3d(
        x=ZERO_RADIUS * np.cos(full_thetas),
        y=ZERO_RADIUS * np.sin(full_thetas),
        z=np.zeros(n_ring),
        mode='lines',
        line=dict(color='black', width=10),
        showlegend=False
    ))
    
    # 0% white background (anti-overlap)
    fig.add_trace(go.Scatter3d(
        x=(ZERO_RADIUS + 0.02) * np.cos(full_thetas),
        y=(ZERO_RADIUS + 0.02) * np.sin(full_thetas),
        z=np.zeros(n_ring),
        mode='lines',
        line=dict(color='white', width=18),
        showlegend=False
    ))
    
    # 100% outer circle
    fig.add_trace(go.Scatter3d(
        x=1.10 * np.cos(full_thetas),
        y=1.10 * np.sin(full_thetas),
        z=np.zeros(n_ring),
        mode='lines',
        line=dict(color='#999999', width=5),
        showlegend=False
    ))

    # Radial scale labels
    fig.add_trace(go.Scatter3d(
        x=[ZERO_RADIUS], y=[0], z=[0],
        mode='text', text=['0%'],
        textfont=dict(size=18, color='black', family='Arial Black'),
        showlegend=False, hoverinfo='skip'
    ))
    
    fig.add_trace(go.Scatter3d(
        x=[1.15], y=[0], z=[0],
        mode='text', text=['100%'],
        textfont=dict(size=18, color='#666666', family='Arial'),
        showlegend=False, hoverinfo='skip'
    ))

    # Gene region surfaces
    r_inner_all = ZERO_RADIUS + 0.05
    fig.update_layout(scene=dict(bgcolor='white'))
    
    if selected_gene and selected_gene != "":
        # Highlight selected gene
        for gname, s, e in genes:
            if gname == selected_gene:
                color = gene_colors[gname]
                npts = 80
                t1 = s / MT_LENGTH * 2 * np.pi
                t2 = e / MT_LENGTH * 2 * np.pi
                ts = np.linspace(t1, t2 if s <= e else t2 + 2 * np.pi, npts)
                fig.add_trace(go.Surface(
                    x=[r_inner_all * np.cos(ts), 1.05 * np.cos(ts)],
                    y=[r_inner_all * np.sin(ts), 1.05 * np.sin(ts)],
                    z=[[0] * npts, [0] * npts],
                    colorscale=[[0, color], [1, color]],
                    showscale=False,
                    opacity=0.25,
                    hoverinfo='skip',
                    showlegend=False
                ))
    else:
        # All genes (faint)
        for gname, s, e in genes:
            color = gene_colors[gname]
            npts = 25
            t1 = s / MT_LENGTH * 2 * np.pi
            t2 = e / MT_LENGTH * 2 * np.pi
            ts = np.linspace(t1, t2 if s <= e else t2 + 2 * np.pi, npts)
            fig.add_trace(go.Surface(
                x=[r_inner_all * np.cos(ts), 0.78 * np.cos(ts)],
                y=[r_inner_all * np.sin(ts), 0.78 * np.sin(ts)],
                z=[[0] * npts, [0] * npts],
                colorscale=[[0, color], [1, color]],
                showscale=False,
                opacity=0.10,
                hoverinfo='skip',
                showlegend=False
            ))

    # Gene labels
    for gname, s, e in genes:
        t_mid = ((s + e) / 2 / MT_LENGTH * 2 * np.pi) % (2 * np.pi)
        size = 15 if selected_gene == gname else 12
        fig.add_trace(go.Scatter3d(
            x=[0.88 * np.cos(t_mid)],
            y=[0.88 * np.sin(t_mid)],
            z=[0.25],
            mode='text',
            text=[gname],
            textfont=dict(size=size, color='black', family='Arial'),
            showlegend=False,
            hoverinfo='skip'
        ))

    # Final layout
    fig.update_layout(
        title=dict(
            text=f"3D mtDNA Heteroplasmy Viewer ({tsv_file})",
            x=0.5,
            y=0.95,
            font=dict(size=24, family='Arial Black')
        ),
        scene=dict(
            xaxis=dict(visible=False, range=[-1.5, 1.5]),
            yaxis=dict(visible=False, range=[-1.5, 1.5]),
            zaxis=dict(
                title='Samples with variant (>1%)',
                range=[-0.5, NUM_SAMPLES * 1.1],
                showgrid=True,
                gridcolor='lightgray',
                showline=False,
                backgroundcolor='white',
                title_font=dict(size=16)
            ),
            camera=dict(eye=dict(x=2.8, y=2.5, z=1.2)),
            aspectmode='cube'
        ),
        width=1200,
        height=850,
        showlegend=False,
        paper_bgcolor='white',
        plot_bgcolor='white',
        margin=dict(l=40, r=40, t=80, b=60, pad=0)
    )

    return fig

# Interactive Dash app
app = Dash(__name__, external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css'])
opts = [{'label': 'All genes', 'value': ''}] + [
    {'label': g[0], 'value': g[0]} for g in genes
]

app.layout = html.Div([
    html.Div([
        html.H1("3D mtDNA Heteroplasmy Viewer", 
                style={'textAlign': 'center', 'margin': '20px 0 10px 0', 
                       'color': '#2c3e50', 'fontWeight': 'bold'}),
        html.Div(f"File: {tsv_file} | Heteroplasmy Ratio (>1%): {total_sites} | Samples: {NUM_SAMPLES}",
                style={'textAlign': 'center', 'fontSize': '18px', 'marginBottom': '30px', 
                       'color': '#7f8c8d'}),
        dcc.Dropdown(
            id='gene_select',
            options=opts,
            value='',
            style={'width': '50%', 'margin': '0 auto 40px 0', 'display': 'block'}
        ),
        dcc.Graph(id='main_plot', 
                 style={'width': '100%', 'height': '850px', 
                       'margin': '0 auto', 'display': 'block'})
    ], style={
        'maxWidth': '1400px', 
        'margin': '0 auto', 
        'padding': '30px 20px',
        'backgroundColor': 'white',
        'boxShadow': '0 4px 20px rgba(0,0,0,0.1)',
        'borderRadius': '12px'
    })
], style={'minHeight': '100vh', 'background': 'linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%)'})

@app.callback(Output('main_plot', 'figure'), Input('gene_select', 'value'))
def update(gene):
    return create_figure(gene)

if __name__ == "__main__":
    print(f"Server running at http://127.0.0.1:8050/")
    print(f"Loaded: {tsv_file} | {total_sites} variant sites | {NUM_SAMPLES} samples")
    print("Hover markers for details. Select gene to highlight region.")
    app.run(debug=True, host='127.0.0.1', port=8050, dev_tools_hot_reload=False)

