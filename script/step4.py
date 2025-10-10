import pandas as pd
import numpy as np
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output
import matplotlib.cm as cm
import matplotlib.colors as mcolors

def clean_and_convert(value):
    if isinstance(value, str):
        return float(value.strip('%')) / 100
    return float(value) / 100 if value <= 100 else float(value) / 10000

def darken_color(hex_color, factor=0.6):
    rgb = np.array(mcolors.to_rgb(hex_color))
    dark_rgb = np.clip(rgb * factor, 0, 1)
    return mcolors.to_hex(dark_rgb)

def position_to_gene(pos, genes, MT_LENGTH):
    for gene, start, end in genes:
        if start <= end:
            if start <= pos <= end:
                return gene
        else:
            if pos >= start or pos <= end:
                return gene
    return None

def gene_center_camera_pos_and_center(gene, genes, MT_LENGTH, distance=0.7, height=0.5):
    if gene is None:
        eye = dict(x=2.2, y=2.2, z=1.5)
        center = dict(x=0, y=0, z=0)
        up = dict(x=0, y=0, z=1)
        return eye, center, up
    for g, start, end in genes:
        if g == gene:
            if start <= end:
                mid = (start + end) / 2
            else:
                mid = ((start + MT_LENGTH + end) / 2) % MT_LENGTH
            theta = mid / MT_LENGTH * 2 * np.pi
            center = dict(x=np.cos(theta), y=np.sin(theta), z=0)
            dir_vec = np.array([center['x'], center['y']])
            norm = np.linalg.norm(dir_vec)
            if norm > 0:
                eye_vec = -dir_vec / norm * distance
            else:
                eye_vec = np.array([0, 0])
            eye = dict(x=eye_vec[0], y=eye_vec[1], z=height)
            up = dict(x=0, y=0, z=1)
            return eye, center, up
    return dict(x=2.2, y=2.2, z=1.5), dict(x=0, y=0, z=0), dict(x=0, y=0, z=1)

def create_figure(selected_gene=None):
    tsv_file = "All_20250117.tsv"  #change Path
    data = pd.read_csv(tsv_file, sep='\t')
    MT_LENGTH = 16569
    samples = data.columns[4:]
    positions = data['Position'].values
    het_matrix = data[samples].applymap(clean_and_convert).fillna(0).values
    thetas = positions / MT_LENGTH * 2 * np.pi

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

    gene_names = [g[0] for g in genes]
    cmap = cm.get_cmap('tab20', len(gene_names))
    gene_colors = {g: mcolors.to_hex(cmap(i)) for i, g in enumerate(gene_names)}

    threshold = 0.01
    count_above_th = (het_matrix > threshold).sum(axis=1)
    max_count = max(count_above_th) if len(count_above_th) > 0 else 1
    desired_max_height = 1.0
    z_curve = (count_above_th / max_count * desired_max_height) if max_count > 0 else count_above_th

    fig = go.Figure()
    min_inner_radius = 0.15

    max_line_width = 6.0
    min_line_width = 0.5

    # 曲線のトレースを作成しつつ、全体のインデックス管理
    scatter_indices = []
    for i, pos in enumerate(positions):
        theta = thetas[i]
        z_height = z_curve[i]
        het_values = np.sort(het_matrix[i, :])
        het_values = het_values[het_values > 0]
        if len(het_values) == 0:
            scatter_indices.append(None)
            continue

        region_gene = position_to_gene(pos, genes, MT_LENGTH)
        base_color = gene_colors.get(region_gene, "#000000") if region_gene else "#000000"

        # 選択された遺伝子かどうか判定
        visible_bool = True if (selected_gene is None or selected_gene == "" or region_gene == selected_gene) else False
        factor = 0.6 if visible_bool else 0.3
        color = darken_color(base_color, factor=factor)

        r_values = het_values * (1 - min_inner_radius) + min_inner_radius
        x_values = r_values * np.cos(theta)
        y_values = r_values * np.sin(theta)
        z_values = np.linspace(0, z_height, len(r_values))

        normalized_z = z_height / desired_max_height if desired_max_height > 0 else 0
        line_width = normalized_z * (max_line_width - min_line_width) + min_line_width

        trace = go.Scatter3d(
            x=x_values, y=y_values, z=z_values,
            mode="lines",
            line=dict(color=color, width=line_width),
            hovertemplate=(
                f"Position: {pos}<br>"
                f"Gene: {region_gene}<br>"
                f"Count: {count_above_th[i]}<br>"
                f"Het (min–max): {r_values.min():.3f}–{r_values.max():.3f}<br>"
                f"Height: {z_height:.2f}<extra></extra>"
            ),
            visible=visible_bool,
            showlegend=False
        )
        fig.add_trace(trace)
        scatter_indices.append(trace)

    # 遺伝子領域のMesh3dトレースを追加
    for gene, start, end in genes:
        start_theta = start / MT_LENGTH * 2 * np.pi
        end_theta = end / MT_LENGTH * 2 * np.pi
        thetas_band = np.linspace(start_theta, end_theta, 100)
        inner_r = min_inner_radius
        outer_r = 1.0
        x_inner = inner_r * np.cos(thetas_band)
        y_inner = inner_r * np.sin(thetas_band)
        x_outer = outer_r * np.cos(thetas_band)
        y_outer = outer_r * np.sin(thetas_band)
        z_zero = np.zeros_like(thetas_band)

        x_band = np.concatenate([x_inner, x_outer])
        y_band = np.concatenate([y_inner, y_outer])
        z_band = np.concatenate([z_zero, z_zero])

        n = len(thetas_band)
        i_idx, j_idx, k_idx = [], [], []
        for idx in range(n - 1):
            i_idx.append(idx)
            j_idx.append(n + idx)
            k_idx.append(n + idx + 1)
            i_idx.append(idx)
            j_idx.append(n + idx + 1)
            k_idx.append(idx + 1)

        visible_bool = True if (selected_gene is None or selected_gene == "" or gene == selected_gene) else False
        opacity = 0.5 if visible_bool else 0.15
        color = gene_colors.get(gene, "lightgrey")

        fig.add_trace(go.Mesh3d(
            x=x_band, y=y_band, z=z_band,
            i=i_idx, j=j_idx, k=k_idx,
            opacity=opacity,
            color=color,
            name=gene,
            hovertext=f"{gene}: {start}-{end}",
            hoverinfo="text",
            showscale=False,
            visible=visible_bool
        ))

        mid_theta = (start_theta + end_theta) / 2
        label_r = outer_r + 0.13
        label_x = label_r * np.cos(mid_theta)
        label_y = label_r * np.sin(mid_theta)
        label_z = 0.09 * desired_max_height

        fig.add_trace(go.Scatter3d(
            x=[label_x], y=[label_y], z=[label_z],
            mode="text",
            text=[gene],
            visible=visible_bool,
            textfont=dict(size=13, color="black", family="Arial"),
            showlegend=False,
            hoverinfo="none"
        ))

    eye, center, up = gene_center_camera_pos_and_center(selected_gene if selected_gene else None, genes, MT_LENGTH)

    fig.update_layout(
        legend=dict(
            title="Genes",
            itemsizing='constant',
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="black",
            borderwidth=1,
        ),
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(title='Sample Count (scaled)', range=[-0.1, desired_max_height + 0.1]),
            aspectmode='data',
            camera=dict(eye=eye, center=center, up=up),
        ),
        title="Radial Heteroplasmy (Internal Curves for Multiple Samples)",
        margin=dict(l=0, r=0, b=0, t=50),
        paper_bgcolor="white",
    )

    return fig


app = Dash(__name__)

gene_options = [{"label": "ALL", "value": ""}]
genes_list = [
    "D-loop", "MT-RNR1", "MT-RNR2", "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2",
    "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5",
    "MT-ND6", "MT-CYB"
]
for g in genes_list:
    gene_options.append({"label": g, "value": g})

app.layout = html.Div([
    html.H2("Mitochondrial Heteroplasmy Radial 3D Plot"),
    dcc.Dropdown(
        id="gene-selector",
        options=gene_options,
        value="",
        clearable=False,
        style={"width": "300px"}
    ),
    dcc.Graph(id="het3d-graph", style={"height": "900px"})
])

@app.callback(
    Output("het3d-graph", "figure"),
    Input("gene-selector", "value")
)
def update_figure(selected_gene):
    gene = selected_gene if selected_gene != "" else None
    return create_figure(gene)

if __name__ == "__main__":
    app.run(debug=True)

