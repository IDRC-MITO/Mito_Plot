import pandas as pd
import numpy as np
import plotly.graph_objects as go
import sys

def clean_and_convert(value):
   
    if isinstance(value, str):
        return float(value.strip('%')) / 100
    return float(value) / 100 if value <= 100 else float(value) / 10000

def plot_heteroplasmy_circular_3d_compact_z(tsv_file, output_html="Mitochondrial_Heteroplasmy_Circular3D_CompactZ.html"):
    data = pd.read_csv(tsv_file, sep='\t')
    MT_LENGTH = 16569

    samples = data.columns[4:]
    positions = data['Position']

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

    fig = go.Figure()
    N = len(samples)

    for i, sample in enumerate(samples):
        heteroplasmies = data[sample].apply(clean_and_convert).fillna(0)
        thetas = positions / MT_LENGTH * 2 * np.pi
        r = heteroplasmies
        x = r * np.cos(thetas)
        y = r * np.sin(thetas)
        # The Z-axis is evenly spaced between 0 and 1
        z = np.full_like(r, 0 if N == 1 else i / (N - 1))

        hover_texts = [
            f"Sample: {sample}<br>Position: {pos}<br>Heteroplasmy: {het*100:.2f}%"
            for pos, het in zip(positions, heteroplasmies)
        ]

        fig.add_trace(go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(size=4, opacity=0.7),
            text=hover_texts,
            hoverinfo='text',
            name=sample
        ))

    # Circular bands in the gene region (fixed at Z=0)
    for gene, start, end in genes:
        start_angle = start / MT_LENGTH * 2 * np.pi
        end_angle = end / MT_LENGTH * 2 * np.pi
        thetas = np.linspace(start_angle, end_angle, 100)

        inner_r = 1.02
        outer_r = 1.08

        x_inner = inner_r * np.cos(thetas)
        y_inner = inner_r * np.sin(thetas)
        z_inner = np.zeros_like(thetas)

        x_outer = outer_r * np.cos(thetas)
        y_outer = outer_r * np.sin(thetas)
        z_outer = np.zeros_like(thetas)

        x_band = np.concatenate([x_inner, x_outer[::-1]])
        y_band = np.concatenate([y_inner, y_outer[::-1]])
        z_band = np.concatenate([z_inner, z_outer[::-1]])

        fig.add_trace(go.Mesh3d(
            x=x_band,
            y=y_band,
            z=z_band,
            opacity=0.3,
            color='lightgrey',
            name=gene,
            hovertext=f"{gene}: {start}-{end}",
            hoverinfo="text",
            showscale=False
        ))

        mid_angle = (start_angle + end_angle) / 2
        label_x = (outer_r + 0.05) * np.cos(mid_angle)
        label_y = (outer_r + 0.05) * np.sin(mid_angle)
        label_z = -0.1  

        fig.add_trace(go.Scatter3d(
            x=[label_x],
            y=[label_y],
            z=[label_z],
            mode="text",
            text=[gene],
            textfont=dict(size=9, color="black"),
            showlegend=False,
            hoverinfo="none"
        ))

    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(title='Sample Index (compressed)', range=[-0.1, 1.1]),
            aspectmode='data',
            camera=dict(eye=dict(x=1.5, y=1.5, z=0.3))
        ),
        showlegend=False,
        title="3D Circular Mitochondrial Heteroplasmy Plot (Compressed Z-axis)",
        margin=dict(l=0, r=0, b=0, t=50),
        paper_bgcolor="white"
    )

    fig.write_html(output_html)
    print(f"✅ 3D circular plot with compressed Z-axis saved as '{output_html}'.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run_3D_circular_compactz.py <heteroplasmy_tsv_file>")
        sys.exit(1)

    tsv_file = sys.argv[1]
    plot_heteroplasmy_circular_3d_compact_z(tsv_file)

