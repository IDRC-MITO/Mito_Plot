import pandas as pd
import numpy as np
import plotly.graph_objects as go
import sys

def clean_and_convert(value):
    if isinstance(value, str):
        return float(value.strip('%')) / 100
    return float(value) / 100 if value <= 100 else float(value) / 10000

def plot_heteroplasmy(tsv_file, output_html="Mitochondrial_Heteroplasmy.html"):
    # Load data
    data = pd.read_csv(tsv_file, sep='\t')
    MT_LENGTH = 16569

    fig = go.Figure()

    samples = data.columns[4:]
    for sample in samples:
        positions = data['Position']
        heteroplasmies = data[sample].apply(clean_and_convert)
        thetas = positions / MT_LENGTH * 360

        hover_texts = [
            f"Sample: {sample}<br>Pos: {pos}<br>Heteroplasmy: {het*100:.2f}%"
            for pos, het in zip(positions, heteroplasmies)
        ]

        fig.add_trace(go.Scatterpolar(
            r=heteroplasmies,
            theta=thetas,
            mode='markers',
            marker=dict(size=6),
            name=sample,
            text=hover_texts,
            hoverinfo='text'
        ))

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

    for gene, start, end in genes:
        start_angle = start / MT_LENGTH * 360
        end_angle = end / MT_LENGTH * 360
        
        thetas = np.linspace(start_angle, end_angle, 100)
        rs = [1.05] * len(thetas)
        
        fig.add_trace(go.Scatterpolar(
            r=rs,
            theta=thetas,
            mode='lines',
            line=dict(width=8, color='lightgrey'),
            name=gene,
            hoverinfo='text',
            text=[f"{gene}: {start}-{end}"] * len(thetas),
            showlegend=False
        ))
        
        mid_angle = (start_angle + end_angle) / 2
        fig.add_trace(go.Scatterpolar(
            r=[1.10],
            theta=[mid_angle],
            mode="text",
            text=[gene],
            textposition="middle center",
            showlegend=False
        ))

    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 1.1],
                tickvals=[0, 0.25, 0.5, 0.75, 1],
                ticktext=['0%', '25%', '50%', '75%', '100%']
            ),
            angularaxis=dict(
                direction="clockwise",
                rotation=90,
                showticklabels=False,
                ticks='',
                showline=False,
                showgrid=True
            ),
        ),
        showlegend=False,
        title="Mitochondrial Heteroplasmy and Gene Locations (Interactive)"
    )

    fig.write_html(output_html)
    print(f"Interactive plot saved as '{output_html}'.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run_step2.py <heteroplasmy_tsv_file>")
        sys.exit(1)

    tsv_file = sys.argv[1]
    plot_heteroplasmy(tsv_file)
