#!/usr/bin/env python3
"""
Interactive PCA scatterplot (PC1/PC2/PC3) with optional missing‑data opacity
---------------------------------------------------------------------------
* Individuals colored by Population
* Dropdown menu to choose PC axes (PC1 vs PC2, PC1 vs PC3, PC2 vs PC3)
* Button to toggle marker opacity inversely scaled by missing‑data proportion
  (low‑missing = fully opaque; high‑missing ≈ 0.25 opacity)
* Hover shows Sample, Population, PC1‑3 values, Missing proportion
* Axis lines at 0, axis labels as PCx
* Legend entries for “Low missing” and “High missing” appear only when
  “Show missing data” is clicked
Usage:
    plot_pca.py \
        --input   <pca_table.tsv> \
        --header  <multiqc_pca.html> \
        --output  snpio_pca_mqc.html
"""
import pandas as pd
import plotly.graph_objects as go
import argparse
import re

def load_pca(path):
    df = pd.read_csv(path, sep=r'\s+', engine='python', comment='#')
    df.columns = [c.strip() for c in df.columns]
    rename_map = {}
    for col in df.columns:
        if re.match(r'Axis\s*1', col, re.I):      rename_map[col] = 'PC1'
        elif re.match(r'Axis\s*2', col, re.I):    rename_map[col] = 'PC2'
        elif re.match(r'Axis\s*3', col, re.I):    rename_map[col] = 'PC3'
        elif re.match(r'Sample', col, re.I):      rename_map[col] = 'Sample'
        elif re.search(r'Population', col, re.I): rename_map[col] = 'Population'
        elif re.search(r'Missing', col, re.I):    rename_map[col] = 'Missing'
        elif re.match(r'Size', col, re.I):        rename_map[col] = 'Size'
    return df.rename(columns=rename_map)

def load_header(path):
    with open(path) as f:
        return f.read().strip()

def main(args):
    df = load_pca(args.input)
    pops = df['Population'].unique().tolist()
    n_traces = len(pops)

    # helper: list of values per population
    def grouped(col):
        return [df[df['Population']==pop][col].tolist() for pop in pops]

    # compute min/max missingness
    all_missing = df['Missing']
    min_miss, max_miss = all_missing.min(), all_missing.max()

    # scale_op: low missing → 1.0; high missing → 0.25
    def scale_op(m):
        if max_miss == min_miss:
            return 1.0
        return 1.0 - ((m - min_miss) / (max_miss - min_miss)) * 0.75

    # per‐population opacity arrays
    missing_opacities = grouped('Missing')
    missing_opacities = [
        [scale_op(m) for m in grp] for grp in missing_opacities
    ]

    # prepare PC axes data
    axes_data = {
        ('PC1','PC2'): (grouped('PC1'), grouped('PC2')),
        ('PC1','PC3'): (grouped('PC1'), grouped('PC3')),
        ('PC2','PC3'): (grouped('PC2'), grouped('PC3')),
    }

    # start figure, add one trace per population
    fig = go.Figure()
    for xi, yi, pop in zip(*axes_data[('PC1','PC2')], pops):
        grp = df[df['Population']==pop]
        cd = list(zip(
            grp['Sample'], grp['Population'],
            grp['PC1'], grp['PC2'], grp['PC3'],
            grp['Missing']
        ))
        fig.add_trace(go.Scatter(
            x=xi, y=yi,
            mode='markers',
            name=pop,
            marker=dict(size=8, opacity=1.0),
            customdata=cd,
            hovertemplate=(
                "Sample: %{customdata[0]}<br>"
                "Population: %{customdata[1]}<br>"
                "PC1: %{customdata[2]:.4f}<br>"
                "PC2: %{customdata[3]:.4f}<br>"
                "PC3: %{customdata[4]:.4f}<br>"
                "Missing: %{customdata[5]:.4f}<extra></extra>"
            )
        ))

    # dummy traces for missingness legend (hidden by default)
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(color='black', size=8, opacity=scale_op(min_miss)),
        name=f'Low missing ({min_miss:.2f})',
        visible=False
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(color='black', size=8, opacity=scale_op(max_miss)),
        name=f'High missing ({max_miss:.2f})',
        visible=False
    ))

    # dropdown for axis selection
    axis_buttons = []
    for (xa, ya), (xs, ys) in axes_data.items():
        axis_buttons.append(dict(
            label=f'{xa} vs {ya}',
            method='update',
            args=[
                {
                    'x': xs + [None, None],
                    'y': ys + [None, None],
                },
                {
                    'xaxis': {'title': xa},
                    'yaxis': {'title': ya},
                }
            ]
        ))

    # buttons for opacity toggle
    opacity_buttons = [
        dict(
            label='Default opacity',
            method='update',
            args=[
                {
                    'marker.opacity': [1.0]*n_traces + [scale_op(min_miss), scale_op(max_miss)],
                    'visible': [True]*n_traces + [False, False],
                },
                {}
            ]
        ),
        dict(
            label='Show missing data',
            method='update',
            args=[
                {
                    'marker.opacity': missing_opacities + [scale_op(min_miss), scale_op(max_miss)],
                    'visible': [True]*n_traces + [True, True],
                },
                {}
            ]
        ),
    ]

    # layout with axis lines at zero
    fig.update_layout(
        updatemenus=[
            dict(type='dropdown',
                 buttons=axis_buttons,
                 direction='down',
                 x=0, y=1.2,
                 xanchor='left', yanchor='top'),
            dict(type='buttons',
                 buttons=opacity_buttons,
                 direction='right',
                 x=0.5, y=1.2,
                 xanchor='left', yanchor='top'),
        ],
        autosize=True,
        margin=dict(l=80, r=20, t=140, b=80),
        paper_bgcolor='white',
        plot_bgcolor='white',
        legend_title_text='Population / Missingness',
        xaxis=dict(
            title='PC1',
            showline=True, linecolor='black', linewidth=1,
            zeroline=True, zerolinecolor='black', zerolinewidth=1,
            mirror=True
        ),
        yaxis=dict(
            title='PC2',
            showline=True, linecolor='black', linewidth=1,
            zeroline=True, zerolinecolor='black', zerolinewidth=1,
            mirror=True
        ),
        shapes=[
            # vertical line at PCx = 0
            dict(type='line',
                 xref='x', yref='paper',
                 x0=0, x1=0,
                 y0=0, y1=1,
                 line=dict(color='black', width=1)),
            # horizontal line at PCy = 0
            dict(type='line',
                 xref='paper', yref='y',
                 x0=0, x1=1,
                 y0=0, y1=0,
                 line=dict(color='black', width=1)),
        ],
        height=700, width=700
    )

    # write header + plot
    header_html = load_header(args.header)
    full_html = fig.to_html(full_html=True, include_plotlyjs='cdn')
    with open(args.output, 'w') as f:
        f.write(header_html + "\n" + full_html)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--input',  required=True, help="PCA table (TSV or whitespace-delimited)")
    p.add_argument('--header', required=True, help="HTML header/template file")
    p.add_argument('--output', required=True, help="Output HTML file")
    args = p.parse_args()
    main(args)
