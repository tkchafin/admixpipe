#!/usr/bin/env python3
import pandas as pd
import plotly.graph_objects as go
import argparse
import re

def load_matrix(path):
    # allow space- or tab-delimited matrices
    df = pd.read_csv(path, sep=r'\s+', index_col=0, engine='python')
    df.index = df.index.astype(str).str.strip()
    df.columns = [str(c).strip() for c in df.columns]
    return df.astype(float)

def load_header(path):
    with open(path) as f:
        return f.read().strip()

def significance_stars(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return ''

def main(args):
    fst_df  = load_matrix(args.fst)
    pval_df = load_matrix(args.pvals)
    pops    = list(fst_df.index)
    n       = len(pops)

    # build data arrays
    z    = [[None]*n for _ in range(n)]
    text = [['']*n   for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i > j:
                fst = fst_df.iat[i, j]
                p   = pval_df.iat[i, j]
                z[i][j]    = fst
                text[i][j] = f"{fst:.3f}{significance_stars(p)}"

    header_html = load_header(args.header)

    # by default show labels
    fig = go.Figure(go.Heatmap(
        z=z,
        x=pops,
        y=pops,
        text=text,
        texttemplate='%{text}',
        textfont=dict(color='white', size=7),
        hovertemplate=(
            "Pop1: %{y}<br>"
            "Pop2: %{x}<br>"
            "F<sub>ST</sub>: %{z:.3f}<br>"
            "p-value: %{customdata:.2e}<extra></extra>"
        ),
        customdata=pval_df.values.tolist(),
        showscale=True,
        colorscale='Viridis'
    ))

    # toggle button
    blank = [['']*n for _ in range(n)]
    fig.update_layout(
        updatemenus=[dict(
            type='buttons',
            direction='right',
            x=0, y=1.15, xanchor='left', yanchor='top',
            active=0,
            buttons=[
                dict(label='Show Labels',
                     method='restyle',
                     args=[{'text': [text]}]
                ),
                dict(label='Hide Labels',
                     method='restyle',
                     args=[{'text': [blank]}]
                )
            ]
        )],
        autosize=True,
        margin=dict(l=80, r=20, t=60, b=80),
        paper_bgcolor='white',
        plot_bgcolor='white',
        height=800,
        width=800
    )

    fig.update_xaxes(tickangle=45, showgrid=False, zeroline=False)
    fig.update_yaxes(autorange='reversed', showgrid=False, zeroline=False)

    full_html = fig.to_html(full_html=True, include_plotlyjs='cdn')
    with open(args.output, 'w') as f:
        f.write(header_html + "\n")
        f.write(full_html)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--fst',    required=True)
    p.add_argument('--pvals',  required=True)
    p.add_argument('--header', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args()
    main(args)
