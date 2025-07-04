#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def parse_clumpp(dir_path):
    """
    Recursively find all ClumppIndFile.output.K files under dir_path,
    parse each line as "idx: q1 q2 ... qK", split off the proportions,
    and return a dict K->DataFrame of floats.
    """
    results = {}
    base = Path(dir_path)
    for qfile in base.rglob('ClumppIndFile.output.*'):
        m = re.search(r'ClumppIndFile\.output\.(\d+)$', qfile.name)
        if not m:
            continue
        K = int(m.group(1))

        # Read into two string columns: idx and proportions
        q_raw = pd.read_csv(
            qfile,
            sep=':',
            header=None,
            names=['idx', 'proportions'],
            dtype=str,
            comment='#',
            engine='python'
        )

        # Split the proportions field on whitespace, convert to float
        q_df = (
            q_raw['proportions']
            .str.strip()
            .str.split(expand=True)
            .astype(float)
        )

        if q_df.shape[1] != K:
            raise ValueError(
                f'For K={K}, expected {K} proportions but parsed {q_df.shape[1]}'
            )

        q_df.columns = [f'Cluster {i+1}' for i in range(K)]
        results[K] = q_df

    return dict(sorted(results.items()))


def make_multiK_plot(data_dict,
                     out_html,
                     template=None,
                     palette='Spectral',
                     sort_pop=False,
                     sort_ind=False,
                     bestk=None,
                     global_inds=None,
                     global_pops=None,
                     min_k=None,
                     max_k=None):
    """
    data_dict: mapping K->DataFrame of cluster proportions
    global_inds, global_pops: lists of length = number of samples
    """
    # 1) Filter K range
    Ks = sorted(
        k for k in data_dict
        if (min_k is None or k >= min_k)
        and (max_k is None or k <= max_k)
    )
    if not Ks:
        raise ValueError('No K values remain after applying --min-k/--max-k')

    # 2) Assign Individual/Population
    for K in Ks:
        df = data_dict[K]
        if len(global_inds) != len(df) or len(global_pops) != len(df):
            raise ValueError(f'Length of inds/pops does not match rows for K={K}')
        df['Individual'] = global_inds
        df['Population'] = global_pops

    # 3) Determine ordering based on reference K
    refK = bestk if (bestk in Ks) else max(Ks)
    ref_df = data_dict[refK]
    cluster_cols_ref = [c for c in ref_df.columns if c.startswith('Cluster')]

    # Population order
    if sort_pop:
        pop_means = ref_df.groupby('Population')[cluster_cols_ref].mean()
        pop_dom = pop_means.idxmax(axis=1)
        sorted_pops = sorted(
            pop_dom.index,
            key=lambda p: int(pop_dom[p].split()[1])
        )
    else:
        sorted_pops = list(ref_df['Population'].unique())

    # Individual order
    if sort_ind:
        pop_dom = ref_df.groupby('Population')[cluster_cols_ref].mean().idxmax(axis=1)
        individual_order = []
        for pop in sorted_pops:
            dom = pop_dom[pop]
            sub = ref_df[ref_df['Population'] == pop].copy()
            sub = sub.sort_values(by=dom, ascending=False)
            individual_order.extend(sub['Individual'].tolist())
    else:
        seen = set()
        individual_order = []
        for ind in ref_df['Individual']:
            if ind not in seen:
                seen.add(ind)
                individual_order.append(ind)

    # 4) Compute population tick positions
    pop_counts = pd.Series(
        {pop: (ref_df['Population'] == pop).sum() for pop in sorted_pops},
        index=sorted_pops
    )
    pop_positions = pop_counts.cumsum() - pop_counts / 2

    # 5) Build full color palette using max clusters
    max_clus = max(Ks)
    if hasattr(px.colors.qualitative, palette):
        base = getattr(px.colors.qualitative, palette)
        if len(base) < max_clus:
            raise ValueError(f"Qualitative palette '{palette}' has only {len(base)} colors but need {max_clus}.")
        colors_full = base[:max_clus]
    elif hasattr(px.colors.sequential, palette):
        scale = getattr(px.colors.sequential, palette)
        colors_full = px.colors.sample_colorscale(
            scale,
            [i / (max_clus - 1) for i in range(max_clus)]
        )
    elif hasattr(px.colors.diverging, palette):
        scale = getattr(px.colors.diverging, palette)
        colors_full = px.colors.sample_colorscale(
            scale,
            [i / (max_clus - 1) for i in range(max_clus)]
        )
    else:
        raise ValueError(f"Palette '{palette}' not found in plotly.colors modules.")

    # 6) Create subplots
    fig = make_subplots(
        rows=len(Ks), cols=1,
        shared_xaxes=True,
        vertical_spacing=0.05,   # more room so "K = x" labels don't overlap
        subplot_titles=[f'K = {K}' for K in Ks]
    )

    # 7) Add bars with a single, complete legend
    legend_shown = set()
    for i, K in enumerate(Ks, start=1):
        df_plot = data_dict[K].set_index('Individual').reindex(individual_order).reset_index()
        cluster_cols = [c for c in df_plot.columns if c.startswith('Cluster')]
        colors = colors_full[:len(cluster_cols)]
        for col, colr in zip(cluster_cols, colors):
            show = col not in legend_shown
            if show:
                legend_shown.add(col)
            fig.add_trace(
                go.Bar(
                    x=df_plot['Individual'],
                    y=df_plot[col],
                    name=col,
                    legendgroup=col,
                    marker_color=colr,
                    marker_line_width=0,    # remove segment borders
                    showlegend=show
                ),
                row=i, col=1
            )

    # 8) Styling & axes
    fig.update_layout(
        barmode='stack',
        height=200 * len(Ks),
        width=900,
        template='simple_white',
        margin=dict(t=100, b=100),
        legend=dict(
            title='Cluster',
            orientation='v',
            x=1.02,
            y=1,
            yanchor='top'
        )
    )

    # remove individual y-axis titles and set range
    fig.update_yaxes(range=[0, 1], title_text='')

    # hide x-ticklabels for all but last row
    for r in range(1, len(Ks)):
        fig.update_xaxes(showticklabels=False, row=r, col=1)

    # last row: population ticks & x-axis title
    fig.update_xaxes(
        tickmode='array',
        tickvals=pop_positions.values,
        ticktext=pop_positions.index,
        showgrid=False,
        title_text='Population',
        row=len(Ks), col=1
    )

    # single shared y-axis label, moved further left and larger
    fig.add_annotation(
        text='Ancestry Proportion',
        xref='paper', yref='paper',
        x=-0.1, y=0.5,                # shift left
        showarrow=False,
        textangle=-90,
        font=dict(size=14)            # larger font
    )

    # 9) Render
    html = fig.to_html(full_html=False, include_plotlyjs='cdn')
    if template:
        header = Path(template).read_text().rstrip('\n') + '\n'
        Path(out_html).write_text(header + html)
    else:
        Path(out_html).write_text(html)

    print(f"✅ Multi-K plot saved to {out_html}")


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Plot multi-K admixture barplots')
    p.add_argument('--indir',   required=True, help='Input directory (recursive)')
    p.add_argument('--out',     required=True, help='Output HTML file')
    p.add_argument('--template', help='HTML template to prepend')
    p.add_argument('--inds',    required=True, help='Global individuals file (one per line)')
    p.add_argument('--pops',    required=True, help='Global populations file (one per line)')
    p.add_argument('--palette', default='Spectral', help='Palette name')
    p.add_argument('--sort_pop', action='store_true', help='Sort populations by dominant cluster')
    p.add_argument('--sort_ind', action='store_true', help='Sort individuals by reference K ancestry')
    p.add_argument('--bestk',    type=int,          help='Reference K (default = max)')
    p.add_argument('--mink',    dest='mink', type=int, help='Minimum K to include')
    p.add_argument('--maxk',    dest='maxk', type=int, help='Maximum K to include')

    args = p.parse_args()
    data = parse_clumpp(args.indir)
    global_inds = pd.read_csv(args.inds, header=None)[0].tolist()
    global_pops = pd.read_csv(args.pops, header=None)[0].tolist()

    make_multiK_plot(
        data, args.out, args.template, args.palette,
        args.sort_pop, args.sort_ind, args.bestk,
        global_inds, global_pops,
        args.mink, args.maxk
    )
