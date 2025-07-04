#!/usr/bin/env python3
import pandas as pd
import plotly.express as px
import plotly.io as pio
import argparse
from pathlib import Path

def load_data(qmat_file, ind_file, pop_file):
    # Parse Q matrix (after ":")
    q_raw = pd.read_csv(qmat_file, sep=":", header=None)
    q_df = q_raw[1].str.strip().str.split(expand=True).astype(float)

    # Load sample IDs and populations
    individuals = pd.read_csv(ind_file, header=None)[0]
    populations = pd.read_csv(pop_file, header=None)[0]

    # Sanity check
    if not (len(individuals) == len(populations) == len(q_df)):
        raise ValueError("Mismatch in number of rows across input files.")

    q_df.columns = [f"Cluster {i+1}" for i in range(q_df.shape[1])]
    q_df["Individual"] = individuals
    q_df["Population"] = populations

    return q_df


def make_plot(df, output_html, template_file=None, palette="Spectral", sort_pop=False, sort_ind=False):
    cluster_cols = [c for c in df.columns if c.startswith("Cluster")]

    # Determine population order
    if sort_pop:
        pop_means = df.groupby('Population')[cluster_cols].mean()
        pop_dom = pop_means.idxmax(axis=1)
        sorted_pops = sorted(pop_dom.index, key=lambda p: int(pop_dom[p].split()[1]))
    else:
        sorted_pops = df['Population'].unique().tolist()

    # Determine individual order
    if sort_ind:
        pop_dom = df.groupby('Population')[cluster_cols].mean().idxmax(axis=1)
        individual_order = []
        for pop in sorted_pops:
            dom = pop_dom[pop]
            sub = df[df['Population'] == pop].copy()
            sub = sub.sort_values(by=dom, ascending=False)
            individual_order.extend(sub['Individual'].tolist())
    else:
        seen = set()
        individual_order = []
        for ind in df['Individual']:
            if ind not in seen:
                seen.add(ind)
                individual_order.append(ind)

    # Melt to long format
    df_long = df.melt(
        id_vars=["Individual", "Population"],
        var_name="Cluster", value_name="Proportion"
    )
    df_long["Individual"] = pd.Categorical(
        df_long["Individual"], categories=individual_order, ordered=True
    )

    # Population tick positions
    pop_counts = df_long.drop_duplicates(subset=["Individual"]) \
                        .groupby('Population').size().reindex(sorted_pops)
    pop_positions = pop_counts.cumsum() - pop_counts / 2

    # Build color sequence from any palette module
    num_clusters = len(cluster_cols)
    # Qualitative: list of discrete colors
    if hasattr(px.colors.qualitative, palette):
        color_seq = getattr(px.colors.qualitative, palette)
    # Sequential: continuous scale sampled
    elif hasattr(px.colors.sequential, palette):
        scale = getattr(px.colors.sequential, palette)
        color_seq = px.colors.sample_colorscale(
            scale, [i / max(1, num_clusters - 1) for i in range(num_clusters)]
        )
    # Diverging: continuous scale sampled
    elif hasattr(px.colors.diverging, palette):
        scale = getattr(px.colors.diverging, palette)
        color_seq = px.colors.sample_colorscale(
            scale, [i / max(1, num_clusters - 1) for i in range(num_clusters)]
        )
    else:
        raise ValueError(
            f"Palette '{palette}' not found in plotly.colors modules (qualitative, sequential, diverging)."
        )

    # Plot
    fig = px.bar(
        df_long,
        x="Individual", y="Proportion", color="Cluster",
        category_orders={"Individual": individual_order},
        color_discrete_sequence=color_seq,
        hover_data=["Individual", "Population", "Cluster", "Proportion"]
    )
    fig.update_traces(marker_line_width=0)

    fig.update_layout(
        barmode="stack",
        xaxis=dict(
            tickmode="array",
            tickvals=pop_positions.values,
            ticktext=pop_positions.index,
            showgrid=False,
            title="Population"
        ),
        yaxis=dict(
            title="Ancestry Proportion",
            range=[0, 1],
            showgrid=False
        ),
        margin=dict(t=60, b=100),
        title=dict(text="ADMIXTURE Ancestry Barplot", x=0.5),
        legend_title="Cluster",
        template="simple_white",
    )

    # Save
    html_body = fig.to_html(full_html=False, include_plotlyjs="cdn")
    if template_file:
        header = Path(template_file).read_text()
        Path(output_html).write_text(header + html_body)
    else:
        Path(output_html).write_text(html_body)

    print(f"✅ Plot saved to: {output_html}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate an interactive ADMIXTURE barplot with Plotly."
    )
    parser.add_argument("--clumpp", required=True, help="Q matrix file (colon-separated ancestry values)")
    parser.add_argument("--inds", required=True, help="File with sample IDs (one per line)")
    parser.add_argument("--pops", required=True, help="File with population IDs (one per line)")
    parser.add_argument("--out", required=True, help="Output HTML file path")
    parser.add_argument("--template", help="Optional HTML template to prepend")
    parser.add_argument(
        "--palette",
        default="Spectral",
        help="Palette name from plotly.colors (qualitative, sequential or diverging)"
    )
    parser.add_argument(
        "--sort_pop", action="store_true",
        help="Sort populations by dominant cluster (default: False)"
    )
    parser.add_argument(
        "--sort_ind", action="store_true",
        help="Sort individuals within populations by dominant cluster ancestry (default: False)"
    )

    args = parser.parse_args()
    df = load_data(args.clumpp, args.inds, args.pops)
    make_plot(
        df, args.out, args.template,
        palette=args.palette,
        sort_pop=args.sort_pop,
        sort_ind=args.sort_ind
    )
