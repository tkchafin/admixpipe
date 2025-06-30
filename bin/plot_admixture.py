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


def make_plot(df, output_html, template_file=None):
    import plotly.express as px
    from pathlib import Path

    # Melt to long format
    df_long = df.melt(
        id_vars=["Individual", "Population"],
        var_name="Cluster",
        value_name="Proportion",
    )

    # Sort ancestry segments by proportion within each individual
    df_long = df_long.sort_values(
        by=["Individual", "Proportion"], ascending=[True, False]
    )

    # Maintain input order of individuals
    df_long["Individual"] = pd.Categorical(
        df_long["Individual"], categories=df["Individual"], ordered=True
    )

    # Create BrBg color scale
    brbg_scale = px.colors.diverging.Spectral
    num_clusters = df.shape[1] - 2  # subtract Individual and Population
    color_seq = px.colors.sample_colorscale(
        brbg_scale, [i / max(1, num_clusters - 1) for i in range(num_clusters)]
    )

    # Plot
    fig = px.bar(
        df_long,
        x="Individual",
        y="Proportion",
        color="Cluster",
        color_discrete_sequence=color_seq,
        hover_data=["Individual", "Population", "Cluster", "Proportion"],
    )

    # Remove spacing and borders between bars
    fig.update_traces(marker_line_width=0)

    # Population tick labels
    pop_counts = df["Population"].value_counts(sort=False)
    pop_positions = pop_counts.cumsum() - pop_counts / 2

    fig.update_layout(
        barmode="stack",
        xaxis=dict(
            tickmode="array",
            tickvals=pop_positions.values,
            ticktext=pop_positions.index,
            showgrid=False,
            title="Population",
        ),
        yaxis=dict(
            title="Ancestry Proportion",
            range=[0, 1],
            showgrid=False,
        ),
        margin=dict(t=60, b=100),
        title=dict(text="ADMIXTURE Ancestry Barplot", x=0.5),
        legend_title="Cluster",
        template="simple_white",
    )

    # Save with optional template header
    html_body = fig.to_html(full_html=False, include_plotlyjs="cdn")
    if template_file:
        with open(template_file) as t:
            header = t.read()
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

    args = parser.parse_args()
    df = load_data(args.clumpp, args.inds, args.pops)
    make_plot(df, args.out, args.template)
