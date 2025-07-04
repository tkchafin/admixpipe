#!/usr/bin/env python3
import re
import argparse
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objs as go
import plotly.io as pio

def parse_template(template_file):
    """Extract metadata from the MultiQC-style HTML comment block."""
    with open(template_file, "r") as f:
        content = f.read()
    match = re.search(r"<!--(.*?)-->", content, re.DOTALL)
    if not match:
        raise ValueError("No metadata block found in template.")
    raw = match.group(1)
    meta = {}
    for line in raw.strip().splitlines():
        if ":" in line:
            key, val = line.strip().split(":", 1)
            meta[key.strip()] = val.strip().strip('"\'')
    return meta

def build_comment(meta):
    lines = ["<!--"]
    for k, v in meta.items():
        lines.append(f'{k}: "{v}"')
    lines.append("-->")
    return "\n".join(lines)

def generate_plot(input_file, output_file, header_comment, bestk=None):
    # 1) load and sort
    df = pd.read_csv(input_file, sep="\t")
    df = df.sort_values("K")

    # 2) build 2×2 subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=("Mean LnP(K)", "L′(K)", "L″(K)", "ΔK"),
        vertical_spacing=0.12, horizontal_spacing=0.1
    )

    # Mean LnP(K)
    fig.add_trace(go.Scatter(
        x=df["K"], y=df["Mean"], mode="lines+markers",
        error_y=dict(type="data", array=df["StDev"], visible=True),
        name="Mean LnP",
        line=dict(color='blue'),
        marker=dict(color='blue')
    ), row=1, col=1)

    # First-order change
    fig.add_trace(go.Scatter(
        x=df["K"], y=df["Lprime"], mode="lines+markers",
        error_y=dict(type="data", array=df["sd_Lprime"], visible=True),
        name="L′(K)",
        line=dict(color='blue'),
        marker=dict(color='blue')
    ), row=1, col=2)

    # Second-order change
    fig.add_trace(go.Scatter(
        x=df["K"], y=df["Lpp"], mode="lines+markers",
        error_y=dict(type="data", array=df["sd_Lpp"], visible=True),
        name="L″(K)",
        line=dict(color='blue'),
        marker=dict(color='blue')
    ), row=2, col=1)

    # ΔK as a line plot
    fig.add_trace(go.Scatter(
        x=df["K"], y=df["DeltaK"], mode="lines+markers",
        name="ΔK",
        line=dict(color='blue'),
        marker=dict(color='blue')
    ), row=2, col=2)

    # Annotate best K if provided
    if bestk is not None:
        fig.add_vline(
            x=bestk,
            line=dict(color='red', dash='dash'),
            annotation_text=f"K={bestk}",
            annotation_position="top right"
        )

    fig.update_layout(
        title_text="Evanno Method Results",
        showlegend=False,
        template="plotly_white",
        height=800, width=800
    )

    # 3) render and write
    html_div = pio.to_html(fig, full_html=False, include_plotlyjs="cdn")
    with open(output_file, "w") as f:
        f.write(header_comment + "\n" + html_div)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate Evanno ΔK multi-panel plot as HTML with MultiQC metadata."
    )
    parser.add_argument(
        "input_file",
        help="Evanno table (K, Mean, StDev, Lprime, sd_Lprime, Lpp, sd_Lpp, DeltaK)"
    )
    parser.add_argument(
        "-o", "--output", default="evanno_mqc.html",
        help="Output HTML file"
    )
    parser.add_argument(
        "--template", required=True,
        help="Path to HTML template containing the metadata comment block"
    )
    parser.add_argument(
        "--bestk", type=int,
        help="Best K value to annotate with a vertical line"
    )
    parser.add_argument("--id", help="Override the 'id' field in the metadata")
    parser.add_argument("--title", help="Override the 'title' field in the metadata")
    parser.add_argument("--section_name", help="Override the 'section_name' field")
    parser.add_argument("--description", help="Override the 'description' field")

    args = parser.parse_args()

    meta = parse_template(args.template)
    for fld in ("id", "title", "section_name", "description"):
        val = getattr(args, fld)
        if val:
            meta[fld] = val

    header_comment = build_comment(meta)
    generate_plot(args.input_file, args.output, header_comment, args.bestk)
