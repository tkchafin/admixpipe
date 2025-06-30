#!/usr/bin/env python3
import re
import argparse
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio

def parse_template(template_file):
    """Extract metadata from the MultiQC-style HTML comment block."""
    with open(template_file, "r") as f:
        content = f.read()
    match = re.search(r"<!--(.*?)-->", content, re.DOTALL)
    if not match:
        raise ValueError("No metadata block found in template.")
    raw_block = match.group(1)
    meta = {}
    for line in raw_block.strip().splitlines():
        if ":" in line:
            key, value = line.strip().split(":", 1)
            meta[key.strip()] = value.strip().strip("\"'")
    return meta

def build_comment(meta_dict):
    lines = ["<!--"]
    for key, value in meta_dict.items():
        lines.append(f'{key}: "{value}"')
    lines.append("-->")
    return "\n".join(lines)

def generate_plot(input_file, output_file, header_comment):
    # 1) read & sort
    df = pd.read_csv(input_file, delim_whitespace=True)
    df['K'] = pd.to_numeric(df['K'], errors='raise')
    df = df.sort_values('K')

    # 2) build the figure
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=df["K"],
            y=df["Mean"],
            error_y=dict(type="data", array=df["StDev"], visible=True),
            mode="lines+markers",
            name="Mean CV Error",
            marker=dict(size=8),
            line=dict(width=2),
        )
    )

    fig.update_layout(
        title="Cross-validation Error by K",
        xaxis_title="K",
        yaxis_title="Mean CV Error",
        template="plotly_white",
    )

    # 3) write out with your MultiQC comment header
    html = pio.to_html(fig, full_html=False, include_plotlyjs="cdn")
    with open(output_file, "w") as f:
        f.write(header_comment + "\n" + html)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate CV error plot as HTML with embedded MultiQC metadata."
    )
    parser.add_argument(
        "input_file", help="Input summary file with K, Mean, StDev columns"
    )
    parser.add_argument(
        "--output", default="cvplot_mqc.html", help="Output HTML filename"
    )
    parser.add_argument(
        "--template",
        required=True,
        help="Path to HTML file containing MultiQC metadata comment block",
    )
    parser.add_argument("--id", help="Override for the 'id' field in the metadata")
    parser.add_argument(
        "--title", help="Override for the 'title' field in the metadata"
    )
    parser.add_argument("--section_name", help="Override for 'section_name'")
    parser.add_argument("--description", help="Override for 'description'")

    args = parser.parse_args()

    # Load and override template
    metadata = parse_template(args.template)
    for fld in ("id", "title", "section_name", "description"):
        val = getattr(args, fld)
        if val:
            metadata[fld] = val

    header_comment = build_comment(metadata)
    generate_plot(args.input_file, args.output, header_comment)
