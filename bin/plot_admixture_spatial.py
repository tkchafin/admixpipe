#!/usr/bin/env python3
"""
Spatial map of ADMIXTURE average ancestry proportions
-----------------------------------------------------
* Folium pie‑charts at sampling sites, sized by sample count
* Mean Q‑matrix proportions per site (population or coordinate)
* Optional GeoJSON / shapefile overlays (with per‑layer style + z‑order)
* Legend for cluster colours and pie‑size scale
* Hover tooltip shows N + mean proportions per cluster

New options
-----------
--basemap     Name of the Leaflet/xyzservices tile layer (default
              'Esri.NatGeoWorldMap'). Any string accepted by
              folium.Map(tiles=...) is valid.
--table_out   Write site summary TSV (default: <out>.tsv)
"""
import geopandas as gpd
import pandas as pd
import folium
from branca.element import Element
import matplotlib.pyplot as plt
import plotly.express as px
import argparse, json, io, base64, re
from pathlib import Path
import numpy as np

# ───────────────────────── helpers ─────────────────────────
def load_data(qmat_file, ind_file, pop_file):
    # Parse Q matrix (after “:”)
    q_raw = pd.read_csv(qmat_file, sep=":", header=None)
    q_df  = q_raw[1].str.strip().str.split(expand=True).astype(float)
    # Load sample IDs and populations
    inds = pd.read_csv(ind_file, header=None, names=["Individual"])
    pops = pd.read_csv(pop_file, header=None, names=["Population"])
    if not (len(inds) == len(pops) == len(q_df)):
        raise ValueError("Row count mismatch among Q matrix, inds, pops")
    q_df.columns = [f"Cluster {i+1}" for i in range(q_df.shape[1])]
    q_df["Individual"] = inds["Individual"]
    q_df["Population"] = pops["Population"]
    return q_df

def load_coords(path):
    return pd.read_csv(path, sep="\t", header=None, names=["ID","Latitude","Longitude"])

def merge_coords(df, coords):
    by_ind = df.merge(coords, left_on="Individual", right_on="ID", how="left")
    by_pop = df.merge(coords, left_on="Population", right_on="ID", how="left")
    best   = by_ind if by_ind["Latitude"].notna().sum() >= by_pop["Latitude"].notna().sum() else by_pop
    return best.dropna(subset=["Latitude","Longitude"])

def aggregate_sites(df, coords, cluster_cols):
    # merge coordinates
    df_coords = merge_coords(df, coords)
    # compute mean proportions per site
    means = df_coords.groupby(["ID","Latitude","Longitude"], as_index=False)[cluster_cols].mean()
    # compute sample count per site
    counts = df_coords.groupby(["ID","Latitude","Longitude"], as_index=False) \
                     .agg(count=("Individual","size"))
    # merge together
    df_sites = means.merge(counts, on=["ID","Latitude","Longitude"])
    return df_sites

def rgb_to_hex(col):
    m = re.match(r"rgb\((\d+),\s*(\d+),\s*(\d+)\)", col)
    return f"#{int(m[1]):02x}{int(m[2]):02x}{int(m[3]):02x}" if m else col

def palette(name, n):
    if hasattr(px.colors.qualitative, name):
        seq = getattr(px.colors.qualitative, name)
    elif hasattr(px.colors.sequential, name):
        seq = px.colors.sample_colorscale(getattr(px.colors.sequential, name),
                                          [i/(n-1) for i in range(n)])
    elif hasattr(px.colors.diverging, name):
        seq = px.colors.sample_colorscale(getattr(px.colors.diverging, name),
                                          [i/(n-1) for i in range(n)])
    else:
        raise ValueError(f"palette '{name}' not found")
    if len(seq) < n:
        seq = (seq * ((n // len(seq)) + 1))[:n]
    return [rgb_to_hex(c) for c in seq]

def pie_icon(vals, cols, size_px):
    fig, ax = plt.subplots(figsize=(1,1), dpi=size_px)
    ax.pie(vals, startangle=90, colors=cols)
    ax.axis("equal")
    buf = io.BytesIO()
    plt.savefig(buf, format="png", transparent=True, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()

def load_overlays(json_path):
    if not json_path:
        return []
    data = json.loads(Path(json_path).read_text())
    if isinstance(data, dict):
        data = [data]
    return sorted(data, key=lambda d: int(d.get("z_order", 0)))

# ───────────────────────── main ─────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qmat",        required=True, help="ADMIXTURE .Q file (colon-separated)")
    ap.add_argument("--inds",        required=True, help="File with individual IDs (one per line)")
    ap.add_argument("--pops",        required=True, help="File with population/site IDs (one per line)")
    ap.add_argument("--site_coords", required=True, help="TSV with ID, Latitude, Longitude (no header)")
    ap.add_argument("--geo_data_json", help="Optional GeoJSON overlay descriptor")
    ap.add_argument("--template",    help="HTML header template to prepend")
    ap.add_argument("--out",         required=True, help="Output HTML file path")
    ap.add_argument("--palette",     default="Spectral", help="Palette name from plotly.colors")
    ap.add_argument("--min_pie_px",  type=int, default=20)
    ap.add_argument("--max_pie_px",  type=int, default=80)
    ap.add_argument("--basemap",     default="Esri.NatGeoWorldMap", help="Leaflet/xyzservices tile layer")
    ap.add_argument("--table_out",   help="Write site summary TSV (default: <out>.tsv)")
    args = ap.parse_args()

    # load and aggregate
    df = load_data(args.qmat, args.inds, args.pops)
    cluster_cols = [c for c in df.columns if c.startswith("Cluster")]
    coords       = load_coords(args.site_coords)
    df_sites     = aggregate_sites(df, coords, cluster_cols)
    if df_sites.empty:
        raise ValueError("No valid site coordinates found")

    # optional TSV
    tsv_path = Path(args.table_out) if args.table_out else Path(args.out).with_suffix(".tsv")
    df_sites[["ID","Latitude","Longitude","count"] + cluster_cols] \
        .to_csv(tsv_path, sep="\t", index=False)
    print(f"✅ Site summary TSV saved to: {tsv_path}")

    # Folium map
    try:
        m = folium.Map(
            location=[df_sites.Latitude.mean(), df_sites.Longitude.mean()],
            zoom_start=6,
            tiles=args.basemap,
            control_scale=True, zoom_control=False,
            width="100%", height="100%"
        )
    except Exception as e:
        print(f"⚠️ Could not load basemap '{args.basemap}' ({e}); falling back.")
        m = folium.Map(
            location=[df_sites.Latitude.mean(), df_sites.Longitude.mean()],
            zoom_start=6, tiles="CartoDB Positron",
            control_scale=True, zoom_control=False,
            width="100%", height="100%"
        )

    # minimal styling tweaks
    m.get_root().header.add_child(Element("""
      <style>
        html,body,#map{height:100%!important;margin:0;}
        .leaflet-control-zoomslider,.leaflet-control-zoom{display:none!important;}
        .leaflet-control-scale{background:rgba(255,255,255,0.9);font-size:10px;}
        #map{border:2px solid #000;}
      </style>"""))
    m.get_root().html.add_child(Element(
        '<div style="position:absolute;top:8px;right:10px;z-index:999;font-size:1.2em;">↑ N</div>'
    ))

    # overlays in z‑order
    for layer in load_overlays(args.geo_data_json):
        style = layer.get("style", {})
        gdf   = gpd.read_file(layer["path"]).to_crs(epsg=4326)
        folium.GeoJson(
            gdf.__geo_interface__,
            name=Path(layer["path"]).stem,
            style_function=lambda f, s=style: s
        ).add_to(m)

    # auto‑zoom
    m.fit_bounds([
        [df_sites.Latitude.min(), df_sites.Longitude.min()],
        [df_sites.Latitude.max(), df_sites.Longitude.max()],
    ])

    # pie markers
    colors = palette(args.palette, len(cluster_cols))
    max_n, min_n = df_sites["count"].max(), df_sites["count"].min()
    for _, row in df_sites.iterrows():
        scale   = np.sqrt(row["count"] / max_n)
        px_size = int(args.min_pie_px + (args.max_pie_px - args.min_pie_px) * scale)
        tooltip = (f"{row.ID} (N={int(row['count'])}): " +
                   ", ".join(f"{c}={row[c]:.2f}" for c in cluster_cols))
        folium.Marker(
            [row.Latitude, row.Longitude],
            icon=folium.CustomIcon(
                pie_icon([row[c] for c in cluster_cols], colors, px_size),
                icon_size=(px_size//2, px_size//2)
            ),
            tooltip=tooltip
        ).add_to(m)

    # legend
    def marker_dim(n):
        return int((args.min_pie_px + (args.max_pie_px - args.min_pie_px) * np.sqrt(n/max_n)) // 2)

    size_html = ""
    for n, lbl in [(min_n, f"{min_n} sample{'s' if min_n>1 else ''}"),
                   (max_n, f"{max_n} samples")]:
        d   = marker_dim(n)
        svg = d + 6
        size_html += (
            f'<div style="display:flex;align-items:center;margin-bottom:2px;">'
            f'<svg width="{svg}" height="{svg}" style="margin-right:6px;">'
            f'<circle cx="{svg//2}" cy="{svg//2}" r="{d//2}" fill="#999"/></svg>{lbl}'
            '</div>'
        )

    cluster_html = "".join(
        f'<div style="margin-bottom:2px;">'
        f'<i style="background:{col};display:inline-block;width:12px;height:12px;'
        f'margin-right:6px;"></i>{c}</div>'
        for c, col in zip(cluster_cols, colors)
    )

    legend_html = (
        f'<div style="border:1px solid #ccc;padding:10px;'
        f'background:#f8f8f8;font-size:0.9em;"><b>Cluster</b>{cluster_html}'
        f'<hr style="margin:6px 0;"><b>Pie size</b>{size_html}</div>'
    )

    wrapped = (
        f'<div style="border:1px solid #ddd;border-radius:4px;'
        f'padding:10px;margin-bottom:1em;display:flex;flex-wrap:wrap;">'
        f'<div style="flex:1 1 500px;min-width:400px;">{m._repr_html_()}</div>'
        f'<div style="flex:0 0 250px;margin-left:20px;">{legend_html}</div>'
        '</div>'
    )

    header = Path(args.template).read_text() if args.template else ""
    Path(args.out).write_text(header + wrapped)
    print(f"✅ Spatial ADMIXTURE map saved → {args.out}")

if __name__ == "__main__":
    main()
