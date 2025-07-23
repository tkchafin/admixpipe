#!/usr/bin/env python3
"""
Spatial map of ADMIXTURE average ancestry proportions for multiple K
-------------------------------------------------------------------
* Folium pie‑charts at sampling sites, sized by sample count
* One layer per K, selectable via LayerControl (radio‐style)
* Colors aligned across K
* Sample‐size legend
* Optional GeoJSON / shapefile overlays (with per‑layer style + z‑order)

Usage:
    plot_admixture_spatial_multiK.py \
      --indir       <clumpp_output_dir> \
      --inds        <individuals file> \
      --pops        <populations file> \
      --site_coords <ID,Lat,Lon TSV> \
      --out         multiK_map.html \
      [--min_k N] [--max_k M] \
      [--palette Spectral] [--basemap Esri.NatGeoWorldMap] \
      [--geo_data_json overlays.json] [--template header.html]
"""
import argparse, json, re, io, base64
from pathlib import Path

import pandas as pd
import geopandas as gpd
import folium
from branca.element import Element
import matplotlib.pyplot as plt
import plotly.express as px
import numpy as np

def parse_clumpp(dir_path):
    results = {}
    for f in Path(dir_path).rglob("ClumppIndFile.output.*"):
        m = re.search(r"\.output\.(\d+)$", f.name)
        if m:
            results[int(m.group(1))] = f
    return dict(sorted(results.items()))

def load_data(qmat_file, ind_file, pop_file):
    q = pd.read_csv(qmat_file, sep=":", header=None)
    df = q[1].str.strip().str.split(expand=True).astype(float)
    inds = pd.read_csv(ind_file, header=None, names=["Individual"])["Individual"]
    pops = pd.read_csv(pop_file, header=None, names=["Population"])["Population"]
    if not (len(df)==len(inds)==len(pops)):
        raise ValueError("Length mismatch among Q, inds, pops")
    df.columns = [f"Cluster {i+1}" for i in range(df.shape[1])]
    df["Individual"] = inds.values
    df["Population"] = pops.values
    return df

def load_coords(path):
    return pd.read_csv(path, sep="\t", header=None,
                       names=["ID","Latitude","Longitude"])

def aggregate_sites(df, coords, cluster_cols):
    by_ind = df.merge(coords, left_on="Individual", right_on="ID", how="left")
    by_pop = df.merge(coords, left_on="Population", right_on="ID", how="left")
    best   = by_ind if by_ind["Latitude"].notna().sum()>=by_pop["Latitude"].notna().sum() else by_pop
    best   = best.dropna(subset=["Latitude","Longitude"])
    best[cluster_cols] = best[cluster_cols].apply(pd.to_numeric, errors="coerce")
    grp    = best.groupby(["ID","Latitude","Longitude"])
    means  = grp[cluster_cols].mean().reset_index()
    counts = grp.size().rename("count").reset_index()
    return means.merge(counts, on=["ID","Latitude","Longitude"])

def rgb_to_hex(col):
    m = re.match(r"rgb\((\d+),\s*(\d+),\s*(\d+)\)", col)
    return f"#{int(m[1]):02x}{int(m[2]):02x}{int(m[3]):02x}" if m else col

def palette(name, n):
    if hasattr(px.colors.qualitative, name):
        seq = getattr(px.colors.qualitative, name)
    elif hasattr(px.colors.sequential, name):
        seq = px.colors.sample_colorscale(getattr(px.colors.sequential,name),
                                          [i/(n-1) for i in range(n)])
    elif hasattr(px.colors.diverging, name):
        seq = px.colors.sample_colorscale(getattr(px.colors.diverging,name),
                                          [i/(n-1) for i in range(n)])
    else:
        raise ValueError(f"palette '{name}' not found")
    if len(seq)<n:
        seq = (seq*((n//len(seq))+1))[:n]
    return [rgb_to_hex(c) for c in seq]

def pie_icon(vals, cols, px):
    fig,ax = plt.subplots(figsize=(1,1), dpi=px)
    ax.pie(vals, startangle=90, colors=cols)
    ax.axis("equal")
    buf = io.BytesIO()
    plt.savefig(buf, format="png", transparent=True,
                bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()

def load_overlays(path):
    if not path:
        return []
    obj = json.loads(Path(path).read_text())
    return [obj] if isinstance(obj,dict) else obj

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--indir",       required=True)
    p.add_argument("--inds",        required=True)
    p.add_argument("--pops",        required=True)
    p.add_argument("--site_coords", required=True)
    p.add_argument("--out",         required=True)
    p.add_argument("--min_k",       type=int)
    p.add_argument("--max_k",       type=int)
    p.add_argument("--palette",     default="Spectral")
    p.add_argument("--basemap",     default="Esri.NatGeoWorldMap")
    p.add_argument("--geo_data_json")
    p.add_argument("--template")
    args = p.parse_args()

    files = parse_clumpp(args.indir)
    Ks    = [k for k in files if
             (args.min_k is None or k>=args.min_k) and
             (args.max_k is None or k<=args.max_k)]
    if not Ks:
        raise ValueError("No K values after filtering")

    # load inds/pops
    inds = pd.read_csv(args.inds, header=None, names=["Individual"])["Individual"].tolist()
    pops = pd.read_csv(args.pops, header=None, names=["Population"])["Population"].tolist()

    coords = load_coords(args.site_coords)
    maxclus= max(Ks)
    cols_all = palette(args.palette, maxclus)

    # compute global sample‐size legend
    # use first K to get counts
    df0    = load_data(str(files[Ks[0]]), args.inds, args.pops)
    sites0 = aggregate_sites(df0, coords, [c for c in df0 if c.startswith("Cluster")])
    min_n, max_n = int(sites0["count"].min()), int(sites0["count"].max())
    # legend HTML
    size_html = ""
    for n,lbl in [(min_n,f"{min_n} sample{'s' if min_n>1 else ''}"),
                  (max_n,f"{max_n} samples")]:
        d = int( (np.sqrt(n/max_n)*(80-20)) + 20 )//2
        svg = d+6
        size_html += (
            f'<div style="display:flex;align-items:center;margin:2px;">'
            f'<svg width="{svg}" height="{svg}" style="margin-right:6px;">'
            f'<circle cx="{svg//2}" cy="{svg//2}" r="{d//2}" fill="#999"/></svg>{lbl}'
            '</div>'
        )
    legend = (
        '<div style="position:fixed;bottom:50px;left:10px;'
        'background:white;padding:8px;border:1px solid #ccc;'
        'font-size:12px;z-index:1000;">'
        '<b>Pie size</b><br/>' + size_html +
        '</div>'
    )

    # init map
    m = folium.Map(
        location=[coords.Latitude.mean(), coords.Longitude.mean()],
        zoom_start=6, tiles=args.basemap,
        control_scale=True, zoom_control=False,
        width="100%", height="100%"
    )
    # add legend
    m.get_root().html.add_child(Element(legend))

    # overlays
    for layer in load_overlays(args.geo_data_json):
        style = layer.get("style",{})
        gdf   = gpd.read_file(layer["path"]).to_crs(epsg=4326)
        folium.GeoJson(
            gdf.__geo_interface__,
            name=Path(layer["path"]).stem,
            style_function=lambda f,s=style: s
        ).add_to(m)

    # one radio‐style layer per K
    for idx,K in enumerate(sorted(Ks)):
        dfk = load_data(str(files[K]), args.inds, args.pops)
        cols= [c for c in dfk if c.startswith("Cluster")]
        sites = aggregate_sites(dfk, coords, cols)
        fg    = folium.FeatureGroup(name=f"K = {K}", show=(idx==0))
        mn    = sites["count"].max()
        for _,r in sites.iterrows():
            scale = np.sqrt(r["count"]/mn)
            px    = int(20 + (80-20)*scale)
            tip   = f"{r.ID} (N={int(r['count'])}): " + ", ".join(f"{c}={r[c]:.2f}" for c in cols)
            fg.add_child(folium.Marker(
                [r.Latitude, r.Longitude],
                icon=folium.CustomIcon(
                    pie_icon([r[c] for c in cols], cols_all[:len(cols)], px),
                    icon_size=(px//2,px//2)
                ),
                tooltip=tip
            ))
        m.add_child(fg)

    control = folium.LayerControl(collapsed=False)
    m.add_child(control)

    # convert checkboxes → radio
    m.get_root().html.add_child(Element("""
<script>
document.querySelectorAll('.leaflet-control-layers-overlays input').forEach((i)=>
  { i.type='radio'; i.name='radios'; }
);
</script>
"""))

    html = m._repr_html_()
    out  = (Path(args.template).read_text() if args.template else "") + html
    Path(args.out).write_text(out)
    print("✅ Multi‑K spatial map saved →", args.out)

if __name__=="__main__":
    main()
