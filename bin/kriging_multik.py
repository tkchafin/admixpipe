#!/usr/bin/env python3
"""
Spatial kriging of ADMIXTURE average ancestry proportions → GeoTIFFs + MultiQC-ready HTML (Multi-K)
---------------------------------------------------------------------------------------------------
* EXACT SAME parsing as your pie-chart script (Q/inds/pops/coords + overlays)
* Reads multiple K solutions from a CLUMPP directory (ClumppIndFile.output.K)
* Kriging per K (PyKrige exponential, theta≈fields::Krig range; falls back to IDW), strictly serial
* Folium HTML maps on Leaflet basemaps (no axes), with:
    - auto-fit to raster extent
    - overlays from --geo_data_json (styled; optional z_order respected)
    - north arrow (↑ N) and styled scale bar
    - drop-down selector for K that toggles radio-style layers via LayerControl
* Outputs (SAME FILENAMES, now true-color RGB to match Plotly):
  - GeoTIFFs: one set per K into --geotiff_dir, each prefixed "K#_"
      - K#_<prefix>_kriging_discrete.tif        (RGB UInt8; colors baked in)
      - K#_<prefix>_kriging_simpson.tif         (RGB UInt8; viridis baked in)
  - HTML (Multi-K in one map each; K chosen via dropdown):
      - <prefix>_kriging_discrete.html
      - <prefix>_kriging_simpson.html
  - <prefix>.tsv                                (site summary; from first K loaded)
"""
import argparse, json, re, io, base64, sys, os, warnings
from pathlib import Path

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.ops import unary_union

# Folium (HTML maps)
import folium
from branca.element import Element

# Colormaps for PNG/arrays
import matplotlib.cm as mplcm
import matplotlib.colors as mplcolors

# GeoTIFF writing
import rasterio
from rasterio.transform import from_bounds

# ───────────────────────── optional deps / shims ─────────────────────────
def try_pykrige():
    try:
        from pykrige.ok import OrdinaryKriging  # noqa
        return True
    except Exception:
        return False

def encode_png(arr_uint8):
    """
    Encode a (H,W,3) RGB or (H,W,4) RGBA uint8 array to base64 PNG.
    Tries Pillow; falls back to matplotlib.
    """
    try:
        from PIL import Image  # type: ignore
        mode = "RGBA" if arr_uint8.shape[-1] == 4 else "RGB"
        im = Image.fromarray(arr_uint8, mode=mode)
        buf = io.BytesIO()
        im.save(buf, format="PNG")
        return base64.b64encode(buf.getvalue()).decode()
    except Exception:
        try:
            import matplotlib.pyplot as plt  # type: ignore
            buf = io.BytesIO()
            plt.imsave(buf, arr_uint8, format="png")
            return base64.b64encode(buf.getvalue()).decode()
        except Exception as e:
            raise RuntimeError("Unable to encode PNG: need Pillow or matplotlib") from e

# ───────────────────────── helpers (parsing retained) ─────────────────────────
def parse_clumpp(dir_path: str, kmin=None, kmax=None):
    """Return {K:int -> Path} for files matching ClumppIndFile.output.K."""
    results = {}
    for f in Path(dir_path).rglob("ClumppIndFile.output.*"):
        m = re.search(r"\.output\.(\d+)$", f.name)
        if m:
            K = int(m.group(1))
            if (kmin is not None and K < kmin) or (kmax is not None and K > kmax):
                continue
            results[K] = f
    return dict(sorted(results.items()))

def load_q_ind_pop(qmat_file, ind_file, pop_file):
    q = pd.read_csv(qmat_file, sep=":", header=None, dtype=str, engine="python")
    df = q[1].astype(str).str.strip().str.split(expand=True).astype(float)
    inds = pd.read_csv(ind_file, header=None, names=["Individual"])["Individual"]
    pops = pd.read_csv(pop_file, header=None, names=["Population"])["Population"]
    if not (len(df) == len(inds) == len(pops)):
        raise ValueError("Length mismatch among Q, inds, pops")
    df.columns = [f"Cluster {i+1}" for i in range(df.shape[1])]
    df["Individual"] = inds.values
    df["Population"] = pops.values
    return df

def load_coords(path):
    return pd.read_csv(path, sep="\t", header=None, names=["ID","Latitude","Longitude"])

def merge_coords(df, coords):
    by_ind = df.merge(coords, left_on="Individual", right_on="ID", how="left")
    by_pop = df.merge(coords, left_on="Population", right_on="ID", how="left")
    best   = by_ind if by_ind["Latitude"].notna().sum() >= by_pop["Latitude"].notna().sum() else by_pop
    return best.dropna(subset=["Latitude","Longitude"])

def aggregate_sites(df, coords, cluster_cols):
    df_coords = merge_coords(df, coords)
    df_coords[cluster_cols] = df_coords[cluster_cols].apply(pd.to_numeric, errors="coerce")
    means = df_coords.groupby(["ID","Latitude","Longitude"], as_index=False)[cluster_cols].mean()
    counts = df_coords.groupby(["ID","Latitude","Longitude"], as_index=False).agg(count=("Individual","size"))
    df_sites = means.merge(counts, on=["ID","Latitude","Longitude"])
    return df_sites

def rgb_to_hex(col):
    m = re.match(r"rgb\((\d+),\s*(\d+),\s*(\d+)\)", col)
    return f"#{int(m[1]):02x}{int(m[2]):02x}{int(m[3]):02x}" if m else col

def palette(name, n):
    if name.lower() in ("spectral", "viridis", "plasma", "inferno", "magma", "cividis"):
        cmap = mplcm.get_cmap(name)
        return [mplcolors.to_hex(cmap(i/(max(n-1,1)))) for i in range(n)]
    cmap = mplcm.get_cmap("tab20")
    return [mplcolors.to_hex(cmap(i % 20)) for i in range(n)]

def load_overlays(path):
    if not path:
        return []
    data = json.loads(Path(path).read_text())
    if isinstance(data, dict):
        data = [data]
    return sorted(data, key=lambda d: int(d.get("z_order", 0)))

# ───────────────────────── kriging & raster helpers ─────────────────────────
def compute_bounds(df_sites, overlays):
    geoms = []
    for layer in overlays:
        gdf = gpd.read_file(layer["path"]).to_crs(epsg=4326)
        if not gdf.empty:
            geoms.append(unary_union(gdf.geometry))
    if geoms:
        xmin, ymin, xmax, ymax = unary_union(geoms).bounds
    else:
        xmin, xmax = float(df_sites["Longitude"].min()), float(df_sites["Longitude"].max())
        ymin, ymax = float(df_sites["Latitude"].min()),  float(df_sites["Latitude"].max())
    dx, dy = (xmax - xmin), (ymax - ymin)
    pad_x, pad_y = (dx * 0.02 if dx > 0 else 0.1), (dy * 0.02 if dy > 0 else 0.1)
    return (xmin - pad_x, ymin - pad_y, xmax + pad_x, ymax + pad_y)

def build_grid(bounds, nx, ny):
    xmin, ymin, xmax, ymax = bounds
    xs = np.linspace(xmin, xmax, nx)
    ys = np.linspace(ymin, ymax, ny)
    xx, yy = np.meshgrid(xs, ys)
    return xs, ys, xx, yy

def krige_one_band(args):
    (name, x, y, z, xs, ys, method) = args
    ny, nx = len(ys), len(xs)
    grid = np.full((ny, nx), np.nan, dtype=float)

    if np.isnan(z).all() or np.nanstd(z) < 1e-9 or len(z) < 3:
        grid[:] = np.nanmean(z)
        return name, grid

    if method == "pykrige":
        try:
            from pykrige.ok import OrdinaryKriging
            # Emulate fields::Krig with exponential covariance; theta ≈ range parameter
            OK = OrdinaryKriging(
                x, y, z,
                variogram_model="exponential",
                variogram_parameters=[1.0, float(args_global.theta), float(args_global.nugget)],  # [sill, range, nugget]
                verbose=False
            )
            zhat, _ = OK.execute("grid", xs, ys)
            grid = np.asarray(zhat, dtype=float)
            np.clip(grid, 0.0, 1.0, out=grid)
            return name, grid
        except Exception:
            pass

    # IDW fallback
    xp = x.reshape(-1, 1, 1)
    yp = y.reshape(-1, 1, 1)
    wp = z.reshape(-1, 1, 1)
    XX = xs.reshape(1, 1, -1)
    YY = ys.reshape(1, -1, 1)
    dx = XX - xp
    dy = YY - yp
    dist2 = dx*dx + dy*dy
    with np.errstate(divide="ignore"):
        w = 1.0 / dist2
    exact = dist2 == 0.0
    if exact.any():
        w[exact] = 1e12
    wsum = np.sum(w, axis=0)
    zsum = np.sum(w * wp, axis=0)
    grid = zsum / wsum
    grid[~np.isfinite(grid)] = np.nanmean(z)
    np.clip(grid, 0.0, 1.0, out=grid)
    return name, grid

def simpson_diversity(arr_stack):
    with np.errstate(invalid="ignore"):
        s2 = np.nansum(np.square(arr_stack), axis=0)
        D = 1.0 - s2
        D[~np.isfinite(D)] = np.nan
    return D

def grid_edges(xs, ys):
    if len(xs) < 2 or len(ys) < 2:
        raise ValueError("grid must be at least 2x2 to compute edges")
    x_edges = np.concatenate([[xs[0] - (xs[1]-xs[0])/2], (xs[:-1]+xs[1:])/2, [xs[-1] + (xs[-1]-xs[-2])/2]])
    y_edges = np.concatenate([[ys[0] - (ys[1]-ys[0])/2], (ys[:-1]+ys[1:])/2, [ys[-1] + (ys[-1]-ys[-2])/2]])
    return x_edges, y_edges

# ─────────────── GeoTIFF writers ───────────────
def write_rgb_geotiff(xs, ys, rgb_uint8, path, mask_alpha=None):
    # rgb_uint8: (ny, nx, 3) uint8 in array orientation; will flip for raster origin
    ny, nx, _ = rgb_uint8.shape
    x_edges, y_edges = grid_edges(xs, ys)
    west, south, east, north = x_edges[0], y_edges[0], x_edges[-1], y_edges[-1]
    transform = from_bounds(west, south, east, north, nx, ny)
    arr = rgb_uint8[::-1, :, :]  # GDAL top-left origin
    r, g, b = arr[..., 0], arr[..., 1], arr[..., 2]
    with rasterio.open(
        path, "w",
        driver="GTiff",
        height=ny, width=nx, count=3,
        dtype="uint8", crs="EPSG:4326", transform=transform,
        tiled=True, compress="deflate", photometric="RGB"
    ) as dst:
        dst.write(r, 1); dst.write(g, 2); dst.write(b, 3)
        if mask_alpha is not None:
            m = (mask_alpha[::-1, :] > 0).astype("uint8") * 255
            dst.write_mask(m)

def write_geotiff_single(xs, ys, band, path, dtype, nodata, colormap=None, descr=None):
    # kept for completeness; not used now that we write RGB to match Plotly exactly
    ny, nx = band.shape
    x_edges, y_edges = grid_edges(xs, ys)
    west, south, east, north = x_edges[0], y_edges[0], x_edges[-1], y_edges[-1]
    transform = from_bounds(west, south, east, north, nx, ny)
    arr = np.where(np.isfinite(band), band, nodata).astype(dtype)
    arr = arr[::-1, :]
    _dtype = np.dtype(dtype).name.lower()
    predictor = 3 if _dtype in ("float32", "float64") else 2
    with rasterio.open(
        path, "w",
        driver="GTiff",
        height=ny, width=nx, count=1,
        dtype=arr.dtype, crs="EPSG:4326", transform=transform,
        nodata=nodata, tiled=True, compress="deflate", predictor=predictor
    ) as dst:
        dst.write(arr, 1)
        if descr:
            dst.set_band_description(1, descr)
        if (colormap is not None) and (_dtype == "uint8"):
            dst.write_colormap(1, colormap)

# ───────────────────────── image encoding + array builders ─────────────────────────
def scalar_to_rgba_array(z, vmin=0.0, vmax=1.0, cmap_name="viridis", alpha=220):
    cmap = mplcm.get_cmap(cmap_name)
    norm = mplcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    z_mask = np.isfinite(z)
    rgba = np.zeros((z.shape[0], z.shape[1], 4), dtype=np.uint8)
    if z_mask.any():
        mapped = (cmap(norm(np.where(z_mask, z, vmin))) * 255).astype(np.uint8)
        rgba[z_mask] = mapped[z_mask]
        rgba[z_mask, 3] = alpha
        rgba[~z_mask, 3] = 0
    return rgba  # not flipped here

def scalar_to_rgba_png(z, vmin=0.0, vmax=1.0, cmap_name="viridis", alpha=220):
    rgba = scalar_to_rgba_array(z, vmin=vmin, vmax=vmax, cmap_name=cmap_name, alpha=alpha)
    return "data:image/png;base64," + encode_png(rgba[::-1, :, :])

def classes_to_rgba_array(classes, class_hex, nodata_val=0, alpha=230):
    ny, nx = classes.shape
    rgba = np.zeros((ny, nx, 4), dtype=np.uint8)
    for idx, hx in enumerate(class_hex, start=1):
        r = int(hx[1:3], 16); g = int(hx[3:5], 16); b = int(hx[5:7], 16)
        mask = classes == idx
        rgba[mask, 0] = r; rgba[mask, 1] = g; rgba[mask, 2] = b; rgba[mask, 3] = alpha
    rgba[classes == nodata_val, 3] = 0
    return rgba  # not flipped here

def classes_to_rgba_png(classes, class_hex, nodata_val=0, alpha=230):
    rgba = classes_to_rgba_array(classes, class_hex, nodata_val, alpha)
    return "data:image/png;base64," + encode_png(rgba[::-1, :, :])

# ───────────────────────── Folium map builders (HTML only) ─────────────────────────
def _styled_map(basemap_name, center_lat, center_lon):
    try:
        m = folium.Map(location=[center_lat, center_lon], zoom_start=6, tiles=basemap_name,
                       control_scale=True, zoom_control=False, width="100%", height="100%")
    except Exception:
        m = folium.Map(location=[center_lat, center_lon], zoom_start=6, tiles="CartoDB Positron",
                       control_scale=True, zoom_control=False, width="100%", height="100%")
    m.get_root().header.add_child(Element("""
      <style>
        html,body,#map{height:100%!important;margin:0;}
        .leaflet-control-zoomslider,.leaflet-control-zoom{display:none!important;}
        .leaflet-control-scale{background:rgba(255,255,255,0.9);font-size:10px;padding:2px 6px;}
        #map{border:2px solid #000;}
        .k-select-wrap{position:absolute;top:8px;left:10px;z-index:1000;background:rgba(255,255,255,0.95);padding:6px;border:1px solid #ccc;border-radius:6px;font-size:12px;}
        .k-select{margin-left:6px;}
      </style>
    """))
    m.get_root().html.add_child(Element('<div style="position:absolute;top:8px;right:10px;z-index:999;font-size:1.2em;">↑ N</div>'))
    return m

def _add_overlays(m, overlays):
    for layer in overlays:
        style = layer.get("style", {})
        gdf   = gpd.read_file(layer["path"]).to_crs(epsg=4326)
        folium.GeoJson(gdf.__geo_interface__, name=Path(layer["path"]).stem,
                       style_function=lambda f, s=style: s).add_to(m)

def _add_samples_black_dots(m, samples_df):
    if samples_df is None or samples_df.empty:
        return
    for _, r in samples_df.iterrows():
        folium.CircleMarker(location=[float(r["Latitude"]), float(r["Longitude"])],
                            radius=3, color="#000000", weight=0.5, fill=True, fill_opacity=1.0).add_to(m)

def _fit_to_bounds(m, bounds):
    west, south, east, north = bounds
    m.fit_bounds([[south, west], [north, east]])

def _base_map_from_bounds(bounds, basemap_name):
    west, south, east, north = bounds
    m = _styled_map(basemap_name, center_lat=(south+north)/2, center_lon=(west+east)/2)
    _fit_to_bounds(m, bounds)
    return m, [[south, west], [north, east]]

def save_folium_with_template(m, template_path, out_path):
    header = Path(template_path).read_text() if template_path else ""
    html = m._repr_html_()
    Path(out_path).write_text(header + html)
    print(f"✅ HTML written with MultiQC header → {out_path}")

def add_dropdown_for_k(m, k_values, label_text="Model (K)"):
    options = "".join([f'<option value="{k}">K = {k}</option>' for k in k_values])
    m.get_root().html.add_child(Element(f"""
<div class="k-select-wrap">
  <b>{label_text}</b>
  <select id="k-select" class="k-select">{options}</select>
</div>
<script>
(function() {{
  const select = document.getElementById('k-select');
  function chooseK(k) {{
    const labels = Array.from(document.querySelectorAll('.leaflet-control-layers-overlays label'));
    const target = labels.find(l => l.textContent.trim() === 'K = ' + k);
    if (!target) return;
    const radio = target.querySelector('input[type="radio"]');
    if (radio && !radio.checked) {{
      radio.click();
    }}
  }}
  chooseK(select.value);
  select.addEventListener('change', (e) => chooseK(e.target.value));
}})();
</script>
"""))

# ───────────────────────── main (serial, multi-K) ─────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir",       required=True)
    ap.add_argument("--min_k",       type=int)
    ap.add_argument("--max_k",       type=int)
    ap.add_argument("--inds",        required=True)
    ap.add_argument("--pops",        required=True)
    ap.add_argument("--site_coords", required=True)
    ap.add_argument("--geo_data_json")
    ap.add_argument("--template_disc")
    ap.add_argument("--template_div")
    ap.add_argument("--out_prefix",  required=True)
    ap.add_argument("--out_dir")
    ap.add_argument("--geotiff_dir")
    ap.add_argument("--table_out")
    ap.add_argument("--palette",     default="Spectral")
    ap.add_argument("--basemap",     default="Esri.NatGeoWorldMap")
    ap.add_argument("--grid_nx",     type=int, default=400)
    ap.add_argument("--grid_ny",     type=int, default=400)
    ap.add_argument("--no_pykrige",  action="store_true")
    # TESS-like kriging controls (good defaults):
    ap.add_argument("--theta",  type=float, default=None,
                    help="Exponential variogram range (≈ fields::Krig theta). If unset, auto = 0.3 * max(map span).")
    ap.add_argument("--nugget", type=float, default=1e-10,
                    help="Variogram nugget (very small keeps surfaces smooth).")
    ap.add_argument("--majority_px", type=int, default=2,
                    help="Post-classification majority filter radius in pixels (0 disables).")
    args = ap.parse_args()

    # expose selected params to krige_one_band / postproc without changing call signatures
    class _ArgsGlobal: pass
    global args_global
    args_global = _ArgsGlobal()
    args_global.theta = args.theta   # finalized after bounds known
    args_global.nugget = args.nugget
    args_global.majority_px = args.majority_px

    prefix = Path(args.out_prefix)
    html_outdir  = Path(args.out_dir) if args.out_dir else prefix.parent
    html_outdir.mkdir(parents=True, exist_ok=True)
    stem = prefix.name

    geotiff_dir = Path(args.geotiff_dir) if args.geotiff_dir else (prefix.parent / "geotiff")
    geotiff_dir.mkdir(parents=True, exist_ok=True)

    def OP_html(suffix): return html_outdir / f"{stem}_{suffix}.html"
    tsv_path  = Path(args.table_out) if args.table_out else (html_outdir / f"{stem}.tsv")

    kfiles = parse_clumpp(args.indir, args.min_k, args.max_k)
    if not kfiles:
        raise ValueError("No ClumppIndFile.output.K found after filtering. Check --indir / K range.")
    Ks = sorted(kfiles.keys())
    Kmax = max(Ks)

    coords = load_coords(args.site_coords)
    overlays = load_overlays(args.geo_data_json) if args.geo_data_json else []
    colors_by_k = {k: palette(args.palette, Kmax)[:k] for k in Ks}

    firstK = Ks[0]
    df0 = load_q_ind_pop(kfiles[firstK], args.inds, args.pops)
    cluster_cols0 = [c for c in df0.columns if c.startswith("Cluster")]
    df_sites0 = aggregate_sites(df0, coords, cluster_cols0)
    if df_sites0.empty:
        raise ValueError("No valid site coordinates found")
    df_sites0[["ID","Latitude","Longitude","count"] + cluster_cols0].to_csv(tsv_path, sep="\t", index=False)
    print(f"✅ Site summary TSV saved to: {tsv_path}")

    # Grid + map bounds
    bounds = compute_bounds(df_sites0, overlays)
    xs, ys, xx, yy = build_grid(bounds, args.grid_nx, args.grid_ny)
    x_edges, y_edges = grid_edges(xs, ys)
    west, south, east, north = x_edges[0], y_edges[0], x_edges[-1], y_edges[-1]
    map_bounds = (west, south, east, north)

    # finalize theta if None: ~30% of max span for smoothness
    if args_global.theta is None:
        xmin, ymin, xmax, ymax = bounds
        span = max(xmax - xmin, ymax - ymin)
        args_global.theta = 0.3 * span

    method = "pykrige" if (not args.no_pykrige and try_pykrige()) else "idw"
    if method == "idw":
        print("ℹ️ PyKrige unavailable/disabled; using IDW fallback.", file=sys.stderr)

    disc_layers, simp_layers = [], []

    # small helper: optional majority filter (no ambiguous/zero labels)
    def apply_majority_filter(labels, K, r):
        if r <= 0:
            return labels
        try:
            from scipy.ndimage import generic_filter
            yy, xx = np.ogrid[-r:r+1, -r:r+1]
            footprint = (xx*xx + yy*yy) <= (r*r)

            def _mode_filter(window):
                vals = window.astype(np.uint16)
                m = np.bincount(vals, minlength=K+1)  # labels 0..K (but we avoid 0 elsewhere)
                if m[1:].max() >= m[0]:
                    return np.argmax(m[1:]) + 1
                return np.argmax(m[1:]) + 1
            return generic_filter(labels.astype(np.uint16), _mode_filter, footprint=footprint, mode="nearest")
        except Exception:
            pad = np.pad(labels, 1, mode="edge")
            out = labels.copy()
            for i in range(labels.shape[0]):
                for j in range(labels.shape[1]):
                    block = pad[i:i+3, j:j+3].ravel()
                    m = np.bincount(block, minlength=K+1)
                    out[i, j] = np.argmax(m[1:]) + 1
            return out

    for K in Ks:
        qpath = kfiles[K]
        dfK = load_q_ind_pop(qpath, args.inds, args.pops)
        cluster_cols = [c for c in dfK.columns if c.startswith("Cluster")]
        df_sites = aggregate_sites(dfK, coords, cluster_cols)

        x = df_sites["Longitude"].to_numpy(float)
        y = df_sites["Latitude"].to_numpy(float)

        tasks = [(name, x, y, df_sites[name].to_numpy(float), xs, ys, method) for name in cluster_cols]

        bands = {}
        for t in tasks:
            nm, grid = krige_one_band(t)
            bands[nm] = grid

        stack = np.stack([bands[n] for n in cluster_cols], axis=0)
        stack = np.clip(stack, 0.0, 1.0)
        sums = np.nansum(stack, axis=0, keepdims=True)
        with np.errstate(invalid="ignore", divide="ignore"):
            over = sums > 1.0001
            if np.any(over):
                stack[:, over[0]] = stack[:, over[0]] / sums[:, over[0]]

        # robust argmax (never produce ambiguous/zero)
        any_valid = np.any(np.isfinite(stack), axis=0)
        tmp = np.where(np.isfinite(stack), stack, -np.inf)
        argmax_idx = np.argmax(tmp, axis=0)          # 0..K-1
        disc = (argmax_idx + 1).astype(np.uint32)    # 1..K
        if not any_valid.all():
            disc[~any_valid] = disc[any_valid].min() if np.any(any_valid) else 1
        if args_global.majority_px > 0:
            disc = apply_majority_filter(disc, K, args_global.majority_px)

        simpson = simpson_diversity(stack)

        base = f"K{K}_{stem}"
        colors = colors_by_k[K]

        # Build the SAME RGBA arrays used for the Plotly/Folium overlays
        disc_rgba = classes_to_rgba_array(disc.astype(int), class_hex=colors, nodata_val=0, alpha=255)
        disc_rgb = disc_rgba[..., :3]
        finite = np.isfinite(simpson)
        vmin = float(np.nanmin(simpson[finite])) if finite.any() else 0.0
        vmax = float(np.nanmax(simpson[finite])) if finite.any() else 1.0
        if not np.isfinite(vmin) or not np.isfinite(vmax) or abs(vmax - vmin) < 1e-9:
            vmin, vmax = 0.0, 1.0
        sim_rgba = scalar_to_rgba_array(simpson, vmin=vmin, vmax=vmax, cmap_name="viridis", alpha=255)
        sim_rgb = sim_rgba[..., :3]

        # Write EXACT visual as RGB GeoTIFFs, keeping SAME FILENAMES
        write_rgb_geotiff(xs, ys, disc_rgb, str(geotiff_dir / f"{base}_kriging_discrete.tif"),
                          mask_alpha=disc_rgba[..., 3])
        print(f"✅ GeoTIFF (discrete argmax, RGB): {geotiff_dir / f'{base}_kriging_discrete.tif'}")

        write_rgb_geotiff(xs, ys, sim_rgb, str(geotiff_dir / f"{base}_kriging_simpson.tif"),
                          mask_alpha=sim_rgba[..., 3])
        print(f"✅ GeoTIFF (Simpson diversity, RGB): {geotiff_dir / f'{base}_kriging_simpson.tif'}")

        # HTML overlays (same visuals)
        url_disc = "data:image/png;base64," + encode_png(disc_rgba[::-1, :, :])
        url_div  = "data:image/png;base64," + encode_png(sim_rgba[::-1, :, :])

        disc_layers.append((K, url_disc))
        simp_layers.append((K, url_div))

    def build_multik_html(layers, layer_label, template_path, out_path):
        m, bounds = _base_map_from_bounds(map_bounds, args.basemap)
        for idx, (K, url) in enumerate(sorted(layers, key=lambda t: t[0])):
            fg = folium.FeatureGroup(name=f"K = {K}", show=(idx == 0))
            folium.raster_layers.ImageOverlay(image=url, bounds=bounds, opacity=1.0,
                                              interactive=False, cross_origin=False).add_to(fg)
            fg.add_to(m)
        _add_overlays(m, overlays)
        _add_samples_black_dots(m, df_sites0[["Longitude","Latitude","ID","count"] + cluster_cols0])
        folium.LayerControl(collapsed=False).add_to(m)
        add_dropdown_for_k(m, [K for (K, _) in sorted(layers, key=lambda t: t[0])], label_text=layer_label)
        save_folium_with_template(m, template_path, out_path)

    build_multik_html(disc_layers, "Model (K) — Discrete", args.template_disc, (html_outdir / f"{stem}_kriging_discrete.html"))
    build_multik_html(simp_layers, "Model (K) — Simpson",  args.template_div,  (html_outdir / f"{stem}_kriging_simpson.html"))

    print("🎉 All per-K rasters written to:", geotiff_dir)
    print("🎉 Multi-K HTML maps:")
    print("   •", html_outdir / f"{stem}_kriging_discrete.html")
    print("   •", html_outdir / f"{stem}_kriging_simpson.html")

if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=UserWarning)
    main()
