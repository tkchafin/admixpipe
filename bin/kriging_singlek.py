#!/usr/bin/env python3
"""
Spatial kriging of ADMIXTURE average ancestry proportions → GeoTIFFs + MultiQC-ready HTML
-----------------------------------------------------------------------------------------
* EXACT SAME parsing as your pie-chart script (Q/inds/pops/coords + overlays)
* Kriging per K (PyKrige spherical if available; falls back to IDW), optionally parallel
* Folium HTML maps on Leaflet basemaps (no axes), with:
    - auto-fit to raster extent
    - overlays from --geo_data_json (styled; optional z_order respected)
    - north arrow (↑ N) and styled scale bar
    - radio-style toggle for multi-K continuous layers
* Outputs:
  - <prefix>_kriging_continuous.tif           (Float32, K bands; nodata=-9999)
  - <prefix>_kriging_discrete.tif             (UInt8/UInt16, 1 band; 0=nodata, 1..K=argmax class)
  - <prefix>_kriging_simpson.tif              (Float32, 1 band; nodata=-9999)
  - <prefix>_kriging_continuous_merged.tif    (Uint8 RGB, 3 bands + mask)
  - <prefix>_kriging_continuous.html          (Folium basemap + one layer per K, radio toggle)
  - <prefix>_kriging_continuous_merged.html   (Folium basemap + dominant-color continuous, gray mixing)
  - <prefix>_kriging_discrete.html            (Folium basemap + single layer)
  - <prefix>_kriging_simpson.html             (Folium basemap + single layer)
  - <prefix>.tsv                              (site summary)

New/changed CLI options
-----------------------
--out_prefix       Prefix for ALL outputs (recommended), e.g. outputs/wtd_admix_kriging
--template_cont    HTML header template to prepend to continuous HTML
--template_merged  HTML header template to prepend to merged continuous HTML
--template_disc    HTML header template to prepend to discrete HTML
--template_div     HTML header template to prepend to Simpson HTML
--basemap          Leaflet/xyzservices tiles name for HTML maps (default: Esri.NatGeoWorldMap)
--grid_nx/ny       Grid resolution in x/y (default 400/400)
--jobs             Parallel workers across K (0=auto, 1=serial)
--no_pykrige       Force IDW instead of PyKrige
--merge_thresh     Min winning proportion to color in merged view (default 0.5)
--merge_margin     Min (winner - runner_up) to color in merged view (default 0.0)
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

# Colormaps for PNG generation
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
    df_coords = merge_coords(df, coords)
    means = df_coords.groupby(["ID","Latitude","Longitude"], as_index=False)[cluster_cols].mean()
    counts = df_coords.groupby(["ID","Latitude","Longitude"], as_index=False) \
                     .agg(count=("Individual","size"))
    df_sites = means.merge(counts, on=["ID","Latitude","Longitude"])
    return df_sites

def rgb_to_hex(col):
    m = re.match(r"rgb\((\d+),\s*(\d+),\s*(\d+)\)", col)
    return f"#{int(m[1]):02x}{int(m[2]):02x}{int(m[3]):02x}" if m else col

def palette(name, n):
    # Use matplotlib colormaps; hex strings
    if name.lower() in ("spectral", "viridis", "plasma", "inferno", "magma", "cividis"):
        cmap = mplcm.get_cmap(name)
        cols = [mplcolors.to_hex(cmap(i/(max(n-1,1)))) for i in range(n)]
        return cols
    # As a robust default, sample tab20
    cmap = mplcm.get_cmap("tab20")
    return [mplcolors.to_hex(cmap(i % 20)) for i in range(n)]

def load_overlays(path):
    if not path:
        return []
    data = json.loads(Path(path).read_text())
    if isinstance(data, dict):
        data = [data]
    # preserve file order; respect optional z_order
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
    xx, yy = np.meshgrid(xs, ys)  # yy rows (lat increasing downwards)
    return xs, ys, xx, yy

def krige_one_band(args):
    """
    (name, x, y, z, xs, ys, method) -> (name, grid[ny,nx])
    """
    (name, x, y, z, xs, ys, method) = args
    ny, nx = len(ys), len(xs)
    grid = np.full((ny, nx), np.nan, dtype=float)

    # degenerate band
    if np.isnan(z).all() or np.nanstd(z) < 1e-9 or len(z) < 3:
        grid[:] = np.nanmean(z)
        return name, grid

    if method == "pykrige":
        try:
            from pykrige.ok import OrdinaryKriging
            OK = OrdinaryKriging(
                x, y, z,
                variogram_model="spherical",
                verbose=False, enable_plotting=False,
                coordinates_type="euclidean"
            )
            zhat, _ = OK.execute("grid", xs, ys)
            grid = np.asarray(zhat, dtype=float)
            np.clip(grid, 0.0, 1.0, out=grid)
            return name, grid
        except Exception:
            pass  # fall back

    # IDW fallback (power=2)
    xp = x.reshape(-1, 1, 1)
    yp = y.reshape(-1, 1, 1)
    wp = z.reshape(-1, 1, 1)
    XX = xs.reshape(1, 1, -1)
    YY = ys.reshape(1, -1, 1)
    dx = XX - xp
    dy = YY - yp
    dist2 = dx*dx + dy*dy
    with np.errstate(divide="ignore"):
        w = 1.0 / dist2  # 1/r^2
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
    # edges from centers
    if len(xs) < 2 or len(ys) < 2:
        raise ValueError("grid must be at least 2x2 to compute edges")
    x_edges = np.concatenate([[xs[0] - (xs[1]-xs[0])/2], (xs[:-1]+xs[1:])/2, [xs[-1] + (xs[-1]-xs[-2])/2]])
    y_edges = np.concatenate([[ys[0] - (ys[1]-ys[0])/2], (ys[:-1]+ys[1:])/2, [ys[-1] + (ys[-1]-ys[-2])/2]])
    return x_edges, y_edges

# ─────────────── GeoTIFF writers (north-up, no basemap or points) ───────────────
def write_geotiff_multiband(xs, ys, stack, band_names, path, nodata=-9999.0):
    """Write K-band Float32 raster (probability surfaces), north-up."""
    bands, ny, nx = stack.shape
    x_edges, y_edges = grid_edges(xs, ys)
    west, south, east, north = x_edges[0], y_edges[0], x_edges[-1], y_edges[-1]
    transform = from_bounds(west, south, east, north, nx, ny)

    # Flip rows so row 0 = north
    data = np.where(np.isfinite(stack), stack, nodata).astype("float32")
    data = data[:, ::-1, :]

    with rasterio.open(
        path, "w",
        driver="GTiff",
        height=ny, width=nx, count=bands,
        dtype="float32", crs="EPSG:4326", transform=transform,
        nodata=nodata, tiled=True, compress="deflate", predictor=3
    ) as dst:
        for i in range(bands):
            dst.write(data[i], i+1)
            dst.set_band_description(i+1, band_names[i])

def write_geotiff_single(xs, ys, band, path, dtype, nodata, colormap=None, descr=None):
    """Write a single-band raster (float or indexed class), north-up."""
    ny, nx = band.shape
    x_edges, y_edges = grid_edges(xs, ys)
    west, south, east, north = x_edges[0], y_edges[0], x_edges[-1], y_edges[-1]
    transform = from_bounds(west, south, east, north, nx, ny)

    arr = np.where(np.isfinite(band), band, nodata).astype(dtype)
    arr = arr[::-1, :]  # flip rows so row 0 = north

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

# ───────────────────────── image encoding for Folium overlays ─────────────────────────
def scalar_to_rgba_png(z, vmin=0.0, vmax=1.0, cmap_name="viridis", alpha=220):
    """Map a 2D float array to RGBA and return base64 PNG data URL; NaN → fully transparent."""
    cmap = mplcm.get_cmap(cmap_name)
    norm = mplcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    z_mask = np.isfinite(z)
    rgba = np.zeros((z.shape[0], z.shape[1], 4), dtype=np.uint8)
    if z_mask.any():
        mapped = (cmap(norm(np.where(z_mask, z, vmin))) * 255).astype(np.uint8)
        rgba[z_mask] = mapped[z_mask]
        rgba[z_mask, 3] = alpha
        rgba[~z_mask, 3] = 0
    # Flip vertically so row 0 = north
    rgba = rgba[::-1, :, :]
    return "data:image/png;base64," + encode_png(rgba)

def classes_to_rgba_png(classes, class_hex, nodata_val=0, alpha=230):
    """Map a 2D int array of classes to RGBA PNG; nodata_val → transparent."""
    ny, nx = classes.shape
    rgba = np.zeros((ny, nx, 4), dtype=np.uint8)
    for idx, hx in enumerate(class_hex, start=1):
        r = int(hx[1:3], 16); g = int(hx[3:5], 16); b = int(hx[5:7], 16)
        mask = classes == idx
        rgba[mask, 0] = r; rgba[mask, 1] = g; rgba[mask, 2] = b; rgba[mask, 3] = alpha
    rgba[classes == nodata_val, 3] = 0
    # Flip vertically so row 0 = north
    rgba = rgba[::-1, :, :]
    return "data:image/png;base64," + encode_png(rgba)

# ───────────────────────── Folium map builders (HTML only) ─────────────────────────
def _styled_map(basemap_name, center_lat, center_lon):
    # Try requested basemap; fallback if needed (matches your example behavior)
    try:
        m = folium.Map(
            location=[center_lat, center_lon],
            zoom_start=6, tiles=basemap_name,
            control_scale=True, zoom_control=False,
            width="100%", height="100%"
        )
    except Exception:
        m = folium.Map(
            location=[center_lat, center_lon],
            zoom_start=6, tiles="CartoDB Positron",
            control_scale=True, zoom_control=False,
            width="100%", height="100%"
        )

    # Inject the same styling, scale bar look, and north arrow as your example
    m.get_root().header.add_child(Element("""
      <style>
        html,body,#map{height:100%!important;margin:0;}
        .leaflet-control-zoomslider,.leaflet-control-zoom{display:none!important;}
        .leaflet-control-scale{background:rgba(255,255,255,0.9);font-size:10px;padding:2px 6px;}
        #map{border:2px solid #000;}
      </style>
    """))
    m.get_root().html.add_child(Element(
        '<div style="position:absolute;top:8px;right:10px;z-index:999;font-size:1.2em;">↑ N</div>'
    ))
    return m

def _add_overlays(m, overlays):
    for layer in overlays:
        style = layer.get("style", {})
        gdf   = gpd.read_file(layer["path"]).to_crs(epsg=4326)
        folium.GeoJson(
            gdf.__geo_interface__,
            name=Path(layer["path"]).stem,
            style_function=lambda f, s=style: s
        ).add_to(m)

def _add_samples_black_dots(m, samples_df):
    if samples_df is None or samples_df.empty:
        return
    for _, r in samples_df.iterrows():
        folium.CircleMarker(
            location=[float(r["Latitude"]), float(r["Longitude"])],
            radius=3, color="#000000", weight=0.5, fill=True, fill_opacity=1.0
        ).add_to(m)

def _fit_to_bounds(m, bounds):
    west, south, east, north = bounds
    m.fit_bounds([[south, west], [north, east]])

def _base_map_from_bounds(bounds, basemap_name):
    west, south, east, north = bounds
    m = _styled_map(basemap_name, center_lat=(south+north)/2, center_lon=(west+east)/2)
    _fit_to_bounds(m, bounds)
    return m, [[south, west], [north, east]]

def build_folium_continuous_map(xs, ys, stack, layer_names, samples_df, overlays, basemap):
    # Compute exact raster bounds from grid edges → fit view exactly to overlay extent
    x_edges, y_edges = grid_edges(xs, ys)
    west, south, east, north = x_edges[0], y_edges[0], x_edges[-1], y_edges[-1]
    m, bounds = _base_map_from_bounds((west, south, east, north), basemap)

    # One FeatureGroup per K; radio-style toggle via JS below
    for i, name in enumerate(layer_names):
        url = scalar_to_rgba_png(stack[i], vmin=0.0, vmax=1.0, cmap_name="viridis", alpha=220)
        fg = folium.FeatureGroup(name=name, show=(i == 0))
        folium.raster_layers.ImageOverlay(
            image=url, bounds=bounds, opacity=1.0, interactive=False, cross_origin=False
        ).add_to(fg)
        fg.add_to(m)

    _add_overlays(m, overlays)
    _add_samples_black_dots(m, samples_df)
    folium.LayerControl(collapsed=False).add_to(m)

    # Convert overlay checkboxes → radio (exactly like your pie-chart script)
    m.get_root().html.add_child(Element("""
<script>
document.querySelectorAll('.leaflet-control-layers-overlays input').forEach((i)=>{ i.type='radio'; i.name='radios'; });
</script>
"""))
    return m

def build_folium_single_map(xs, ys, image_url, samples_df, overlays, basemap, layer_name):
    x_edges, y_edges = grid_edges(xs, ys)
    west, south, east, north = x_edges[0], y_edges[0], x_edges[-1], y_edges[-1]
    m, bounds = _base_map_from_bounds((west, south, east, north), basemap)
    fg = folium.FeatureGroup(name=layer_name, show=True)
    folium.raster_layers.ImageOverlay(
        image=image_url, bounds=bounds, opacity=1.0, interactive=False, cross_origin=False
    ).add_to(fg)
    fg.add_to(m)
    _add_overlays(m, overlays)
    _add_samples_black_dots(m, samples_df)
    folium.LayerControl(collapsed=False).add_to(m)
    return m

def save_folium_with_template(m, template_path, out_path):
    header = Path(template_path).read_text() if template_path else ""
    html = m._repr_html_()
    Path(out_path).write_text(header + html)
    print(f"✅ HTML written with MultiQC header → {out_path}")

# ───────────────────────── main ─────────────────────────
def main():
    ap = argparse.ArgumentParser()
    # same inputs
    ap.add_argument("--qmat",        required=True, help="ADMIXTURE .Q file (colon-separated)")
    ap.add_argument("--inds",        required=True, help="File with individual IDs (one per line)")
    ap.add_argument("--pops",        required=True, help="File with population/site IDs (one per line)")
    ap.add_argument("--site_coords", required=True, help="TSV with ID, Latitude, Longitude (no header)")
    ap.add_argument("--geo_data_json", help="Optional GeoJSON/shapefile overlay descriptor")
    # MultiQC templates (per output)
    ap.add_argument("--template_disc",   help="MultiQC header template for discrete (argmax) HTML")
    ap.add_argument("--template_div",    help="MultiQC header template for Simpson diversity HTML")
    # Outputs
    ap.add_argument("--out_prefix",    required=True, help="Prefix for all outputs, e.g., outputs/wtd_admix_kriging")
    ap.add_argument("--out_dir",       help="Optional directory for outputs (defaults to prefix parent)")
    ap.add_argument("--table_out",     help="Write site summary TSV (default: <prefix>.tsv)")
    # Display & compute options
    ap.add_argument("--palette",       default="Spectral", help="Palette name for discrete/merged")
    ap.add_argument("--basemap",       default="Esri.NatGeoWorldMap", help="Leaflet/xyzservices tiles name")
    ap.add_argument("--grid_nx",       type=int, default=400)
    ap.add_argument("--grid_ny",       type=int, default=400)
    ap.add_argument("--no_pykrige",    action="store_true", help="Force IDW instead of PyKrige")
    ap.add_argument("--jobs",          type=int, default=0, help="Parallel workers across K (0=auto, 1=serial)")
    # Merged view thresholds
    ap.add_argument("--merge_thresh",  type=float, default=0.5, help="Min winning proportion to color (else gray)")
    ap.add_argument("--merge_margin",  type=float, default=0.0, help="Min (winner - runner_up) to color (else gray)")
    args = ap.parse_args()

    prefix = Path(args.out_prefix)
    outdir  = Path(args.out_dir) if args.out_dir else prefix.parent
    outdir.mkdir(parents=True, exist_ok=True)
    stem = prefix.name

    def OP(suffix, ext):
        return outdir / (f"{stem}{('_' + suffix) if suffix else ''}{ext}")

    # Output paths
    tiff_disc = OP("kriging_discrete",   ".tif")
    tiff_div  = OP("kriging_simpson",    ".tif")
    html_disc = OP("kriging_discrete",   ".html")
    html_div  = OP("kriging_simpson",    ".html")
    tsv_path  = Path(args.table_out) if args.table_out else OP("", ".tsv")

    # Load and aggregate
    df = load_data(args.qmat, args.inds, args.pops)
    cluster_cols = [c for c in df.columns if c.startswith("Cluster")]
    coords       = load_coords(args.site_coords)
    df_sites     = aggregate_sites(df, coords, cluster_cols)
    if df_sites.empty:
        raise ValueError("No valid site coordinates found")

    # TSV
    df_sites[["ID","Latitude","Longitude","count"] + cluster_cols] \
        .to_csv(tsv_path, sep="\t", index=False)
    print(f"✅ Site summary TSV saved to: {tsv_path}")

    overlays = load_overlays(args.geo_data_json) if args.geo_data_json else []

    # Grid
    bounds = compute_bounds(df_sites, overlays)
    xs, ys, xx, yy = build_grid(bounds, args.grid_nx, args.grid_ny)

    # Kriging inputs
    x = df_sites["Longitude"].to_numpy(float)
    y = df_sites["Latitude"].to_numpy(float)
    method = "pykrige" if (not args.no_pykrige and try_pykrige()) else "idw"
    if method == "idw":
        print("ℹ️ PyKrige unavailable/disabled; using IDW fallback.", file=sys.stderr)

    tasks = []
    for name in cluster_cols:
        z = df_sites[name].to_numpy(float)
        tasks.append((name, x, y, z, xs, ys, method))

    # SERIAL kriging across K (no multiprocessing)
    bands = {}
    for t in tasks:
        nm, grid = krige_one_band(t)
        bands[nm] = grid

    # Stack and light renorm
    stack = np.stack([bands[n] for n in cluster_cols], axis=0)  # (K, ny, nx)
    stack = np.clip(stack, 0.0, 1.0)
    sums = np.nansum(stack, axis=0, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        over = sums > 1.0001
        if np.any(over):
            stack[:, over[0]] = stack[:, over[0]] / sums[:, over[0]]

    # GeoTIFFs (rasters only; no basemap/overlays/points)
    # Discrete (argmax)
    K = len(cluster_cols)
    argmax_idx = np.nanargmax(stack, axis=0)  # 0..K-1
    all_nan = np.all(~np.isfinite(stack), axis=0)
    disc = (argmax_idx + 1).astype(np.uint32)   # temp wide int
    disc[all_nan] = 0

    colors = palette(args.palette, K)

    def _hex_to_rgb(h):
        return tuple(int(h[i:i+2], 16) for i in (1,3,5))

    if K <= 255:
        disc8 = disc.astype("uint8")
        cmap = {i: _hex_to_rgb(colors[i-1]) for i in range(1, K+1)}  # RGB, no alpha
        write_geotiff_single(xs, ys, disc8, str(tiff_disc),
                             dtype="uint8", nodata=0, colormap=cmap, descr="Kmax (1..K)")
    else:
        disc16 = disc.astype("uint16")
        write_geotiff_single(xs, ys, disc16, str(tiff_disc),
                             dtype="uint16", nodata=0, colormap=None, descr="Kmax (1..K)")
    print(f"✅ GeoTIFF (discrete argmax): {tiff_disc}")

    # Simpson diversity
    simpson = simpson_diversity(stack)
    write_geotiff_single(xs, ys, simpson, str(tiff_div), dtype="float32", nodata=-9999.0, descr="Simpson diversity")
    print(f"✅ GeoTIFF (Simpson diversity): {tiff_div}")

    # Sample points dataframe for HTML (black dots on top, like original)
    samples_df = df_sites[["Longitude","Latitude","ID","count"] + cluster_cols].copy()

    # HTML (Folium) — basemap under, overlays (geo_data_json) in z-order, points on top
    # Discrete single
    url_disc = classes_to_rgba_png(disc.astype(int), class_hex=colors, nodata_val=0, alpha=230)
    m_disc = build_folium_single_map(xs, ys, url_disc, samples_df, overlays, args.basemap, "Discrete (argmax)")
    save_folium_with_template(m_disc, args.template_disc, html_disc)

    # Simpson single
    url_div  = scalar_to_rgba_png(simpson, vmin=np.nanmin(simpson), vmax=np.nanmax(simpson),
                                  cmap_name="viridis", alpha=220)
    m_div = build_folium_single_map(xs, ys, url_div, samples_df, overlays, args.basemap, "Simpson diversity")
    save_folium_with_template(m_div, args.template_div, html_div)

    print("🎉 All rasters and MultiQC-ready Folium HTMLs written.")

if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=UserWarning)
    main()
