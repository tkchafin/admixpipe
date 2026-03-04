#!/usr/bin/env python3
"""
plot_evaladmix.py

Generate evalAdmix-style "above/below" residual-correlation heatmaps in Plotly.

Outputs:
  1) A single "all K" HTML file with a dropdown selector (optional)
     - Implementation uses a lightweight HTML+JS wrapper that swaps pre-rendered
       Plotly <div> blocks (fast, robust, no trace/shape visibility gymnastics).
  2) A single "best K" HTML file (optional; requires --best_k)

Both outputs prepend a MultiQC-style metadata HTML comment block extracted
from a provided template file.

Notes:
  - Keeps your exact ordering logic (orderInds(pop=fam_col2, q=Q)).
  - Averages corres over replicates per K with plain sum (propagates NaN like R).
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio


# ----------------------------
# MultiQC-style template helpers
# ----------------------------

def parse_template(template_file: str) -> dict:
    """Extract metadata from the first HTML comment block <!-- ... -->."""
    with open(template_file, "r") as f:
        content = f.read()
    match = re.search(r"<!--(.*?)-->", content, re.DOTALL)
    if not match:
        raise ValueError(f"No metadata block found in template: {template_file}")
    raw = match.group(1)
    meta: dict = {}
    for line in raw.strip().splitlines():
        if ":" in line:
            key, val = line.strip().split(":", 1)
            meta[key.strip()] = val.strip().strip('"\'')
    return meta


def build_comment(meta: dict) -> str:
    lines = ["<!--"]
    for k, v in meta.items():
        lines.append(f'{k}: "{v}"')
    lines.append("-->")
    return "\n".join(lines)


def header_from_template(
    template_path: str,
    override_id: Optional[str] = None,
    override_title: Optional[str] = None,
    override_section_name: Optional[str] = None,
    override_description: Optional[str] = None,
) -> str:
    meta = parse_template(template_path)
    if override_id:
        meta["id"] = override_id
    if override_title:
        meta["title"] = override_title
    if override_section_name:
        meta["section_name"] = override_section_name
    if override_description:
        meta["description"] = override_description
    return build_comment(meta)


# ----------------------------
# IO helpers
# ----------------------------

def read_corres_matrix(path: str) -> np.ndarray:
    """Read evalAdmix *.corres into NxN float matrix (R read.table -> as.matrix)."""
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    mat = df.to_numpy(dtype=float)
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise ValueError(f"corres matrix must be square; got {mat.shape} from {path}")
    return mat


def read_fam_pop(prefix: str) -> Tuple[List[str], List[str], List[str]]:
    """
    Read <prefix>.fam and return (fid, iid, pop) where pop == column 2 (IID),
    EXACTLY like the wrapper uses pop.rx(True,2).
    """
    famf = f"{prefix}.fam"
    if not os.path.isfile(famf):
        raise FileNotFoundError(f"Missing FAM: {famf}")
    fam = pd.read_csv(famf, sep=r"\s+", header=None, engine="python")
    if fam.shape[1] < 2:
        raise ValueError(f"FAM has <2 columns: {famf}")
    fid = fam.iloc[:, 0].astype(str).tolist()
    iid = fam.iloc[:, 1].astype(str).tolist()
    pop = fam.iloc[:, 1].astype(str).tolist()  # wrapper behavior
    return fid, iid, pop


_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")

def read_q_matrix(path: str) -> np.ndarray:
    """
    Read either:
      - CLUMPP IndFile.output.K lines like: '1  1  (0)  1  :  0.7582 0.2418'
        -> parse floats AFTER ':' per line
      - ADMIXTURE .Q (pure numeric whitespace table)
    Returns: N x K float numpy array.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing Q file: {path}")

    # Fast path: try pure numeric table (ADMIXTURE .Q)
    try:
        df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
        q = df.to_numpy(dtype=float)
        return q
    except Exception:
        pass

    # CLUMPP IndFile path: parse line by line
    rows: List[List[float]] = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if ":" in line:
                rhs = line.split(":", 1)[1]
                floats = _FLOAT_RE.findall(rhs)
                if floats:
                    rows.append([float(x) for x in floats])
                continue

            floats = _FLOAT_RE.findall(line)
            if floats:
                rows.append([float(x) for x in floats])

    if not rows:
        raise ValueError(f"Could not parse any Q rows from: {path}")

    lens = {len(r) for r in rows}
    if len(lens) != 1:
        raise ValueError(
            f"Inconsistent number of Q columns parsed from {path}: lengths={sorted(lens)}. "
            f"Check file format."
        )

    return np.asarray(rows, dtype=float)


# ----------------------------
# ports of R logic from visFunctions.R in evaladmix
# ----------------------------

def unique_preserve_order(x: Sequence[str]) -> List[str]:
    """Match R unique(): first occurrence order."""
    return list(pd.unique(np.asarray(x, dtype=str)))


def order_inds_exact(q: Optional[np.ndarray] = None,
                     pop: Optional[Sequence[str]] = None,
                     popord: Optional[Sequence[str]] = None) -> np.ndarray:
    """
    Port of R orderInds().
    Returns 0-based indices.
    """
    if pop is not None:
        pop_arr = np.asarray(pop, dtype=str)
        if popord is None:
            popord = unique_preserve_order(pop_arr)
        else:
            popord = list(popord)

        if q is not None:
            q = np.asarray(q, dtype=float)
            ord_all: List[int] = []

            for x in popord:
                idx = np.where(pop_arr == x)[0]
                if idx.size == 0:
                    continue
                means = np.nanmean(q[idx, :], axis=0)
                main_k = int(np.argmax(means))  # which.max (first max)
                within = np.argsort(q[idx, main_k], kind="mergesort")  # stable like R
                ord_all.extend(idx[within].tolist())

            return np.asarray(ord_all, dtype=int)

        # pop but no q
        ord_all = []
        for x in popord:
            ord_all.extend(np.where(pop_arr == x)[0].tolist())
        return np.asarray(ord_all, dtype=int)

    if pop is None and q is not None:
        q = np.asarray(q, dtype=float)
        main_k = np.argmax(q, axis=1)
        main_q = q[np.arange(q.shape[0]), main_k]
        # order(main_k, main_q)
        return np.lexsort((main_q, main_k)).astype(int)

    raise ValueError("Need at least an argument to order (pop or q).")


@dataclass
class PlotCorResResult:
    z_plot: np.ndarray           # reordered, upper replaced, diag=10
    z_raw: np.ndarray            # reordered original (individual) matrix
    mean_cors: pd.DataFrame      # pop x pop means
    pop_ord: np.ndarray          # ordered pop labels
    ord0: np.ndarray             # ordering indices (0-based)
    zmin: float
    zmax: float


def plotcorres_prepare_exact(cor_mat: np.ndarray,
                             pop: Optional[Sequence[str]] = None,
                             ord0: Optional[Sequence[int]] = None,
                             min_z: float = np.nan,
                             max_z: float = np.nan) -> PlotCorResResult:
    """Port of R plotCorRes up to image()."""
    cor_mat = np.asarray(cor_mat, dtype=float)
    if cor_mat.ndim != 2 or cor_mat.shape[0] != cor_mat.shape[1]:
        raise ValueError("cor_mat must be square NxN")
    N = cor_mat.shape[0]

    if ord0 is None and pop is not None:
        pop_arr0 = np.asarray(pop, dtype=str)
        ord0 = np.argsort(pop_arr0, kind="mergesort")  # R order(pop)
    elif ord0 is None:
        ord0 = np.arange(N, dtype=int)
    else:
        ord0 = np.asarray(ord0, dtype=int)

    if pop is None:
        pop_arr = np.array([" "] * N, dtype=str)
    else:
        pop_arr = np.asarray(pop, dtype=str)

    pop_ord = pop_arr[ord0]
    cor_ord = cor_mat[np.ix_(ord0, ord0)]
    z_raw = cor_ord.copy()

    uniq_pop = unique_preserve_order(pop_ord)
    mean_cors = pd.DataFrame(index=uniq_pop, columns=uniq_pop, dtype=float)

    for i1, p1 in enumerate(uniq_pop):
        idx1 = np.where(pop_ord == p1)[0]
        for i2, p2 in enumerate(uniq_pop):
            idx2 = np.where(pop_ord == p2)[0]
            block = cor_ord[np.ix_(idx1, idx2)]
            mean_cors.iat[i1, i2] = float(np.nanmean(block))  # mean(., na.rm=TRUE)

    z = cor_ord.copy()
    for i1 in range(N - 1):
        for i2 in range(i1 + 1, N):
            z[i1, i2] = mean_cors.loc[pop_ord[i2], pop_ord[i1]]

    z_lims = np.array([min_z, max_z], dtype=float)

    if np.all(np.isnan(z_lims)):
        finite = z[np.isfinite(z)]
        m = float(np.max(np.abs(finite))) if finite.size else 0.0
        z_lims = np.array([-m, m], dtype=float)

    if np.any(np.isnan(z_lims)):
        v = float(z_lims[~np.isnan(z_lims)][0])
        z_lims = np.array([-v, v], dtype=float)

    zmin2, zmax2 = float(z_lims[0]), float(z_lims[1])

    np.fill_diagonal(z, 10.0)

    return PlotCorResResult(
        z_plot=z,
        z_raw=z_raw,
        mean_cors=mean_cors,
        pop_ord=pop_ord,
        ord0=ord0,
        zmin=zmin2,
        zmax=zmax2,
    )


# ----------------------------
# Plotly figure building
# ----------------------------

def build_plotly_figure(res: PlotCorResResult,
                        title: str,
                        fid: Optional[Sequence[str]] = None,
                        iid: Optional[Sequence[str]] = None,
                        color_palette: Sequence[str] = ("#001260", "#EAEDE9", "#601200"),
                        show_pop_boundaries: bool = True) -> go.Figure:
    """
    Plotly heatmap matching base R image(t(mat)) orientation:
      - Use transpose for z
      - Reverse y axis
    Hover indicates whether cell is individual (lower) or pop-mean (upper).
    """
    z = res.z_plot.copy()
    z_raw = res.z_raw
    pop = res.pop_ord
    N = z.shape[0]

    fid_ord = np.asarray(fid, dtype=str)[res.ord0] if fid is not None else None
    iid_ord = np.asarray(iid, dtype=str)[res.ord0] if iid is not None else None

    hover = np.empty((N, N), dtype=object)
    for i in range(N):
        for j in range(N):
            if i == j:
                hover[i, j] = "Diagonal"
                continue
            if j > i:
                val = z[i, j]
                hover[i, j] = (
                    f"Type: pop-mean (upper triangle)<br>"
                    f"Row pop: {pop[i]}<br>Col pop: {pop[j]}<br>"
                    f"Value: {val:.6g}"
                )
            else:
                val = z_raw[i, j]
                extra = ""
                if fid_ord is not None and iid_ord is not None:
                    extra = (
                        f"<br>Row: {fid_ord[i]}/{iid_ord[i]}"
                        f"<br>Col: {fid_ord[j]}/{iid_ord[j]}"
                    )
                hover[i, j] = (
                    f"Type: individual (lower triangle)<br>"
                    f"Row pop: {pop[i]}<br>Col pop: {pop[j]}<br>"
                    f"Value: {val:.6g}{extra}"
                )

    # Mask diagonal for heatmap; overlay it as black points
    z_heat = z.copy()
    np.fill_diagonal(z_heat, np.nan)

    # Match R image(t(z)):
    z_plot = z_heat.T
    hover_plot = hover.T

    colorscale = [
        (0.0, color_palette[0]),
        (0.5, color_palette[1]),
        (1.0, color_palette[2]),
    ]

    fig = go.Figure()

    fig.add_trace(
        go.Heatmap(
            z=z_plot,
            zmin=res.zmin,
            zmax=res.zmax,
            colorscale=colorscale,
            hovertext=hover_plot,
            hoverinfo="text",
            colorbar=dict(title="corr"),
        )
    )

    # Diagonal in black
    fig.add_trace(
        go.Scatter(
            x=np.arange(N),
            y=np.arange(N),
            mode="markers",
            marker=dict(size=5, color="black"),
            hoverinfo="skip",
            showlegend=False,
        )
    )

    if show_pop_boundaries:
        uniq = unique_preserve_order(pop)
        counts = [int(np.sum(pop == u)) for u in uniq]
        cuts = np.cumsum(counts)
        for c in cuts[:-1]:
            pos = c - 0.5
            fig.add_shape(type="line", x0=pos, x1=pos, y0=-0.5, y1=N - 0.5,
                          line=dict(color="black", width=1))
            fig.add_shape(type="line", x0=-0.5, x1=N - 0.5, y0=pos, y1=pos,
                          line=dict(color="black", width=1))

    fig.update_layout(
        title=title,
        width=900,
        height=900,
        xaxis=dict(showticklabels=False, range=[-0.5, N - 0.5]),
        yaxis=dict(showticklabels=False, range=[N - 0.5, -0.5]),  # reverse like R image
        margin=dict(l=40, r=40, t=60, b=40),
    )

    return fig


# ----------------------------
# Discovery: group corres by K from filename
# ----------------------------

_CORRES_RE = re.compile(r"^(?P<prefix>.+)\.(?P<K>\d+)_(?P<rep>\d+)\.corres$")


def find_corres_by_k(prefix: str, workdir: str) -> Dict[str, List[str]]:
    """Discover files matching <prefix>.<K>_<rep>.corres in workdir and group by K."""
    pattern = os.path.join(workdir, f"{prefix}.*_*.corres")
    files = [os.path.basename(p) for p in glob.glob(pattern)]
    out: Dict[str, List[str]] = {}
    for fn in sorted(files):
        m = _CORRES_RE.match(fn)
        if not m:
            continue
        if m.group("prefix") != prefix:
            continue
        k = m.group("K")
        out.setdefault(k, []).append(os.path.join(workdir, fn))
    return out


def average_corres(mats: List[np.ndarray]) -> np.ndarray:
    """
    Wrapper behavior: Reduce('+', matrixList) / nRuns.
    Plain sum propagates NaN like R's NA propagation.
    """
    stack = np.stack(mats, axis=0)
    return np.sum(stack, axis=0) / float(stack.shape[0])


def choose_q_for_k(prefix: str,
                   k: str,
                   workdir: str,
                   qfilepaths_json: Optional[str] = None) -> Tuple[str, np.ndarray]:
    """
    Prefer CLUMPP consensus Q if qfilePaths.json is present.
    Otherwise fall back to first replicate .Q found: <prefix>.<k>_<rep>.Q
    """
    if qfilepaths_json and os.path.isfile(qfilepaths_json):
        with open(qfilepaths_json) as fh:
            mp = json.load(fh)
        if k in mp and os.path.isfile(mp[k]):
            return mp[k], read_q_matrix(mp[k])

    cand = sorted(glob.glob(os.path.join(workdir, f"{prefix}.{k}_*.Q")))
    if not cand:
        raise FileNotFoundError(
            f"No Q found for K={k}. Provide qfilePaths.json or ensure {prefix}.{k}_rep.Q exists."
        )
    return cand[0], read_q_matrix(cand[0])


# ----------------------------
# HTML writers (improved)
# ----------------------------

def _sanitize_id(s: str) -> str:
    s2 = re.sub(r"[^a-zA-Z0-9_\-]+", "_", s)
    if re.match(r"^[0-9]", s2):
        s2 = "id_" + s2
    return s2


def write_allk_dropdown_html(figs_by_k: Dict[int, go.Figure],
                             out_html: str,
                             header_comment: str,
                             page_title: str = "evalAdmix residual correlations (all K)") -> None:
    """
    Robust all-K HTML: we render one Plotly div per K and use a dropdown + tiny JS
    to show/hide the corresponding div. No merging traces/shapes, no Plotly layout hacks.

    Uses Plotly CDN *once* via the first rendered div and removes it from others.
    """
    if not figs_by_k:
        raise ValueError("No figures provided for all-K output.")

    ks = sorted(figs_by_k.keys())
    div_ids = {k: _sanitize_id(f"evaladmix_k{k}") for k in ks}

    # Render all figures as divs. Include plotlyjs only for first, omit for others.
    div_blocks: List[str] = []
    for i, k in enumerate(ks):
        include_js = "cdn" if i == 0 else False
        div = pio.to_html(
            figs_by_k[k],
            full_html=False,
            include_plotlyjs=include_js,
            config={"responsive": True},
        )
        # Wrap div so we can toggle visibility
        style = "" if i == 0 else "display:none;"
        wrapped = (
            f"<div id='{div_ids[k]}' class='kpanel' style='{style}'>\n"
            f"{div}\n"
            f"</div>\n"
        )
        div_blocks.append(wrapped)

    # Simple dropdown + JS
    options = "\n".join([f"<option value='{k}'>K={k}</option>" for k in ks])
    first_k = ks[0]

    js = f"""
<script>
(function() {{
  const panels = document.querySelectorAll('.kpanel');
  function showK(k) {{
    panels.forEach(p => p.style.display = 'none');
    const el = document.getElementById('{div_ids[first_k]}'.replace('{first_k}', k));
  }}
}})();
</script>
"""

    # The above attempted replacement is awkward in JS; do it cleanly:
    js = f"""
<script>
(function() {{
  const idByK = {{"{str(ks[0])}": "{div_ids[ks[0]]}"{"".join([f', "{k}": "{div_ids[k]}"' for k in ks[1:]])}}};
  const panels = document.querySelectorAll('.kpanel');

  function showK(k) {{
    panels.forEach(p => p.style.display = 'none');
    const id = idByK[String(k)];
    if (id) {{
      const el = document.getElementById(id);
      if (el) el.style.display = '';
    }}
  }}

  const sel = document.getElementById('kselector');
  if (sel) {{
    sel.addEventListener('change', (e) => showK(e.target.value));
  }}
  // initial
  showK({first_k});
}})();
</script>
"""

    controls = f"""
<div style="display:flex; align-items:center; gap:10px; margin: 8px 0 14px 0;">
  <label for="kselector" style="font-family: sans-serif;">Select K:</label>
  <select id="kselector" style="font-size: 14px; padding: 4px 6px;">
    {options}
  </select>
</div>
"""

    html = (
        header_comment
        + "\n"
        + "<!DOCTYPE html>\n<html>\n<head>\n"
        + f"<meta charset='utf-8'>\n<title>{page_title}</title>\n"
        + "<meta name='viewport' content='width=device-width, initial-scale=1'>\n"
        + "</head>\n<body>\n"
        + controls
        + "\n".join(div_blocks)
        + js
        + "\n</body>\n</html>\n"
    )

    os.makedirs(os.path.dirname(out_html) or ".", exist_ok=True)
    with open(out_html, "w") as f:
        f.write(html)


def write_singlek_html(fig: go.Figure,
                       out_html: str,
                       header_comment: str,
                       page_title: str = "evalAdmix residual correlations") -> None:
    div = pio.to_html(fig, full_html=False, include_plotlyjs="cdn", config={"responsive": True})
    html = (
        header_comment
        + "\n"
        + "<!DOCTYPE html>\n<html>\n<head>\n"
        + f"<meta charset='utf-8'>\n<title>{page_title}</title>\n"
        + "<meta name='viewport' content='width=device-width, initial-scale=1'>\n"
        + "</head>\n<body>\n"
        + div
        + "\n</body>\n</html>\n"
    )
    os.makedirs(os.path.dirname(out_html) or ".", exist_ok=True)
    with open(out_html, "w") as f:
        f.write(html)


# ----------------------------
# Main run logic
# ----------------------------

def build_all_figs(prefix: str,
                   workdir: str,
                   qfilepaths_json: Optional[str],
                   min_z: float,
                   max_z: float) -> Tuple[Dict[int, go.Figure], List[str]]:
    by_k = find_corres_by_k(prefix, workdir)
    if not by_k:
        raise SystemExit(f"No corres files found for prefix={prefix} in {workdir}")

    fid, iid, pop = read_fam_pop(prefix)

    figs_by_k: Dict[int, go.Figure] = {}
    loglines: List[str] = []

    for k_str in sorted(by_k.keys(), key=lambda x: int(x)):
        k = int(k_str)
        corres_files = by_k[k_str]
        mats = [read_corres_matrix(p) for p in corres_files]
        cor_mean = average_corres(mats)

        q_path, q = choose_q_for_k(prefix, k_str, workdir, qfilepaths_json=qfilepaths_json)

        ord0 = order_inds_exact(q=q, pop=pop)
        res = plotcorres_prepare_exact(cor_mean, pop=pop, ord0=ord0, min_z=min_z, max_z=max_z)

        title = f"K={k}"
        fig = build_plotly_figure(res, title=title, fid=fid, iid=iid)

        figs_by_k[k] = fig
        loglines.append(f"[K={k}] averaged {len(corres_files)} corres; Q={os.path.basename(q_path)}")

    return figs_by_k, loglines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True, help="PLINK prefix (expects <prefix>.fam)")
    ap.add_argument("--workdir", default=".", help="Directory containing *.corres and *.Q")
    ap.add_argument("--outdir", default=".", help="Output directory (for HTML files)")
    ap.add_argument("--qfilePaths", default=None,
                    help="Optional qfilePaths.json for CLUMPP consensus Q per K (preferred).")
    ap.add_argument("--min_z", type=float, default=-0.1, help="min_z as in plotCorRes")
    ap.add_argument("--max_z", type=float, default=0.1, help="max_z as in plotCorRes")

    # All-K output
    ap.add_argument("--allk_html", default=None,
                    help="If set, write a single HTML containing all K with a dropdown selector (basename or path).")
    ap.add_argument("--template_allk", default=None,
                    help="Template HTML containing MultiQC metadata comment block for all-K output.")
    ap.add_argument("--allk_id", default=None, help="Override metadata 'id' for all-K output.")
    ap.add_argument("--allk_title", default=None, help="Override metadata 'title' for all-K output.")
    ap.add_argument("--allk_section_name", default=None, help="Override metadata 'section_name' for all-K output.")
    ap.add_argument("--allk_description", default=None, help="Override metadata 'description' for all-K output.")

    # Best-K output
    ap.add_argument("--best_k", type=int, default=None,
                    help="If set, also write a single HTML for this K only.")
    ap.add_argument("--bestk_html", default="evaladmix_bestk.html",
                    help="Output HTML filename for best-K plot (fixed default).")
    ap.add_argument("--template_bestk", default=None,
                    help="Template HTML containing MultiQC metadata comment block for best-K output.")
    ap.add_argument("--bestk_id", default=None, help="Override metadata 'id' for best-K output.")
    ap.add_argument("--bestk_title", default=None, help="Override metadata 'title' for best-K output.")
    ap.add_argument("--bestk_section_name", default=None, help="Override metadata 'section_name' for best-K output.")
    ap.add_argument("--bestk_description", default=None, help="Override metadata 'description' for best-K output.")

    args = ap.parse_args()

    figs_by_k, loglines = build_all_figs(
        prefix=args.prefix,
        workdir=args.workdir,
        qfilepaths_json=args.qfilePaths,
        min_z=args.min_z,
        max_z=args.max_z,
    )

    # All-K HTML
    if args.allk_html:
        if not args.template_allk:
            raise SystemExit("--allk_html requires --template_allk")

        header_allk = header_from_template(
            args.template_allk,
            override_id=args.allk_id,
            override_title=args.allk_title,
            override_section_name=args.allk_section_name,
            override_description=args.allk_description,
        )

        out_allk = args.allk_html
        if not os.path.isabs(out_allk):
            out_allk = os.path.join(args.outdir, out_allk)

        write_allk_dropdown_html(
            figs_by_k=figs_by_k,
            out_html=out_allk,
            header_comment=header_allk,
            page_title="evalAdmix residual correlations (all K)",
        )
        print(f"[allK] wrote: {out_allk}")

    # Best-K HTML
    if args.best_k is not None:
        k = int(args.best_k)
        if k not in figs_by_k:
            have = ", ".join(str(x) for x in sorted(figs_by_k.keys()))
            raise SystemExit(f"--best_k={k} not found. Available K: {have}")

        if not args.template_bestk:
            raise SystemExit("--best_k requires --template_bestk")

        header_bestk = header_from_template(
            args.template_bestk,
            override_id=args.bestk_id,
            override_title=args.bestk_title,
            override_section_name=args.bestk_section_name,
            override_description=args.bestk_description,
        )

        out_bestk = args.bestk_html
        if not os.path.isabs(out_bestk):
            out_bestk = os.path.join(args.outdir, out_bestk)

        write_singlek_html(
            fig=figs_by_k[k],
            out_html=out_bestk,
            header_comment=header_bestk,
            page_title=f"evalAdmix residual correlations (K={k})",
        )
        print(f"[bestK] wrote: {out_bestk}")

    # Always print build info
    for ln in loglines:
        print(ln)


if __name__ == "__main__":
    main()
