process BESTK {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/btmartin721/snpio:1.3.21'

    input:
        tuple val(meta), path(cv_file)
        tuple val(meta2), path(evanno_metrics)
        tuple val(meta3), path(best_results)

    output:
        tuple val(meta), path('bestK.txt')                , emit: bestK_file
        tuple val(meta), path('best_clumpp_indfile.out') , emit: bestK_clumpp
        path "versions.yml"                               , emit: versions

    script:
    """
    python <<'PY'
    import csv
    import math

    method = "${params.bestk_method}".strip().lower()
    cv_file = "${cv_file}"
    evanno_file = "${evanno_metrics}"

    def read_tsv(path):
        with open(path, "r", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\\t")
            return [row for row in reader]

    def to_float(x):
        if x is None:
            return None
        x = str(x).strip()
        if x == "" or x.lower() in {"na", "nan", "none"}:
            return None
        return float(x)

    def elbow_k(points):
        '''
        points: list of (k, y), sorted by k
        returns K at maximum perpendicular distance from the line
        joining first and last point; ties go to smaller K
        '''
        if len(points) == 0:
            raise ValueError("No valid points available for elbow detection")
        if len(points) == 1:
            return points[0][0]
        if len(points) == 2:
            return points[0][0]

        # prefer interior points only for elbow detection
        interior = points[1:-1]
        if not interior:
            return points[0][0]

        x1, y1 = points[0]
        x2, y2 = points[-1]

        denom = math.hypot(y2 - y1, x2 - x1)
        if denom == 0:
            return points[0][0]

        best_k = None
        best_dist = -1.0

        for k, y in interior:
            dist = abs((y2 - y1) * k - (x2 - x1) * y + x2 * y1 - y2 * x1) / denom
            if dist > best_dist or (math.isclose(dist, best_dist) and (best_k is None or k < best_k)):
                best_dist = dist
                best_k = k

        return best_k

    def choose_min_metric(path, metric_col):
        rows = read_tsv(path)
        best_k = None
        best_val = None

        for row in rows:
            k = to_float(row.get("K"))
            val = to_float(row.get(metric_col))

            if k is None or val is None:
                continue

            k = int(k)
            if k <= 1:
                continue

            if best_val is None or val < best_val or (math.isclose(val, best_val) and k < best_k):
                best_val = val
                best_k = k

        if best_k is None:
            raise ValueError(f"Could not determine best K from metric '{metric_col}' in file: {path}")
        return best_k

    def choose_max_metric(path, metric_col):
        rows = read_tsv(path)
        best_k = None
        best_val = None

        for row in rows:
            k = to_float(row.get("K"))
            val = to_float(row.get(metric_col))

            if k is None or val is None:
                continue

            k = int(k)
            if k <= 1:
                continue

            if best_val is None or val > best_val or (math.isclose(val, best_val) and k < best_k):
                best_val = val
                best_k = k

        if best_k is None:
            raise ValueError(f"Could not determine best K from metric '{metric_col}' in file: {path}")
        return best_k

    def choose_elbow(path, metric_col, use_abs=False):
        rows = read_tsv(path)
        points = []

        for row in rows:
            k = to_float(row.get("K"))
            y = to_float(row.get(metric_col))

            if k is None or y is None:
                continue

            k = int(k)
            if k <= 1:
                continue

            if use_abs:
                y = abs(y)

            points.append((k, y))

        points.sort(key=lambda x: x[0])

        if not points:
            raise ValueError(f"No valid points found for metric '{metric_col}' in {path}")

        return elbow_k(points)

    if method == "cv":
        best_k = choose_min_metric(cv_file, "Mean")
    elif method == "evanno":
        best_k = choose_max_metric(evanno_file, "DeltaK")
    elif method == "lnl":
        best_k = choose_elbow(evanno_file, "Mean", use_abs=False)
    elif method == "l1":
        best_k = choose_elbow(evanno_file, "Lprime", use_abs=False)
    elif method == "l2":
        best_k = choose_elbow(evanno_file, "Lpp", use_abs=True)
    else:
        raise ValueError(f"Unsupported bestk_method: {method}")

    with open("bestK.txt", "w") as out:
        out.write(f"{best_k}\\n")
    PY

    K=\$(cat bestK.txt)

    MATCHED_FILE=\$(find "${best_results}/" -type f -name "ClumppIndFile.output.\$K" -print -quit)

    if [ ! -f "\$MATCHED_FILE" ]; then
        echo "ERROR: Could not find ClumppIndFile.output.\$K in ${best_results}" >&2
        exit 1
    fi

    MATCHED_FILE=\$(realpath "\$MATCHED_FILE")
    ln -s "\$MATCHED_FILE" best_clumpp_indfile.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
