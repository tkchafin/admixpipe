#!/usr/bin/env python3
"""
evanno.py: Compute Evanno et al. (2005) "Delta K" statistics from STRUCTURE-like log-likelihood output.

Reads a tab-delimited table with columns:
    K    Mean    StDev    [other columns...]

Outputs a table with, for each K:
    K           (as in input)
    Mean        L(K)
    StDev       σ(K)
    Lprime      L'(K) = L(K) - L(K-1)
    sd_Lprime   √[σ(K)^2 + σ(K-1)^2]
    Lpp         L''(K) = L(K+1) - 2 L(K) + L(K-1)
    sd_Lpp      √[σ(K+1)^2 + 4 σ(K)^2 + σ(K-1)^2]
    DeltaK      |Lpp| / σ(K)

Note: for K=1 or K=max, some fields will be NaN.
"""
import sys
import argparse
import pandas as pd
import numpy as np

def compute_evanno(df):
    # Ensure sorted by K
    df = df.sort_values('K').reset_index(drop=True)

    # Extract L(K) and σ(K)
    L = df['Mean'].astype(float)
    sd = df['StDev'].astype(float)

    # First-order change L'(K) and its propagated error
    Lprime = L - L.shift(1)
    sd_Lprime = np.sqrt(sd.pow(2) + sd.shift(1).pow(2))

    # Second-order change L''(K) and its propagated error
    Lpp = L.shift(-1) - 2 * L + L.shift(1)
    sd_Lpp = np.sqrt(
        sd.shift(-1).pow(2)
        + 4 * sd.pow(2)
        + sd.shift(1).pow(2)
    )

    # Delta K = |L''(K)| / σ(K)
    DeltaK = Lpp.abs() / sd

    # Assemble output
    out = pd.DataFrame({
        'K': df['K'],
        'Mean': L,
        'StDev': sd,
        "Lprime": Lprime,
        "sd_Lprime": sd_Lprime,
        "Lpp": Lpp,
        "sd_Lpp": sd_Lpp,
        "DeltaK": DeltaK
    })
    return out


def main():
    p = argparse.ArgumentParser(
        description="Compute Evanno methods metrics (L\', L\'', 'ΔK') for ADMIXTURE output"
    )
    p.add_argument(
        'infile',
        help='Input file: tab-delimited, must contain columns K, Mean, StDev'
    )
    p.add_argument(
        '-o', '--outfile',
        help='Output TSV file (default: stdout)',
        default=None
    )
    args = p.parse_args()

    # Read input
    try:
        df = pd.read_csv(args.infile, sep=r'\s+|\t', engine='python')
    except Exception:
        sys.exit(f"Error reading input file {args.infile}")

    # Check required columns
    for col in ('K', 'Mean', 'StDev'):
        if col not in df.columns:
            sys.exit(f"Input file must contain column '{col}'")

    # Compute metrics
    result = compute_evanno(df)

    # Write output
    if args.outfile:
        result.to_csv(args.outfile, sep='\t', index=False, float_format='%.6f')
    else:
        result.to_csv(sys.stdout, sep='\t', index=False, float_format='%.6f')

if __name__ == '__main__':
    main()
