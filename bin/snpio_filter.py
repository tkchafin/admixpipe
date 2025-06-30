#!/usr/bin/env python3
import argparse
import os
from snpio import NRemover2, VCFReader


def get_prefix_from_vcf_path(vcf_path):
    basename = os.path.basename(vcf_path)
    for ext in [".vcf.gz", ".vcf"]:
        if basename.endswith(ext):
            return basename[: -len(ext)]
    return basename  # fallback if no known extension


def main():
    parser = argparse.ArgumentParser(
        description="Run SNPio to filter individuals and generate missingness reports"
    )
    parser.add_argument("--vcf", required=True, help="Path to the VCF file.")
    parser.add_argument(
        "--popmap", required=True, help="Path to the population map file."
    )
    parser.add_argument(
        "--ind_cov",
        type=float,
        default=0.75,
        help="Maximum allowed missingness per individual (default: 0.75)",
    )
    parser.add_argument(
        "--min_maf",
        type=float,
        default=0.05,
        help="Maximum MAF to retain a SNP (default: 0.05)",
    )
    parser.add_argument(
        "--snp_cov",
        type=float,
        default=0.9,
        help="Maximum missing data to retain a SNP (default: 0.9)",
    )
    args = parser.parse_args()

    # extract prefix from VCF filename
    prefix = get_prefix_from_vcf_path(args.vcf)

    # read data
    gd = VCFReader(
        filename=args.vcf,
        popmapfile=args.popmap,
        force_popmap=True,
        verbose=True,
        plot_format="png",
        plot_fontsize=8,
        plot_dpi=300,
        prefix=prefix,
    )

    # generate missingness reports
    gd.missingness_reports()

    # Filter VCF
    nrm = NRemover2(gd)
    gd_filt = (
        nrm.filter_missing_sample(args.ind_cov)
        .filter_monomorphic(exclude_heterozygous=False)
        .filter_missing(args.snp_cov)
        .filter_maf(args.min_maf)
        .resolve()
    )
    nrm.plot_sankey_filtering_report()

    # Write the filtered VCF using the modified prefix
    output_vcf = f"{prefix}.filter.vcf"
    gd_filt.write_vcf(output_vcf)


if __name__ == "__main__":
    main()
