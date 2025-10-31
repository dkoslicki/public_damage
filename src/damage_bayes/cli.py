"""Command line interface for the Damage Bayes toolkit."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from . import __version__
from .counting import write_frequency_tables
from .modeling import analyze_frequency_tables


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Estimate ancient DNA damage patterns using a Bayesian model.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    count_parser = subparsers.add_parser(
        "count",
        help="Generate nucleotide frequency tables from FASTQ input.",
    )
    count_parser.add_argument(
        "fastq",
        nargs="+",
        help="Path(s) to FASTQ files to process.",
    )
    count_parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory where frequency tables will be written. Defaults to the current directory.",
    )
    count_parser.add_argument(
        "--n-bp",
        type=int,
        default=25,
        help="Number of terminal bases to analyse. Default: 25.",
    )
    count_parser.add_argument(
        "--trim",
        type=int,
        default=1,
        help="Number of bases to trim before analysis. Default: 1.",
    )

    analyze_parser = subparsers.add_parser(
        "analyze",
        help="Run the Bayesian model on paired frequency tables.",
    )
    analyze_parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing *_5_end_freq and *_3_end_freq tables.",
    )
    analyze_parser.add_argument(
        "--output-pdf",
        type=Path,
        default=Path("mcmc_damage_summary_combined.pdf"),
        help="Path to the PDF summary report. Default: mcmc_damage_summary_combined.pdf",
    )
    analyze_parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("global-damage-results.csv"),
        help="Path to the CSV summary table. Default: global-damage-results.csv",
    )
    analyze_parser.add_argument(
        "--model-start",
        type=int,
        default=2,
        help="First position (inclusive) to include in the model. Default: 2.",
    )
    analyze_parser.add_argument(
        "--model-end",
        type=int,
        default=10,
        help="Last position (inclusive) to include in the model. Default: 10.",
    )
    analyze_parser.add_argument(
        "--draws",
        type=int,
        default=1000,
        help="Number of posterior samples to draw. Default: 1000.",
    )
    analyze_parser.add_argument(
        "--tune",
        type=int,
        default=1000,
        help="Number of tuning iterations. Default: 1000.",
    )
    analyze_parser.add_argument(
        "--chains",
        type=int,
        default=4,
        help="Number of MCMC chains. Default: 4.",
    )
    analyze_parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="Number of CPU cores to use. Default: 1.",
    )
    analyze_parser.add_argument(
        "--target-accept",
        type=float,
        default=0.95,
        help="Target acceptance rate for the sampler. Default: 0.95.",
    )
    analyze_parser.add_argument(
        "--progressbar",
        action="store_true",
        help="Display progress bars during sampling.",
    )
    analyze_parser.add_argument(
        "--hdi-prob",
        type=float,
        default=0.95,
        help="Probability mass for highest density intervals. Default: 0.95.",
    )
    analyze_parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible sampling.",
    )
    analyze_parser.add_argument(
        "--five-suffix",
        type=str,
        default="_5_end_freq",
        help="Suffix identifying 5′ frequency tables. Default: _5_end_freq.",
    )
    analyze_parser.add_argument(
        "--three-suffix",
        type=str,
        default="_3_end_freq",
        help="Suffix identifying 3′ frequency tables. Default: _3_end_freq.",
    )

    return parser


def _run_count(args: argparse.Namespace) -> None:
    for fastq in args.fastq:
        out5, out3 = write_frequency_tables(fastq, args.output_dir, n_bp=args.n_bp, trim=args.trim)
        print(f"Wrote {out5}")
        print(f"Wrote {out3}")


def _run_analyze(args: argparse.Namespace) -> None:
    results = analyze_frequency_tables(
        input_dir=args.input_dir,
        output_pdf=args.output_pdf,
        output_csv=args.output_csv,
        five_suffix=args.five_suffix,
        three_suffix=args.three_suffix,
        model_start=args.model_start,
        model_end=args.model_end,
        draws=args.draws,
        tune=args.tune,
        chains=args.chains,
        cores=args.cores,
        target_accept=args.target_accept,
        progressbar=args.progressbar,
        hdi_prob=args.hdi_prob,
        random_seed=args.seed,
    )
    print(f"PDF saved to: {args.output_pdf}")
    print(f"CSV saved to: {args.output_csv}")
    print(results)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "count":
        _run_count(args)
    elif args.command == "analyze":
        _run_analyze(args)
    else:  # pragma: no cover - defensive
        parser.error(f"Unknown command: {args.command}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
