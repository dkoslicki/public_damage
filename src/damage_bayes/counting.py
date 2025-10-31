"""Compute nucleotide frequency tables at the termini of sequencing reads."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Iterable, Tuple

import pandas as pd

from .fastq import read_fastq

__all__ = [
    "count_base_frequencies",
    "write_frequency_tables",
]


def _initialize_counter() -> dict[int, dict[str, int]]:
    return defaultdict(lambda: defaultdict(int))


def count_base_frequencies(
    fastq_path: str | Path,
    n_bp: int = 25,
    trim: int = 1,
    allowed_bases: Iterable[str] = ("A", "T", "C", "G"),
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Count nucleotides at both ends of reads in ``fastq_path``.

    Parameters
    ----------
    fastq_path:
        FASTQ file to process.
    n_bp:
        Number of terminal bases to inspect. Defaults to 25, matching the
        notebook implementation.
    trim:
        Number of positions to skip from the termini before counting. Defaults
        to 1.
    allowed_bases:
        Iterable of bases that should be tabulated. Any other characters are
        ignored to mirror the behaviour of the original AWK routine.

    Returns
    -------
    tuple of :class:`pandas.DataFrame`
        Two dataframes corresponding to the 5′ and 3′ ends respectively. Each
        dataframe contains columns ``Position_from_5end``/``Position_from_3end``,
        the nucleotide counts, and the total observations per position.
    """

    if n_bp <= 0:
        raise ValueError("n_bp must be a positive integer")
    if trim < 0:
        raise ValueError("trim must be non-negative")
    if trim >= n_bp:
        raise ValueError("trim must be less than n_bp")

    fastq_path = Path(fastq_path)
    counts5 = _initialize_counter()
    counts3 = _initialize_counter()
    totals5 = defaultdict(int)
    totals3 = defaultdict(int)

    allowed = tuple(allowed_bases)
    start_pos = trim + 1

    for record in read_fastq(str(fastq_path)):
        seq = record.sequence.upper()
        length = len(seq)
        if length == 0:
            continue

        for pos in range(start_pos, n_bp + 1):
            totals5[pos] += 1
            totals3[pos] += 1

            if length < pos:
                continue

            base5 = seq[pos - 1]
            base3 = seq[-pos]

            if base5 in allowed:
                counts5[pos][base5] += 1
            if base3 in allowed:
                counts3[pos][base3] += 1

    def build_dataframe(counts, totals, label):
        rows = []
        for pos in range(start_pos, n_bp + 1):
            row = {label: pos}
            for base in allowed:
                row[f"{base}_freq"] = counts[pos][base]
            row["Total"] = totals[pos]
            rows.append(row)
        return pd.DataFrame(rows)

    df5 = build_dataframe(counts5, totals5, "Position_from_5end")
    df3 = build_dataframe(counts3, totals3, "Position_from_3end")
    return df5, df3


def write_frequency_tables(
    fastq_path: str | Path,
    output_dir: str | Path,
    n_bp: int = 25,
    trim: int = 1,
) -> Tuple[Path, Path]:
    """Write 5′ and 3′ nucleotide frequency tables for ``fastq_path``.

    The output filenames follow the convention ``<name>_5_end_freq`` and
    ``<name>_3_end_freq`` where ``<name>`` is derived from the FASTQ filename
    without its extension.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df5, df3 = count_base_frequencies(fastq_path, n_bp=n_bp, trim=trim)

    base_name = Path(fastq_path).name
    if base_name.endswith(".fastq"):
        base_name = base_name[:-6]
    elif base_name.endswith(".fastq.gz"):
        base_name = base_name[:-9]

    out5 = output_dir / f"{base_name}_5_end_freq"
    out3 = output_dir / f"{base_name}_3_end_freq"

    df5.to_csv(out5, sep="\t", index=False)
    df3.to_csv(out3, sep="\t", index=False)

    return out5, out3
