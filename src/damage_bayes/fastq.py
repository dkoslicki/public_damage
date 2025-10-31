"""Utilities for streaming FASTQ records."""
from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass


@dataclass(frozen=True)
class FastqRecord:
    """A single FASTQ record."""

    header: str
    sequence: str
    separator: str
    quality: str


def _read_nonempty_line(handle) -> str:
    line = handle.readline()
    if not line:
        return ""
    return line.rstrip("\n\r")


def read_fastq(path: str) -> Iterator[FastqRecord]:
    """Yield :class:`FastqRecord` objects from ``path``.

    The parser is resilient to sequences and quality strings that span multiple
    lines. It yields each record as soon as the quality string for a header has
    been fully collected.
    """

    with open(path, "r", encoding="utf-8") as handle:
        while True:
            header = _read_nonempty_line(handle)
            if not header:
                return
            if not header.startswith("@"):
                raise ValueError(f"Invalid FASTQ header: {header!r}")

            # Collect sequence lines until the separator is reached.
            seq_lines: list[str] = []
            while True:
                line = handle.readline()
                if not line:
                    raise ValueError("Unexpected end of file before '+' separator")
                if line.startswith("+"):
                    separator = line.rstrip("\n\r")
                    break
                seq_lines.append(line.strip())
            sequence = "".join(seq_lines)

            # Collect quality lines until the full length of the sequence is covered.
            quality_chunks: list[str] = []
            remaining = len(sequence)
            while remaining > 0:
                chunk = _read_nonempty_line(handle)
                if not chunk:
                    raise ValueError("Unexpected end of file in quality scores")
                quality_chunks.append(chunk)
                remaining -= len(chunk)

            quality = "".join(quality_chunks)
            if len(quality) != len(sequence):
                raise ValueError(
                    "Quality string length does not match sequence length for "
                    f"record {header!r}"
                )

            yield FastqRecord(header=header, sequence=sequence, separator=separator, quality=quality)
