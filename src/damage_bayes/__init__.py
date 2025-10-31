"""Damage Bayes: Bayesian modeling of ancient DNA damage patterns."""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("damage-bayes")
except PackageNotFoundError:  # pragma: no cover - during local development
    __version__ = "0.0.0"

__all__ = ["__version__"]
