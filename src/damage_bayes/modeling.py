"""Bayesian modelling utilities for ancient DNA damage patterns."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc as pm
from matplotlib.backends.backend_pdf import PdfPages

__all__ = [
    "DamageModelResult",
    "run_damage_model",
    "analyze_frequency_tables",
]


@dataclass
class DamageModelResult:
    """Summary of posterior damage estimates for a single strand direction."""

    positions: np.ndarray
    mean_damage: np.ndarray
    hdi_lower: np.ndarray
    hdi_upper: np.ndarray
    alpha_mean: float
    beta_mean: float


_DEF_POSITIONS = np.arange(2, 26)


def run_damage_model(
    positions: Sequence[int] | np.ndarray,
    success_counts: Sequence[int] | np.ndarray,
    fail_counts: Sequence[int] | np.ndarray,
    draws: int = 1000,
    tune: int = 1000,
    chains: int = 4,
    cores: int = 1,
    target_accept: float = 0.95,
    progressbar: bool = False,
    hdi_prob: float = 0.95,
    random_seed: int | None = None,
) -> DamageModelResult:
    """Run the Bayesian damage model described in the notebook."""

    positions = np.asarray(positions)
    success_counts = np.asarray(success_counts)
    fail_counts = np.asarray(fail_counts)

    total_counts = success_counts + fail_counts

    with pm.Model() as model:
        alpha = pm.Beta("alpha", alpha=1, beta=1)
        beta = pm.HalfNormal("beta", sigma=1)
        p = alpha * pm.math.exp(-beta * (positions - 1))
        pm.Binomial("obs", n=total_counts, p=p, observed=success_counts)

        trace = pm.sample(
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            chains=chains,
            cores=cores,
            progressbar=progressbar,
            return_inferencedata=True,
            random_seed=random_seed,
        )

    summary = az.summary(trace, var_names=["alpha", "beta"])
    alpha_mean = float(summary.loc["alpha", "mean"])
    beta_mean = float(summary.loc["beta", "mean"])

    posterior_samples = trace.posterior.stack(samples=("chain", "draw"))
    alpha_samples = posterior_samples["alpha"].values
    beta_samples = posterior_samples["beta"].values

    pos_range = _DEF_POSITIONS
    mean_p = []
    hdi_lower = []
    hdi_upper = []

    for pos in pos_range:
        probs = alpha_samples * np.exp(-beta_samples * (pos - 1))
        mean_p.append(np.mean(probs))
        hdi = az.hdi(probs, hdi_prob=hdi_prob)
        hdi_lower.append(hdi[0])
        hdi_upper.append(hdi[1])

    return DamageModelResult(
        positions=pos_range,
        mean_damage=np.array(mean_p),
        hdi_lower=np.array(hdi_lower),
        hdi_upper=np.array(hdi_upper),
        alpha_mean=alpha_mean,
        beta_mean=beta_mean,
    )


def _subset_positions(df: pd.DataFrame, column: str, start: int, end: int) -> pd.DataFrame:
    return df[df[column].between(start, end)].copy()


@dataclass
class PairedSample:
    name: str
    five_prime: Path
    three_prime: Path


def _pair_frequency_tables(
    input_dir: Path,
    five_suffix: str,
    three_suffix: str,
) -> list[PairedSample]:
    files_5 = sorted(input_dir.glob(f"*{five_suffix}"))
    files_3 = sorted(input_dir.glob(f"*{three_suffix}"))

    pairing = {}
    for path in files_5:
        pairing[path.name.replace(five_suffix, "")] = {"five": path}
    for path in files_3:
        pairing.setdefault(path.name.replace(three_suffix, ""), {}).update({"three": path})

    paired = []
    for sample, sides in pairing.items():
        if "five" in sides and "three" in sides:
            paired.append(PairedSample(sample, sides["five"], sides["three"]))
    return paired


def analyze_frequency_tables(
    input_dir: str | Path,
    output_pdf: str | Path,
    output_csv: str | Path,
    five_suffix: str = "_5_end_freq",
    three_suffix: str = "_3_end_freq",
    model_start: int = 2,
    model_end: int = 10,
    draws: int = 1000,
    tune: int = 1000,
    chains: int = 4,
    cores: int = 1,
    target_accept: float = 0.95,
    progressbar: bool = False,
    hdi_prob: float = 0.95,
    random_seed: int | None = None,
) -> pd.DataFrame:
    """Analyse paired 5′/3′ frequency tables and persist summary outputs."""

    input_dir = Path(input_dir)
    output_pdf = Path(output_pdf)
    output_csv = Path(output_csv)

    paired_samples = _pair_frequency_tables(input_dir, five_suffix, three_suffix)

    results = []

    with PdfPages(output_pdf) as pdf:
        for sample in paired_samples:
            alpha_ct = beta_ct = alpha_ga = beta_ga = None
            try:
                df_5 = pd.read_csv(sample.five_prime, sep="\t")
                df_3 = pd.read_csv(sample.three_prime, sep="\t")

                df5sub = _subset_positions(df_5, "Position_from_5end", model_start, model_end)
                df3sub = _subset_positions(df_3, "Position_from_3end", model_start, model_end)

                result_ct = run_damage_model(
                    df5sub["Position_from_5end"].values,
                    df5sub["T_freq"].values,
                    df5sub["C_freq"].values,
                    draws=draws,
                    tune=tune,
                    chains=chains,
                    cores=cores,
                    target_accept=target_accept,
                    progressbar=progressbar,
                    hdi_prob=hdi_prob,
                    random_seed=random_seed,
                )

                result_ga = run_damage_model(
                    df3sub["Position_from_3end"].values,
                    df3sub["A_freq"].values,
                    df3sub["G_freq"].values,
                    draws=draws,
                    tune=tune,
                    chains=chains,
                    cores=cores,
                    target_accept=target_accept,
                    progressbar=progressbar,
                    hdi_prob=hdi_prob,
                    random_seed=random_seed,
                )

                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
                clean_name = sample.name.replace("filtered_25bp_100k_", "")

                ax1.plot(result_ct.positions, result_ct.mean_damage, color="firebrick", label="C→T (5′)")
                ax1.fill_between(result_ct.positions, result_ct.hdi_lower, result_ct.hdi_upper, color="firebrick", alpha=0.3)
                ax1.set_title(f"5′ Damage: {clean_name}")
                ax1.set_xlabel("Position from 5′")
                ax1.set_ylabel("Damage Rate")
                ax1.grid(True)

                ax2.plot(result_ga.positions, result_ga.mean_damage, color="darkblue", label="G→A (3′)")
                ax2.fill_between(result_ga.positions, result_ga.hdi_lower, result_ga.hdi_upper, color="darkblue", alpha=0.3)
                ax2.set_title(f"3′ Damage: {clean_name}")
                ax2.set_xlabel("Position from 3′")
                ax2.set_ylabel("Damage Rate")
                ax2.set_xlim(ax2.get_xlim()[::-1])
                ax2.grid(True)

                plt.suptitle(clean_name, fontsize=14, fontweight="bold")
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

                alpha_ct = result_ct.alpha_mean
                beta_ct = result_ct.beta_mean
                alpha_ga = result_ga.alpha_mean
                beta_ga = result_ga.beta_mean

            except Exception as exc:  # pragma: no cover - defensive logging
                print(f"Error processing {sample.name}: {exc}")
            finally:
                results.append(
                    {
                        "Sample": sample.name.replace("filtered_25bp_100k_", ""),
                        "Alpha_CtoT_Mean": round(alpha_ct, 3) if alpha_ct is not None else None,
                        "Beta_CtoT_Mean": round(beta_ct, 3) if beta_ct is not None else None,
                        "Alpha_GtoA_Mean": round(alpha_ga, 3) if alpha_ga is not None else None,
                        "Beta_GtoA_Mean": round(beta_ga, 3) if beta_ga is not None else None,
                    }
                )

    df_results = pd.DataFrame(results)
    df_results.to_csv(output_csv, index=False)
    return df_results
