import numpy as np

from damage_bayes.modeling import DamageModelResult, run_damage_model


def test_run_damage_model_returns_expected_shapes():
    positions = np.array([2, 3, 4])
    successes = np.array([30, 20, 10])
    failures = np.array([70, 80, 90])

    result = run_damage_model(
        positions,
        successes,
        failures,
        draws=50,
        tune=50,
        chains=1,
        cores=1,
        progressbar=False,
        random_seed=123,
    )

    assert isinstance(result, DamageModelResult)
    assert result.positions.shape == (24,)  # positions 2-25 inclusive
    assert result.mean_damage.shape == result.positions.shape
    assert np.all(result.hdi_upper >= result.hdi_lower)
    assert result.alpha_mean > 0
    assert result.beta_mean >= 0
