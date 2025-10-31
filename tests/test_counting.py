import pandas as pd

from damage_bayes.counting import count_base_frequencies


def test_counts_match_expected_values():
    df5, df3 = count_base_frequencies("data/subsample_1.fastq", n_bp=5, trim=1)

    expected_5 = pd.DataFrame(
        {
            "Position_from_5end": [2, 3, 4, 5],
            "A_freq": [2358, 2404, 2370, 2249],
            "T_freq": [2804, 2604, 2821, 2659],
            "C_freq": [2137, 2158, 2305, 2467],
            "G_freq": [2701, 2834, 2504, 2625],
            "Total": [10000, 10000, 10000, 10000],
        }
    )
    expected_3 = pd.DataFrame(
        {
            "Position_from_3end": [2, 3, 4, 5],
            "A_freq": [2939, 2730, 2843, 2832],
            "T_freq": [2365, 2418, 2368, 2240],
            "C_freq": [2658, 2628, 2478, 2515],
            "G_freq": [2038, 2224, 2311, 2413],
            "Total": [10000, 10000, 10000, 10000],
        }
    )

    pd.testing.assert_frame_equal(df5.reset_index(drop=True), expected_5)
    pd.testing.assert_frame_equal(df3.reset_index(drop=True), expected_3)
