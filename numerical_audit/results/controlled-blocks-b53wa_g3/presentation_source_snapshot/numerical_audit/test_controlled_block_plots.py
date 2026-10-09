"""Regression for a below-floor first case and zero errors in log plots."""
import matplotlib.pyplot as plt
import numpy as np
import pytest

from qhat.numerical_audit.plot_controlled_block_orders import STYLES, plot_summary


def synthetic_summary():
    comparisons = []
    for total_time in (.25, 1.):
        comparisons.extend([
            {'case': 'first_below_floor', 'role': 'primary', 'T': total_time,
             'random_median_E': 0., 'random_q25_E': 0., 'random_q75_E': 5e-14,
             'deterministic_E': {name: 0. for name, *_ in STYLES}},
            {'case': 'second_visible', 'role': 'primary', 'T': total_time,
             'random_median_E': 1e-6, 'random_q25_E': 5e-7, 'random_q75_E': 2e-6,
             'deterministic_E': {name: (0. if name == 'block_signed' else 1e-6) for name, *_ in STYLES}},
        ])
    return {'analysis_floor': 1e-12, 'comparisons': comparisons}


@pytest.mark.parametrize('normalized', [False, True])
def test_first_omitted_case_keeps_legend_and_zero_errors_have_explicit_floor_markers(tmp_path, normalized):
    summary = synthetic_summary()
    fig, axes, checks = plot_summary(summary, normalized=normalized)
    labels = [text.get_text() for text in axes[0].get_legend().get_texts()]
    assert all(label in labels for _, label, *_ in STYLES)
    assert any('display floor' in label for label in labels)
    assert len(checks['omitted_ratio_cases']) == (2 if normalized else 0)
    assert len(checks['capped_deterministic_points']) == (2 if normalized else 10)
    if normalized:
        assert all(any(text.get_text() == 'ratio omitted' for text in ax.texts) for ax in axes)
    # Force a complete render and ensure no deterministic point gets an invalid
    # zero coordinate on the logarithmic axis.
    fig.savefig(tmp_path/f'edge_{normalized}.png')
    for ax in axes:
        for collection in ax.collections:
            offsets = collection.get_offsets()
            if len(offsets) and collection.__class__.__name__ == 'PathCollection':
                assert np.all(offsets[:, 1] > 0)
    assert summary == synthetic_summary()  # presentation never mutates raw data
    plt.close(fig)
