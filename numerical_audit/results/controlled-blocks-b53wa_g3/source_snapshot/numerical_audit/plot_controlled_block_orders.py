#!/usr/bin/env python3
"""Export descriptive plots of a completed fixed-block permutation campaign."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/qhat-audit-mpl')
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    args = parser.parse_args()
    completion = json.loads((args.run/'completion.json').read_text())
    summary = json.loads((args.run/'summary.json').read_text())
    if completion.get('status') != 'complete' or completion['rows'] != summary['rows']:
        raise RuntimeError('Only a completed campaign can be plotted.')
    styles = [
        ('block_signed', 'Signed parent', '#0072B2', 'o', -.20),
        ('block_parent_magnitude', 'Parent magnitude', '#D55E00', '^', -.065),
        ('block_l1', 'Block coefficient l1', '#009E73', 's', .065),
        ('block_reverse_signed', 'Reverse signed', '#CC79A7', 'x', .20),
    ]
    for normalized in (False, True):
        fig, axes = plt.subplots(2, 1, figsize=(13, 8.2), sharex=True, constrained_layout=True)
        for ax, total_time in zip(axes, (.25, 1.)):
            records = [c for c in summary['comparisons'] if c['role'] == 'primary' and c['T'] == total_time]
            for x, c in enumerate(records):
                denominator = c['random_median_E'] if normalized else 1.
                if normalized and denominator <= summary['analysis_floor']:
                    continue
                q25, q75 = c['random_q25_E']/denominator, c['random_q75_E']/denominator
                ax.vlines(x, q25, q75, color='.75', linewidth=11, zorder=1)
                ax.plot(x, c['random_median_E']/denominator, marker='_', color='.35', ms=16, zorder=2)
                for name, label, color, marker, offset in styles:
                    ax.scatter(x+offset, c['deterministic_E'][name]/denominator,
                               color=color, marker=marker, s=42, label=label if x == 0 else None, zorder=3)
            ax.set_yscale('log')
            ax.grid(axis='y', alpha=.18)
            ax.set_title(f'T = {total_time:g}, r = 100; fixed ownership, Hartree–Fock state', fontsize=11)
            ax.set_ylabel('E / median random E' if normalized else 'Full-state overlap error E')
            if normalized:
                ax.axhline(1., color='.4', linestyle=':', linewidth=1.)
            else:
                ax.axhline(summary['analysis_floor'], color='.4', linestyle=':', linewidth=.8)
            ax.set_xticks(range(len(records)), [c['case'].replace('_', '\n') for c in records], fontsize=9)
        axes[0].legend(loc='best', ncols=2, fontsize=9)
        fig.suptitle('Whole-block sorting with unchanged block membership\nGray bars: random-order interquartile ranges; gray ticks: random medians (20 seeds)', fontsize=12)
        stem = 'block_errors_relative_to_random' if normalized else 'block_errors_absolute'
        fig.savefig(args.run/f'{stem}.pdf', bbox_inches='tight')
        fig.savefig(args.run/f'{stem}.png', dpi=180, bbox_inches='tight')
        plt.close(fig)
    print('Plots saved in', args.run)


if __name__ == '__main__':
    main()
