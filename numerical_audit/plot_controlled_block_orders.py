#!/usr/bin/env python3
"""Export descriptive plots of a completed fixed-block permutation campaign.

Values at/below the declared interpretation floor remain unchanged in the data.
Only their log-plot positions are capped and marked by downward triangles;
ratios with a random median at/below the floor are explicitly omitted.
"""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/qhat-audit-mpl')
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

STYLES = (
    ('block_signed', 'Signed parent', '#0072B2', 'o', -.20),
    ('block_parent_magnitude', 'Parent magnitude', '#D55E00', '^', -.065),
    ('block_l1', 'Block coefficient l1', '#009E73', 's', .065),
    ('block_reverse_signed', 'Reverse signed', '#CC79A7', 'x', .20),
)


def plot_summary(summary, *, normalized=False):
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.8), sharex=True, constrained_layout=True)
    floor = summary['analysis_floor']
    metadata = {'normalized': normalized, 'omitted_ratio_cases': [], 'capped_deterministic_points': []}
    for ax, total_time in zip(axes, (.25, 1.)):
        records = [c for c in summary['comparisons'] if c['role'] == 'primary' and c['T'] == total_time]
        for x, case in enumerate(records):
            denominator = case['random_median_E'] if normalized else 1.
            if normalized and denominator <= floor:
                metadata['omitted_ratio_cases'].append({'case': case['case'], 'T': total_time})
                ax.plot(x, .045, marker='x', color='.45', transform=ax.get_xaxis_transform(), clip_on=False)
                ax.text(x, .075, 'ratio omitted', ha='center', va='bottom', rotation=90,
                        fontsize=9.5, color='.4', transform=ax.get_xaxis_transform())
                continue
            q25 = max(case['random_q25_E'], floor)/denominator
            q75 = max(case['random_q75_E'], floor)/denominator
            ax.vlines(x, q25, q75, color='.75', linewidth=11, zorder=1)
            ax.plot(x, max(case['random_median_E'], floor)/denominator, marker='_', color='.35', ms=16, zorder=2)
            if case['random_q25_E'] <= floor:
                ax.plot(x, floor/denominator, marker='v', color='.55', ms=6, zorder=2)
            for name, label, color, marker, offset in STYLES:
                raw_error = case['deterministic_E'][name]
                if raw_error <= floor:
                    metadata['capped_deterministic_points'].append({'case': case['case'], 'T': total_time, 'ordering': name, 'raw_E': raw_error})
                    ax.scatter(x+offset, floor/denominator, edgecolor=color, facecolor='none', marker='v', s=45, zorder=3)
                else:
                    ax.scatter(x+offset, raw_error/denominator, color=color, marker=marker, s=42, zorder=3)
        ax.set_yscale('log')
        ax.grid(axis='y', alpha=.18)
        ax.set_title(f'T = {total_time:g}, r = 100; fixed ownership, Hartree–Fock state', fontsize=11)
        ax.set_ylabel('E / median random E' if normalized else 'Full-state overlap error E')
        ax.axhline(1. if normalized else floor, color='.4', linestyle=':', linewidth=1.)
        ax.set_xticks(range(len(records)), [c['case'].replace('_', '\n') for c in records], fontsize=10)
    handles = [Line2D([], [], color=color, marker=marker, linestyle='none', label=label, markersize=6)
               for _, label, color, marker, _ in STYLES]
    handles += [Line2D([], [], color='.65', linewidth=7, label='Random IQR (20 seeds)'),
                Line2D([], [], color='.3', marker='v', markerfacecolor='none', linestyle='none', label=f'E ≤ {floor:g} (display floor)')]
    axes[0].legend(handles=handles, loc='best', ncols=2, fontsize=10)
    fig.suptitle('Whole-block sorting with unchanged block membership\nGray bars: random IQR; ticks: random median (20 seeds)', fontsize=11)
    return fig, axes, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    args = parser.parse_args()
    completion = json.loads((args.run/'completion.json').read_text())
    summary = json.loads((args.run/'summary.json').read_text())
    if completion.get('status') != 'complete' or completion['rows'] != summary['rows']:
        raise RuntimeError('Only a completed campaign can be plotted.')
    checks = []
    for normalized in (False, True):
        fig, _, metadata = plot_summary(summary, normalized=normalized)
        checks.append(metadata)
        stem = 'block_errors_relative_to_random' if normalized else 'block_errors_absolute'
        fig.savefig(args.run/f'{stem}.pdf', bbox_inches='tight')
        fig.savefig(args.run/f'{stem}.png', dpi=180, bbox_inches='tight')
        plt.close(fig)
    with (args.run/'presentation_plot_checks.json').open('x') as stream:
        json.dump(checks, stream, indent=2)
        stream.write('\n')
    print('Plots saved in', args.run)


if __name__ == '__main__':
    main()
