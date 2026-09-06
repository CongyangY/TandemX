"""Six-panel diagnostic of observed joint-read interval calibration."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from benchmarks.challenge.schema import digest_file


def plot(run: Path, outdir: Path) -> None:
    validation = json.loads((run/'validation.json').read_text())
    if not validation['complete'] or validation['heldout_used']:
        raise ValueError('Require completed development calibration')
    rows = list(csv.DictReader((run/'metrics.tsv').open(), delimiter='\t'))
    if len(rows) != validation['family_conditions']:
        raise ValueError('Calibration row-count mismatch')
    outdir.mkdir(parents=True, exist_ok=False)
    inputs = {name: digest_file(run/name) for name in ('metrics.tsv', 'validation.json', 'inputs.tsv')}
    def available(r): return r['sampling_interval_low'] != 'NA'
    def covered(r): return r['truth_covered'] == 'True'
    def subset(**fields): return [r for r in rows if all(r[k] == str(v) for k, v in fields.items())]
    plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.size': 8,
                         'axes.spines.top': False, 'axes.spines.right': False})
    colors = ['#0072B2', '#E69F00', '#009E73']
    coverages = ['1', '5', '20']
    seeds = sorted({r['seed'] for r in rows})
    errors = ['error_free', 'iid_low_indel', 'iid_high_indel']
    fig, grid = plt.subplots(2, 3, figsize=(12, 7.3))
    axes = grid.ravel()
    for i, ax in enumerate(axes):
        ax.text(-.16, 1.08, 'ABCDEF'[i], transform=ax.transAxes, fontweight='bold', fontsize=13)
    statuses = ['approximate_95pct_joint_read_log_normal', 'insufficient_effective_reads', 'no_finite_extrapolation']
    bottoms = np.zeros(3)
    for status, color, label in zip(statuses, [colors[0], '#bdbdbd', '#D55E00'], ['Interval available', 'Sparse read support', 'No finite fit']):
        counts = [sum(r['interval_status'] == status for r in subset(coverage=c)) for c in coverages]
        axes[0].bar(coverages, counts, bottom=bottoms, color=color, label=label)
        bottoms += counts
    if set(r['interval_status'] for r in rows)-set(statuses):
        raise ValueError('Unrepresented interval status; update plot explicitly')
    axes[0].set(title='Availability across read depths', xlabel='Nominal simulated coverage (×)', ylabel='Family conditions', ylim=(0, max(bottoms)*1.4))
    axes[0].legend(fontsize=7, frameon=False, loc='upper left')
    for seed, color in zip(seeds, colors):
        rates = []
        for c in coverages:
            a = [r for r in subset(seed=seed, coverage=c) if available(r)]
            rates.append(sum(covered(r) for r in a)/len(a) if a else np.nan)
        axes[1].plot(coverages, rates, 'o-', color=color, label='Genome '+seed)
    axes[1].axhline(.95, ls='--', color='black', lw=.8)
    axes[1].set(title='Coverage conditional on availability', xlabel='Nominal simulated coverage (×)', ylabel='Fraction containing truth', ylim=(.75,1.025))
    axes[1].legend(frameon=False, fontsize=7)
    labels = []
    for pos, (div, error) in enumerate((d,e) for d in ['0.0','0.02'] for e in errors):
        labels.append(('0%' if div == '0.0' else '2%')+'\n'+{'error_free':'none','iid_low_indel':'low','iid_high_indel':'high'}[error])
        for j, (seed, color) in enumerate(zip(seeds, colors)):
            a = [r for r in subset(seed=seed, coverage=20, array_scope='factorial', unit_substitution_rate=div, error_model=error) if available(r)]
            axes[2].scatter(pos+(j-1)*.15, sum(covered(r) for r in a)/len(a) if a else np.nan, color=color, s=22)
    axes[2].axhline(.95, ls='--', color='black', lw=.8)
    axes[2].set(xticks=range(6), xticklabels=labels, title='20×: divergence and read errors', ylabel='Fraction containing truth', xlabel='Unit divergence / read-error tier', ylim=(.7,1.025))
    for error, color, label in zip(errors, colors, ['No read errors', 'Low errors', 'High errors']):
        g = [r for r in subset(coverage=20, error_model=error) if r['estimated_copy_number'] != 'NA']
        axes[3].scatter([float(r['truth_copies']) for r in g], [float(r['estimated_copy_number']) for r in g], s=12, alpha=.35, color=color, label=label)
    axes[3].plot([10,10000], [10,10000], '--', color='black', lw=.8)
    axes[3].set(xscale='log', yscale='log', xlabel='Planted copy number', ylabel='Estimated copy number', title='20× point estimates (all finite fits)')
    axes[3].legend(fontsize=7, frameon=False)
    for c, color in zip(coverages, colors):
        a = [r for r in subset(coverage=c) if available(r)]
        axes[4].scatter([float(r['minimum_effective_reads']) for r in a], [float(r['relative_interval_width']) for r in a], s=12, alpha=.5, color=color, label=c+'×')
    axes[4].set(xscale='log', yscale='log', xlabel='Minimum effective reads across k', ylabel='Interval width / planted copies', title='Sampling support and interval width')
    axes[4].legend(fontsize=7, frameon=False)
    copies = sorted({float(r['truth_copies']) for r in rows})
    for j, (c, color) in enumerate(zip(coverages, colors)):
        fractions = []
        for n in copies:
            g = [r for r in subset(coverage=c) if float(r['truth_copies']) == n]
            fractions.append(sum(not available(r) for r in g)/len(g))
        axes[5].bar(np.arange(len(copies))+(j-1)*.25, fractions, width=.23, color=color, label=c+'×')
    axes[5].set(xticks=range(len(copies)), xticklabels=[str(int(n)) for n in copies], xlabel='Planted copy number', ylabel='Fraction without interval', title='Missingness retained by abundance', ylim=(0,1.05))
    axes[5].legend(fontsize=7, frameon=False)
    fig.suptitle('Joint-read multi-k sampling uncertainty: three development genomes, 27 conditions', fontsize=12)
    fig.text(.5,.012,'Known founder catalogue; k = 15, 21, 27, 31; ≥20 effective reads. Families and nested conditions are dependent. Held-out data unused.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0,.035,1,.95), h_pad=2, w_pad=2)
    for extension in ('pdf','svg','png'):
        fig.savefig(outdir/f'joint_uncertainty.{extension}', dpi=180)
    plt.close(fig)
    shutil.copyfile(run/'metrics.tsv', outdir/'source_data.tsv')
    shutil.copyfile(__file__, outdir/Path(__file__).name)
    if any(digest_file(run/name) != value for name, value in inputs.items()):
        raise ValueError('Inputs changed while plotting')
    (outdir/'figure_receipt.json').write_text(json.dumps(dict(complete=True, inputs=inputs,
        source=str(run), panels=6, script_sha256=digest_file(Path(__file__)),
        files={p.name:digest_file(p) for p in outdir.iterdir() if p.is_file()}), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    plot(args.run, args.outdir)
