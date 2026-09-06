"""Four-panel complete-library and nested-sample distribution audit."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from benchmarks.challenge.schema import digest_file, read_table, write_table


def histogram(path: Path, value_field: str, expected_reads: int) -> Counter:
    result = Counter()
    for row in read_table(path):
        value, count = int(row[value_field]), int(row['read_count'])
        if value < 0 or count < 0:
            raise ValueError('Negative QC histogram value/count')
        result[value] += count
    if sum(result.values()) != expected_reads:
        raise ValueError('Histogram total differs from receipt: '+str(path))
    return result


def plot(qcdir: Path, samplesdir: Path, outdir: Path, label: str) -> None:
    qc = json.loads((qcdir/'qc.json').read_text())
    sampling = json.loads((samplesdir/'sampling_receipt.json').read_text())
    if not qc.get('complete') or not sampling.get('complete') or qc['input_sha256'] != sampling['source_sha256']:
        raise ValueError('Require complete matching QC and sample receipts')
    inputs = [qcdir/'qc.json', samplesdir/'sampling_receipt.json', samplesdir/'sampling_plan.json']
    datasets = [dict(name=f'Full library (n={qc["read_count"]:,})', prefix=None, read_count=qc['read_count'],
                     total_bases=qc['total_bases'], expected_total_bases=qc['total_bases'], color='#222222', probability=1.0)]
    colors = ['#4477AA', '#CC6677', '#228833', '#AA3377', '#66CCEE', '#EE7733']
    for i, row in enumerate(sampling['samples']):
        probability = float(row['inclusion_probability'])
        datasets.append(dict(name=f'{probability:.3g} fraction (n={row["read_count"]:,})', prefix=row['sample_id'],
                             read_count=row['read_count'], total_bases=row['total_bases'],
                             expected_total_bases=row['expected_total_bases'], color=colors[i % len(colors)], probability=probability))
    plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.family': 'DejaVu Sans',
                         'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.4), layout='constrained')
    source = []
    for dataset in datasets:
        prefix, count = dataset['prefix'], dataset['read_count']
        lengths = qcdir/'length_histogram.tsv' if prefix is None else samplesdir/(prefix+'.length_histogram.tsv')
        joint = qcdir/'joint_distribution.tsv' if prefix is None else samplesdir/(prefix+'.joint_distribution.tsv')
        inputs.extend([lengths, joint])
        for panel, path, field, axis, scale, cumulative in (
                ('A', lengths, 'length_bp', axes[0, 0], .001, True),
                ('B', joint, 'gc_bin_percent', axes[0, 1], 1, False),
                ('C', joint, 'mean_quality_bin_phred', axes[1, 0], 1, True)):
            hist = histogram(path, field, count)
            total = 0
            x, y = [], []
            for value, number in sorted(hist.items()):
                total += number
                fraction = (total if cumulative else number)/count if count else None
                source.append(dict(panel=panel, dataset=dataset['name'], metric=field,
                                   x=value*scale, value=fraction, read_count=number))
                x.append(value*scale); y.append(fraction)
            if count:
                axis.step(x, y, where='post' if cumulative else 'mid', color=dataset['color'],
                          label=dataset['name'], lw=1.4 if prefix else 2.2, alpha=.85)
    for axis, title, xlabel, ylabel in (
            (axes[0, 0], 'A   Read length', 'Read length (kb)', 'Cumulative read fraction'),
            (axes[0, 1], 'B   Read GC content', 'GC bin lower bound (%)', 'Read fraction per 1% bin'),
            (axes[1, 0], 'C   Reported read quality', 'Read Phred bin lower bound (5-unit bins)', 'Cumulative read fraction')):
        axis.set_title(title, loc='left', fontweight='bold')
        axis.set(xlabel=xlabel, ylabel=ylabel)
        axis.grid(axis='y', alpha=.15)
    axes[0, 0].set_xlim(left=0)
    axes[0, 0].set_ylim(0, 1.02)
    axes[1, 0].set_ylim(0, 1.02)
    axes[0, 0].legend(fontsize=7.8, loc='lower right', frameon=False)
    axis = axes[1, 1]
    labels = []
    for i, dataset in enumerate(datasets):
        observed = dataset['total_bases']/1e9
        expected = dataset['expected_total_bases']/1e9
        labels.append('Full' if dataset['prefix'] is None else f'{dataset["probability"]:.3g} fraction')
        axis.scatter(observed, i, color=dataset['color'], s=65, alpha=.8)
        axis.scatter(expected, i, color='black', marker='|', s=90, zorder=3)
        if observed:
            axis.annotate(f'{observed:.4g} Gb', (observed, i), xytext=(5, 0), textcoords='offset points', va='center', fontsize=8)
        source.extend([dict(panel='D', dataset=dataset['name'], metric=metric, x=i, value=value,
                            read_count=dataset['read_count']) for metric, value in [('observed_Gb', observed), ('expected_Gb', expected)]])
    axis.set_yticks(range(len(labels)), labels)
    axis.invert_yaxis()
    axis.set_xscale('log')
    maximum = max(d['total_bases'] for d in datasets)/1e9
    positive = [d['total_bases']/1e9 for d in datasets if d['total_bases']]
    axis.set_xlim(min(positive)/2, maximum*3)
    axis.set_xlabel('Selected sequence bases (Gb; log scale)')
    axis.set_title('D   Observed data volume', loc='left', fontweight='bold')
    fig.suptitle(label+' — complete input and nested sampling QC', fontsize=12, fontweight='bold')
    fig.supxlabel('D: coloured dots show observed bases; black ticks show expected bases. Fractions sample whole reads, not fixed base totals.\n'
                  'Read quality derives from reported error probabilities; neither quality nor nominal coverage establishes biological accuracy.', fontsize=8)
    outdir.mkdir(parents=True, exist_ok=False)
    write_table(outdir/'panel_source.tsv', source, ['panel', 'dataset', 'metric', 'x', 'value', 'read_count'])
    outputs = []
    for extension in ('pdf', 'svg', 'png'):
        path = outdir/('input_qc.'+extension)
        fig.savefig(path, dpi=180)
        outputs.append(path)
    plt.close(fig)
    provenance = dict(inputs={str(p.resolve()): digest_file(p) for p in inputs},
                      source_sha256=digest_file(outdir/'panel_source.tsv'), script_sha256=digest_file(Path(__file__)),
                      outputs={p.name: digest_file(p) for p in outputs},
                      interpretation='QC of one complete included library; samples are not biological replicates')
    (outdir/'figure_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qcdir', type=Path, required=True)
    parser.add_argument('--samplesdir', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--label', required=True)
    args = parser.parse_args()
    plot(args.qcdir, args.samplesdir, args.outdir, args.label)
