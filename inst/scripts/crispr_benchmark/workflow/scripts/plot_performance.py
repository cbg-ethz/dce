from pathlib import Path

import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

import statannot
from natsort import natsorted


sns.set_context('talk')


def main(dname, out_dir):
    # prepare
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(dname / 'measures.csv')

    # print statistics
    print(df.groupby('method').count())
    print(df.groupby(['method'])['roc_auc'].median())
    print(df.groupby(['method'])['roc_auc'].std())

    # aggregated plot
    fig, ax = plt.subplots(figsize=(8, 6))

    sns.boxplot(data=df, x='method', y='roc_auc', order=['dce', 'cor', 'pcor'])
    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0.3))

    sns.stripplot(data=df, x='method', y='roc_auc', order=['dce', 'cor', 'pcor'])

    statannot.add_stat_annotation(
        ax,
        data=df,
        x='method',
        y='roc_auc',
        order=['dce', 'cor', 'pcor'],
        box_pairs=[('dce', 'cor'), ('dce', 'pcor')],
        test='Wilcoxon',
        text_format='simple',
        loc='outside',
        verbose=2,
    )

    ax.set_xlabel('Method')
    ax.set_ylabel('ROC-AUC')

    fig.tight_layout()
    fig.savefig(out_dir / 'method_comparison.pdf')

    # stratified plot
    g = sns.catplot(
        data=df,
        x='method',
        y='roc_auc',
        hue='perturbed_gene',
        row='treatment',
        kind='box',
        hue_order=natsorted(df['perturbed_gene'].unique()),
        aspect=2,
    )

    g.set_axis_labels('Method', 'ROC-AUC')
    g._legend.set_title('Perturbed gene(s)')

    g.savefig(out_dir / 'method_comparison_stratified.pdf')


if __name__ == '__main__':
    main(Path(snakemake.input.dname), Path(snakemake.output.out_dir))
