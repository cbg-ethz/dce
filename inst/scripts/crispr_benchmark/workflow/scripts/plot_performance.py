from pathlib import Path

import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

import statannot


sns.set_context('talk')


def main(dname, out_dir):
    # prepare
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(dname / 'measures.csv').set_index(
        ['treatment', 'perturbed_gene', 'method', 'study', 'pathway']
    )

    # print statistics
    print(df.groupby('method').count())
    print(df.groupby(['method'])['roc_auc'].median())
    print(df.groupby(['method'])['roc_auc'].std())

    # create plots
    fig, ax = plt.subplots(figsize=(8, 6))

    sns.boxplot(
        data=df.reset_index(), x='method', y='roc_auc', order=['dce', 'cor', 'pcor']
    )
    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0.3))

    sns.stripplot(
        data=df.reset_index(), x='method', y='roc_auc', order=['dce', 'cor', 'pcor']
    )

    statannot.add_stat_annotation(
        ax,
        data=df.reset_index(),
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


if __name__ == '__main__':
    main(Path(snakemake.input.dname), Path(snakemake.output.out_dir))
