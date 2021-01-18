from pathlib import Path

import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt


def main(dname, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(dname / 'measures.csv').set_index(['method', 'study', 'treatment', 'perturbed_gene', 'pathway'])

    # general overview of DCE
    plt.figure(figsize=(8, 6))

    sns.boxplot(
        data=pd.melt(df.loc['dce'].drop(columns=['optimal_roc_threshold'])),
        x='variable', y='value')

    plt.xlabel('Performance measure')
    plt.ylabel('Score')

    plt.tight_layout()
    plt.savefig(out_dir / 'overview_boxplot.pdf')

    # compare methods
    df_sub = pd.melt(df.reset_index(level=['method', 'treatment'])[['method', 'treatment', 'roc_auc']], id_vars=['method', 'treatment'])

    plt.figure(figsize=(8, 6))

    sns.boxplot(
        data=df_sub,
        x='method', y='value', hue='treatment')

    plt.xlabel('Method')
    plt.ylabel('ROC-AUC')

    plt.tight_layout()
    plt.savefig(out_dir / 'method_comparison.pdf')


if __name__ == '__main__':
    main(Path(snakemake.input.dname), Path(snakemake.output.out_dir))
