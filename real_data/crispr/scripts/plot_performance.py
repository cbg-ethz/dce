from pathlib import Path

import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt


def main(dname, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(dname / 'measures.csv').set_index(['study', 'treatment', 'perturbed_gene', 'pathway'])

    # general overview
    plt.figure(figsize=(8, 6))

    sns.boxplot(
        data=pd.melt(df.drop(columns=['optimal_roc_threshold'])),
        x='variable', y='value')

    plt.xlabel('Performance measure')
    plt.ylabel('Score')

    plt.tight_layout()
    plt.savefig(out_dir / 'overview_boxplot.pdf')


if __name__ == '__main__':
    main(Path(snakemake.input.dname), Path(snakemake.output.out_dir))
