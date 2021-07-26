from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx

import seaborn as sns
import matplotlib.pyplot as plt

from tqdm import tqdm


sns.set_context('talk')


def read_data(dir_list):
    tmp = []
    for dir_ in tqdm(dir_list):
        # parse input
        *_, pathway, gene, study, treatment, deconf_param = dir_.split('/')
        gene_list = gene.split(',')

        pathway_fname = f'results/pathways/csv_files/{pathway}.csv'
        data_wt_fname = f'resources/data/{study}/Counts_Ctrl_{treatment}.csv'
        data_mt_fname = f'resources/data/{study}/Counts_{gene}_{treatment}.csv'

        # read data
        graph = nx.from_pandas_edgelist(pd.read_csv(pathway_fname), 'source', 'sink')
        df_wt = pd.read_csv(data_wt_fname, index_col=0)
        df_mt = pd.read_csv(data_mt_fname, index_col=0)
        df = pd.concat([df_wt, df_mt], axis=1)

        # compute statistics
        for g in gene_list:
            tmp.extend(
                [
                    {
                        'type': 'pathway_degree',
                        'gene': gene,
                        'source': f'{pathway} -- {g}',
                        'value': graph.degree[g] if g in graph else np.nan,
                    },
                    {
                        'type': 'mean_expression_count',
                        'gene': gene,
                        'source': f'{treatment} -- {g}',
                        'value': df.loc[g].mean(),
                    },
                ]
            )

    return pd.DataFrame(tmp)


def main(dir_list, fname, out_dir):
    # prepare environment
    out_dir.mkdir(parents=True, exist_ok=True)

    # read data
    df = read_data(dir_list)
    df.to_csv(fname, index=False)

    # create plots
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=df[df['type'] == 'pathway_degree'], x='gene', y='value')
    plt.xlabel('Perturbed gene(s)')
    plt.ylabel('Pathway degree')
    plt.tight_layout()
    plt.savefig(out_dir / 'degrees.pdf')

    plt.figure(figsize=(8, 6))
    sns.boxplot(
        data=df[df['type'] == 'mean_expression_count'].assign(
            treatment=lambda x: x['source'].str.split(' -- ').str[0]
        ),
        x='gene',
        y='value',
        hue='treatment',
    )
    plt.xlabel('Perturbed gene(s)')
    plt.ylabel('Mean expression')
    plt.tight_layout()
    plt.savefig(out_dir / 'counts.pdf')


if __name__ == '__main__':
    main(
        snakemake.input.dir_list, snakemake.output.fname, Path(snakemake.output.out_dir)
    )
