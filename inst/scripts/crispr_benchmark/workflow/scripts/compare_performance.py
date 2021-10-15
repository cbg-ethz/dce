"""
Compare DCE estimation performance (ROC-AUC) between different
deconfounding methods.
"""


from pathlib import Path

import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

from bioinf_common.plotting import add_identity


def main(dname_list, outdir):
    outdir.mkdir(parents=True, exist_ok=True)

    # read data
    df_list = []
    for dname in dname_list:
        tmp = pd.read_csv(dname / 'measures.csv')

        deconf_param = str(dname).split('/')[-2].split('~')[1]
        tmp['deconf_param'] = deconf_param

        df_list.append(tmp)
    df = pd.concat(df_list)

    # transform data
    df_prep = df.pivot(index=['pathway', 'perturbed_gene', 'study', 'treatment', 'method'], columns=['deconf_param'], values=['roc_auc'])['roc_auc'].reset_index('method')

    # plot data
    g = sns.pairplot(data=df_prep, hue='method')

    g.map_offdiag(lambda *args, **kwargs: add_identity(plt.gca(), color='grey', ls='dashed'))
    g.map_offdiag(lambda *args, **kwargs: plt.xlim(0, 1))
    g.map_offdiag(lambda *args, **kwargs: plt.ylim(0, 1))

    g.savefig(outdir / 'pairplot.pdf')


if __name__ == '__main__':
    main([Path(d) for d in snakemake.input.dname_list], Path(snakemake.output.dname))
