from pathlib import Path

import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.metrics import (
    precision_recall_curve, roc_curve, auc,
    average_precision_score, f1_score,
    confusion_matrix
)

from tqdm import tqdm


def main(fname, out_dir):
    df = pd.read_csv(fname)
    # df['true_effect'] = (df['perturbed_gene'] == df['source']) | (df['perturbed_gene'] == df['target'])
    df['true_effect'] = df.apply(lambda x: x['source'] in x['perturbed_gene'].split(',') or x['target'] in x['perturbed_gene'].split(','), axis=1)

    plot_dir = out_dir / 'plots'
    plot_dir.mkdir()

    tmp = []
    for idx, group in tqdm(df.groupby([
        'study', 'treatment', 'perturbed_gene', 'pathway'
    ])):
        # sanity checks
        if group['true_effect'].sum() == 0:
            print(f'Skipping {idx}, no true positives...')
            continue

        # TODO: find better way of handling NAs
        group.dropna(inplace=True)

        # compute performance measures
        precision_list, recall_list, pr_thresholds = precision_recall_curve(group['true_effect'], group['dce_pvalue'])
        fpr_list, tpr_list, roc_thresholds = roc_curve(group['true_effect'], group['dce_pvalue'])

        roc_auc = auc(fpr_list, tpr_list)
        pr_auc = auc(recall_list, precision_list)

        ap_score = average_precision_score(group['true_effect'], group['dce_pvalue'])

        f1_score_ = f1_score(group['true_effect'], group['dce_pvalue'] < .05)

        # plots
        app = '_'.join(str(e) for e in idx)

        plt.figure()
        cm = confusion_matrix(group['true_effect'], group['dce_pvalue'] < .05)
        sns.heatmap(cm, annot=True, fmt='d', square=True)
        plt.xlabel('Predicted label')
        plt.ylabel('True label')
        plt.tight_layout()
        plt.savefig(plot_dir / f'confusion_matrix_{app}.pdf')

        plt.figure(figsize=(8, 6))
        plt.plot(fpr_list, tpr_list)
        plt.plot([0, 1], [0, 1], color='grey', ls='dashed')
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.title(f'ROC-AUC: {roc_auc:.2}')
        plt.tight_layout()
        plt.savefig(plot_dir / f'roc_curce_{app}.pdf')

        plt.figure(figsize=(8, 6))
        plt.plot(recall_list, precision_list)
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title(f'PR-AUC: {pr_auc:.2}')
        plt.tight_layout()
        plt.savefig(plot_dir / f'pr_curve_{app}.pdf')

        plt.close('all')

        # store results
        study, treatment, perturbed_gene, pathway = idx
        tmp.append({
            'study': study,
            'treatment': treatment,
            'perturbed_gene': perturbed_gene,
            'pathway': pathway,

            'roc_auc': roc_auc,
            'pr_auc': pr_auc,
            'ap_score': ap_score,
            'f1_score': f1_score_

        })

    df_perf = pd.DataFrame(tmp)
    df_perf.to_csv(out_dir / 'measures.csv', index=False)


if __name__ == '__main__':
    main(snakemake.input.fname, Path(snakemake.output.out_dir))
