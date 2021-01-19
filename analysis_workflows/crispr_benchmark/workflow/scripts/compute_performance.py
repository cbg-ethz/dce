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


def compute_roc(class_prob, true_class):
    fpr_list, tpr_list, roc_thresholds = roc_curve(true_class, class_prob)
    roc_auc = auc(fpr_list, tpr_list)

    return fpr_list, tpr_list, roc_auc


def main(fname, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(fname)
    # df['true_effect'] = (df['perturbed_gene'] == df['source']) | (df['perturbed_gene'] == df['target'])
    df['true_effect'] = df.apply(lambda x: x['source'] in x['perturbed_gene'].split(',') or x['target'] in x['perturbed_gene'].split(','), axis=1)

    plot_dir = out_dir / 'plots'
    plot_dir.mkdir(parents=True, exist_ok=True)

    tmp = []
    for pathway, group_pathway in tqdm(df.groupby('pathway')):
        s = 2
        fig_roc, ax_roc = plt.subplots(figsize=(s * 8, s * 6))
        fig_pr, ax_pr = plt.subplots(figsize=(s * 8, s * 6))

        for idx, group in tqdm(group_pathway.groupby([
            'study', 'treatment', 'perturbed_gene'
        ]), leave=False):
            study, treatment, perturbed_gene = idx

            # TODO: find better way of handling NAs
            group.dropna(inplace=True)

            # sanity checks
            if group['true_effect'].sum() == 0:
                print(f'[{pathway}] Skipping {idx}, no true positives...')
                continue

            # compute performance measures
            #edge_score = group['dce_pvalue']
            edge_score = abs(group['dce'])

            precision_list, recall_list, pr_thresholds = precision_recall_curve(group['true_effect'], edge_score)
            fpr_list, tpr_list, roc_thresholds = roc_curve(group['true_effect'], edge_score)

            roc_auc = auc(fpr_list, tpr_list)
            pr_auc = auc(recall_list, precision_list)

            ap_score = average_precision_score(group['true_effect'], edge_score)

            # choose some "optimal" threshold
            thres_opt = roc_thresholds[abs(tpr_list - (1 - fpr_list)).argmin()]

            f1_score_ = f1_score(group['true_effect'], edge_score >= thres_opt)

            # plots
            app = '_'.join(str(e) for e in idx)

            plt.figure()
            cm = confusion_matrix(group['true_effect'], edge_score >= thres_opt)
            sns.heatmap(cm, annot=True, fmt='d', square=True)
            plt.xlabel('Predicted label')
            plt.ylabel('True label')
            plt.tight_layout()
            plt.savefig(plot_dir / f'confusion_matrix_{pathway}_{app}.pdf')

            base_line, = ax_roc.plot(
                fpr_list, tpr_list,
                '-o',
                label=f'{app} ({roc_auc:.2})')
            ax_pr.plot(
                recall_list, precision_list,
                '-o',
                label=f'{app} ({pr_auc:.2})')

            # competing methods
            method_columns = ['cor', 'pcor']

            for method in method_columns:
                fpr, tpr, auc_val = compute_roc(group[method], group['true_effect'])

                cur_color = base_line.get_color()
                ax_roc.plot(
                    fpr, tpr,
                    's:',
                    label=f'[{method}] {app} ({auc_val:.2})',
                    color=cur_color)

                tmp.append({
                    'method': method,

                    'pathway': pathway,
                    'perturbed_gene': perturbed_gene,
                    'study': study,
                    'treatment': treatment,

                    'roc_auc': auc_val
                })

            # store results
            tmp.append({
                'method': 'dce',

                'pathway': pathway,
                'perturbed_gene': perturbed_gene,
                'study': study,
                'treatment': treatment,

                'roc_auc': roc_auc,
                'pr_auc': pr_auc,
                'ap_score': ap_score,
                'f1_score': f1_score_,

                'optimal_roc_threshold': thres_opt
            })

        # finalize figures
        ax_roc.plot([0, 1], [0, 1], color='grey', ls='dashed')
        ax_roc.set_xlabel('FPR')
        ax_roc.set_ylabel('TPR')
        ax_roc.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        fig_roc.tight_layout()
        fig_roc.savefig(plot_dir / f'roc_curce_{pathway}.pdf')

        ax_pr.set_xlabel('Recall')
        ax_pr.set_ylabel('Precision')
        ax_pr.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        fig_pr.tight_layout()
        fig_pr.savefig(plot_dir / f'pr_curve_{pathway}.pdf')

        plt.close('all')

    df_perf = pd.DataFrame(tmp)
    df_perf.to_csv(out_dir / 'measures.csv', index=False)


if __name__ == '__main__':
    main(snakemake.input.fname, Path(snakemake.output.out_dir))
