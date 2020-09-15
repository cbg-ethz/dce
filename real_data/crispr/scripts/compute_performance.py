from pathlib import Path

import pandas as pd

from sklearn.metrics import (
    precision_recall_curve, roc_curve, auc,
    average_precision_score, f1_score
)


def main(fname, out_dir):
    df = pd.read_csv(fname)
    # df['true_effect'] = (df['perturbed_gene'] == df['source']) | (df['perturbed_gene'] == df['target'])
    df['true_effect'] = df.apply(lambda x: x['source'] in x['perturbed_gene'].split(',') or x['target'] in x['perturbed_gene'].split(','), axis=1)

    tmp = []
    for idx, group in df.groupby([
        'study', 'treatment', 'perturbed_gene', 'pathway'
    ]):
        # sanity checks
        if group['true_effect'].sum() == 0:
            print(f'Skipping {idx}, no true positives...')
            continue

        # compute performance measures
        precision_list, recall_list, pr_thresholds = precision_recall_curve(group['true_effect'], group['dce_pvalue'])
        fpr_list, tpr_list , roc_thresholds = roc_curve(group['true_effect'], group['dce_pvalue'])

        roc_auc = auc(fpr_list, tpr_list)
        pr_auc = auc(recall_list, precision_list)

        ap_score = average_precision_score(group['true_effect'], group['dce_pvalue'])

        f1_score_ = f1_score(group['true_effect'], group['dce_pvalue'] < .05)

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
