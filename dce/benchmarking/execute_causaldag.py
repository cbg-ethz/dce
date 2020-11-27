import sys
import itertools

import numpy as np
import pandas as pd

import causaldag


def main(X_wt, X_mt, fname_out):
    # transform data
    X_wt_trans = np.log(X_wt + 1).to_numpy()
    X_mt_trans = np.log(X_mt + 1).to_numpy()

    # inference
    p = X_wt_trans.shape[1]

    # difference_matrix = causaldag.dci(
    difference_matrix, _ = causaldag.dci_stability_selection(
        X_wt_trans, X_mt_trans,
        difference_ug=list(itertools.combinations(range(p), 2))
    )

    # extract differential edges
    # ddag_edges = set(zip(*np.where(difference_matrix != 0)))

    m = pd.DataFrame(difference_matrix, index=X_wt.columns, columns=X_wt.columns)

    m.to_csv(fname_out)


if __name__ == '__main__':
    main(
        pd.read_csv(sys.argv[1], index_col=0),
        pd.read_csv(sys.argv[2], index_col=0),
        sys.argv[3]
    )
