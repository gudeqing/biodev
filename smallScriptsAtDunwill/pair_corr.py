from scipy import stats
import pandas as pd
import os


def pair_corr(exp, pair_info, pval_cutoff=0.05, corr_cutoff=0.3, out=None):
    pairs = [x.strip().split()[:2] for x in open(pair_info)]
    data = pd.read_csv(exp, index_col=0, header=0, sep=None, engine='python')
    dup_ind = data.index.duplicated()
    if sum(dup_ind) > 1:
        print(data.index[dup_ind], 'are duplicated', 'and keep first one')
        data = data.loc[~dup_ind]
    if out is None:
        out = os.path.basename(pair_info)
        out = f'P{pval_cutoff}.C{corr_cutoff}.' + out
    with open(out, 'w') as f:
        for p1, p2 in pairs:
            judge = p1 in data.index
            judge2 = p2 in data.index
            if not judge:
                print(p1, 'not in exp')
            if not judge2:
                print(p2, 'not in exp')
            if judge and judge2:
                try:
                    corr, pval = stats.pearsonr(data.loc[p1], data.loc[p2])
                    direct = 'positive' if corr > 0 else 'negative'
                    if abs(corr) >= corr_cutoff and pval <= pval_cutoff:
                        f.write(f'{p1}\t{p2}\t{corr}\t{direct}\t{pval}\n')
                except Exception as e:
                    print(e)
                    print(p1, p2)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pair_corr'])
