import pandas as pd
import os
from scipy import stats


def set_node_attrs(edges, diff_info, RPGs, HISTs, TFs, edge_column='gene_symbol', out=None):
    if out is None:
        out = edges + '.nodes'
    edges = pd.read_csv(edges, header=0, sep=None, engine='python', index_col=False)
    nodes = set(edges.iloc[:, 0]) | set(edges.iloc[:, 1])
    print('node number:', len(nodes))
    data = pd.read_csv(diff_info, header=0, sep=None, engine='python')
    if edge_column not in data.columns:
        raise Exception(f'edge_column must be in header of {data.columns}')
    data = data.set_index(edge_column)
    dup_ind = data.index.duplicated()
    if sum(dup_ind) > 1:
        print(data.index[dup_ind], 'are duplicated', 'and keep first one')
        data = data.loc[~dup_ind]
    rpgs = [x.strip().split()[0] for x in open(RPGs)]
    hists = [x.strip().split()[0] for x in open(HISTs)]
    tfs = [x.strip().split()[0] for x in open(TFs)]
    with open(out, 'w') as f:
        f.write('node\ttype\t{}\n'.format('\t'.join(data.columns)))
        for s in nodes:
            s_type = 'Gene'
            if s in rpgs:
                s_type = 'RPG'
            if s in hists:
                s_type = 'HIST'
            if s in tfs:
                s_type = 'TF'
            f.write('{}\t{}\t{}\n'.format(s, s_type, '\t'.join(data.loc[s].apply(lambda x:str(x)))))


def pair_corr(exp, pair_info, reg_direct=None, pval_cutoff=0.05, corr_cutoff=0.3, out=None):
    direct_dict = dict()
    if reg_direct is not None:
        direct = pd.read_csv(reg_direct, header=0, sep=None, engine='python')
        direct_dict = dict(zip(zip(direct.iloc[:, 0], direct.iloc[:, 3]), direct.iloc[:, 5]))
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
        if direct_dict:
            f.write('source\ttarget\tcorrelation\tcorr_type\tpvalue\tregulation\n')
        else:
            f.write('source\ttarget\tcorrelation\tcorr_type\tpvalue\n')
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
                        if reg_direct:
                            regulation = 'unknown'
                            if (p1, p2) in direct_dict:
                                regulation = direct_dict[(p1, p2)]
                                if regulation == '-->':
                                    regulation = 'up'
                                else:
                                    'down'
                            elif (p2, p1) in direct_dict:
                                regulation = direct_dict[(p2, p1)]
                                if regulation == '-->':
                                    regulation = 'down'
                                else:
                                    'up'
                            f.write(f'{p1}\t{p2}\t{corr}\t{direct}\t{pval}\t{regulation}\n')
                        else:
                            f.write(f'{p1}\t{p2}\t{corr}\t{direct}\t{pval}\n')
                except Exception as e:
                    print(e)
                    print(p1, p2)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['set_node_attrs', 'pair_corr'])

