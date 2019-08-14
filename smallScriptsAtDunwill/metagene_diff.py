import os
import pandas as pd
import scipy.stats as stats
from functools import partial
import statistics


def read_metagene_group(group_txt, format='gene2group'):
    if format != 'gene2group':
        result_dict = dict()
        with open(group_txt) as f:
            for line in f:
                metagene, genes = line.strip().split('\t', 1)
                genes = [y.strip() for x in genes.split() for y in x.strip().split(',')]
                if metagene in result_dict:
                    raise Exception(f"{metagene} duplicated, please check it!")
                else:
                    result_dict[metagene] = genes
    else:
        result_dict = read_sample_group(group_txt)
    return result_dict


def read_sample_group(group_txt):
    group_df = pd.read_csv(group_txt, sep=None, engine='python', header=0, index_col=0)
    group_dict = dict()
    for scheme in group_df:
        tmp_dict = dict(list(group_df.loc[:, [scheme]].groupby(scheme)))
        for group, df_val in tmp_dict.items():
            if df_val.shape[0] == group_df.shape[0]:
                raise Exception('In column of {}, only one group was found!'.format(scheme))
            group_dict[group] = sorted(df_val.index)
    return group_dict


def read_compare_info(cmp_info, group_dict):
    with open(cmp_info) as f:
        cmp_list = list()
        error_names = list()
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            tmp_ctrl, tmp_test = line.strip().split()
            if tmp_ctrl not in group_dict:
                error_names.append(tmp_ctrl)
            if tmp_test not in group_dict:
                error_names.append(tmp_test)
            cmp_list.append((tmp_ctrl, tmp_test))
        if error_names:
            print(f'Please be aware that Each group name of {error_names} is not in group dict!')
    cmp_list = sorted(list(set(cmp_list)))
    return cmp_list


def score_metagene(exp_matrix, metagene_dict:dict, score='geometric_mean'):
    exp_df = pd.read_csv(exp_matrix, header=0, index_col=0, sep=None, engine='python')
    score_df = pd.DataFrame(columns=exp_df.columns)
    valid_member_list = list()
    for metagene, genes in metagene_dict.items():
        intersect = set(genes) & set(exp_df.index)
        if intersect:
            valid_members = ','.join(intersect)
            valid_member_list.append(valid_members)
            meta_exp = exp_df.loc[intersect]
            if score == 'geometric_mean':
                meta_score = meta_exp.apply(stats.gmean, axis=0)
            else:
                meta_score = meta_exp.mean(axis=0)
            score_df.loc[metagene] = meta_score
        else:
            print(f'None member of metagene "{metagene}" is in expression matrix')
    score_df['members'] = valid_member_list
    score_df.set_index('members', append=True, inplace=True)
    return score_df


def diff_test(score_df, group_dict, cmp_list, method='mannwhitneyu', equal_var=True, prefix=''):
    for ctrl, test in cmp_list:
        if not(ctrl in group_dict and test in group_dict):
            print(f'skip {ctrl} and {test} for both or one of them are/is not in group info dict!')
            continue
        ctrl_samples = group_dict[ctrl]
        test_samples = group_dict[test]
        ctrl_num = len(ctrl_samples)
        target_data = score_df[ctrl_samples+test_samples]
        if method == 'ranksums':
            test_func = stats.ranksums
        elif method == 'mannwhitneyu':
            test_func = stats.mannwhitneyu
        elif method == "wilcoxon":
            test_func = stats.wilcoxon
        elif method == 'ttest_ind':
            if equal_var:
                test_func = stats.ttest_ind
            else:
                test_func = partial(stats.ttest_ind, equal_var=False)
        test_df = pd.DataFrame()
        test_df['pvalue'] = target_data.apply(lambda x:test_func(x[:ctrl_num], x[ctrl_num:])[1], axis=1)
        ctrl_median_exp = target_data[ctrl_samples].apply(statistics.median, axis=1)
        test_median_exp = target_data[test_samples].apply(statistics.median, axis=1)
        test_df[ctrl+'_median'] = ctrl_median_exp
        test_df[test+'_median'] = test_median_exp
        test_df['median_log2FC'] = test_median_exp - ctrl_median_exp
        test_df['regulation'] = 'down'
        test_df.loc[test_df['median_log2FC']>0, 'regulation'] = 'up'
        test_df = test_df.join(target_data)
        test_df.sort_values(by='pvalue', inplace=True)
        test_df.to_csv(f'{prefix}{ctrl}_vs_{test}.metagene.{method}.xls', sep='\t')


def metagene_diff(exp_matrix, sample_group, metagene_group, compare, prefix='',
                  metagene_group_format='gene2group', score='geometric_mean',
                  method='mannwhitneyu', equal_var=True):
    metagene_group_dict = read_metagene_group(metagene_group, format=metagene_group_format)
    metagene_score_df = score_metagene(exp_matrix, metagene_group_dict, score=score)
    sample_group_dict = read_sample_group(sample_group)
    compare_list = read_compare_info(compare, sample_group_dict)
    diff_test(metagene_score_df, sample_group_dict, compare_list,
              method=method, equal_var=equal_var, prefix=prefix)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['metagene_diff'])

