import pandas as pd


def top_rank_gene(exp_matrix, top=(2000., 0.75), index_col=(0, 1), header=0):
    """
    挑选在N个样本中表达量排名均在前T的基因, N=top[1], T=top[0]
    :param exp_matrix:
    :param top: 需要给两个值, 空格分隔即可. 若 top[0]<=1，则根据比例计算, 否则直接用该数;第二个数top[1]与第一个同理。
    :param index_col: 行索引号, 可以指定多列
    :param header: 列索引号, 可以指定多行
    :return:
    """
    exp = pd.read_csv(exp_matrix, index_col=index_col, header=header,
                      sep=None, engine='python')
    exp = exp.loc[exp.sum(axis=1) > 0]
    rexp = exp.rank(ascending=False)
    select = list(top)

    if float(top[0]) <= 1:
        select[0] = int(exp.shape[0]*top[0])
    if float(top[1]) <= 1:
        select[1] = int(exp.shape[1]*top[1])
    select = [int(x) for x in select]
    ind = rexp.apply(lambda x: sum(y <=select[0] for y in x) >= select[1], axis=1)
    top_exp = exp.loc[ind]
    out_name = 'top{}.exp.{}_in_{}samples.csv'.format(top_exp.shape[0], select[0], select[1])
    top_exp.to_csv(out_name)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['top_rank_gene'])
