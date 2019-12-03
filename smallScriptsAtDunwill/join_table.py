#! /data/users/dqgu/anaconda3/bin/python
import pandas as pd


def merge_data(files:list, index_cols:list=None, new_col_name=None, out='merged.table.xls', how='left'):
    """
    使用pandas的join函数合并表格, 要求所有table的第一行为header, 注意：如果第一列对应的header为空，则命名为index
    :param files: table路径, 如table1 table2 table3
    :param index_cols: 指定一列或几列作为索引，后续合并将基于指定的索引进行合并, 默认按行号合并。
    :param new_col_name: 至少两列，前两列用来构建字典，对合并后的表格的列名进行重新命名。
    :param out: 合并后的输出文件名
    :param how: How to handle the operation of the two objects.
        left: use calling frame’s index (or column if on is specified)
        right: use other’s index. 即每次合并以最后一个表格的index为主
        outer: form union of calling frame’s index with other's index and sort it. lexicographically.
        inner: form intersection of calling frame’s index with other's index, preserving the order of the calling’s one.
    :return: 输出合并后的表格
    """
    # index_cols = ('Target', 'Gene', 'ENTREZ_GENE_ID', 'NCBI_NAME', 'NCBI_ACCESSION', 'GENE_FUNCTION')
    table = pd.read_csv(files[0], index_col=None, header=0, sep=None, engine='python')
    if table.columns[0] == 'Unnamed: 0':
        table.columns = ['index'] + list(table.columns[1:])
    # print(table.head())
    if index_cols:
        table.set_index(index_cols, inplace=True)
    for each in files[1:]:
        each_table = pd.read_csv(each, index_col=None, header=0, sep=None, engine='python')
        if each_table.columns[0] == 'Unnamed: 0':
            each_table.columns = ['index'] + list(each_table.columns[1:])
        if index_cols:
            each_table.set_index(index_cols, inplace=True)
        table = table.join(each_table, how=how)
    table.columns = [x.strip() for x in table.columns]
    if new_col_name is not None:
        new_name_df = pd.read_csv(new_col_name, index_col=None, header=None, sep=None, engine='python')
        new_name_dict = dict(zip(new_name_df.iloc[0], new_name_df.iloc[1]))
        table.columns = [new_name_dict[x] for x in table.columns]
    out_index = True
    if not index_cols:
        out_index = False
    if not out.endswith('.csv'):
        table.to_csv(out, sep='\t', index=out_index)
    else:
        table.to_csv(out, out_index)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['merge_data'])

