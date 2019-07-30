import pandas as pd


def table_knife(table, index_cols:list=None, header:list=[0],
                target_cols=None, target_rows=None,
                sum_table=False, describe_table=False,
                sort_by:list=None, axis=0, descending=False,
                new_col_name=None, new_row_name=None,
                out_name='new_table.xls'):
    """
    功能: 对数据表进行切割, 排序, 求和, 统计, 修改列名或行名
    :param table: 数据表路径
    :param index_cols: 指定哪一列或哪些列作为行索引, 默认用行的序号作为索引，如果第一列没有header，默认用"index"作为列名
    :param header: 指定那一行或哪些行作为表头header
    :param target_cols: 文件路径，包含一列，指定要提取的列, 当header为多列时，无效
    :param target_rows: 文件路径，包含一列，指定要提取的行, 当行索引为多列时，无效
    :param sum_table: 打印出每一列的和或每一行的和, 默认统计行
    :param describe_table: 打印统计表格信息
    :param sort_by: 指定要排序依据的列, 可以是多列
    :param axis: 指定按行处理(=0) 还是按列处理(=1)
    :param descending: 若使用，则降序
    :param new_col_name: 文件路径，两列，第一列为旧的列名，第二列为新的列名，tab分割, 无header, 当header为多列时，无效
    :param new_row_name: 文件路径，两列，第一列为旧的行名，第二列为新的行名，tab分割, 无header, 当行索引为多列时，无效
    :param out_name: 新的表名
    :return: None
    """
    table = pd.read_csv(table, index_col=None, header=header, sep=None, engine='python')
    if table.columns[0] == 'Unnamed: 0':
        table.columns = ['index'] + list(table.columns[1:])
    if index_cols is not None:
        table.set_index(index_cols, inplace=True)
    if target_cols:
        table = table[[x.strip() for x in open(target_cols)]]
    if target_rows:
        table = table.loc[[x.strip() for x in open(target_rows)]]
    if new_col_name:
        new_name_dict = dict(x.strip().split('\t')[:2] for x in open(new_col_name))
        table.columns = [new_name_dict[x] if x in new_name_dict else x for x in table.columns]
    if new_row_name:
        new_name_dict = dict(x.strip().split('\t')[:2] for x in open(new_row_name))
        table.index = [new_name_dict[x] if x in new_name_dict else x for x in table.index]
    if sum_table:
        table = table.sum(axis=axis)
        print(table)
    if describe_table:
        table = table.describe()
        print(table)
    if sort_by:
        table.sort_values(by=sort_by, inplace=True, axis=axis, ascending=not descending)
    if out_name.endswith('.csv'):
        table.to_csv(out_name)
    else:
        table.to_csv(out_name, sep='\t')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['table_knife'])
