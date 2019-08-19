#! /data/users/dqgu/anaconda3/bin/python
import pandas as pd


def merge_data(files:list, index_cols:list=None, new_col_name=None, out='merged.table.xls', how='outer'):
    if index_cols is None:
        index_cols = ('Target', 'Gene', 'ENTREZ_GENE_ID', 'NCBI_NAME', 'NCBI_ACCESSION', 'GENE_FUNCTION')
    table = pd.read_csv(files[0], index_col=None, header=0, sep=None, engine='python')
    print(table.head())
    table.set_index(index_cols, inplace=True)
    for each in files[1:]:
        each_table = pd.read_csv(each, index_col=None, header=0, sep=None, engine='python')
        each_table.set_index(index_cols, inplace=True)
        table = table.join(each_table, how=how)
    table.columns = [x.strip() for x in table.columns]
    if new_col_name is not None:
        new_name_df = pd.read_csv(new_col_name, index_col=None, header=0, sep=None, engine='python')
        new_name_dict = dict(zip(new_name_df.iloc[0], new_name_df.iloc[1]))
        table.columns = [new_name_dict[x] for x in table.columns]
    if not out.endswith('.csv'):
        table.to_csv(out, sep='\t')
    else:
        table.to_csv(out)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['merge_data'])

