# coding=utf-8
import pandas as pd


def exp_calculator_with_count(count_table_file):
    """
    calculate fpkm and tpm based on count table with second column containing gene length.
    :param count_table_file: example:
    -----------
    gene_id gene_length sample1 sample2
    gene1   1001    29  50
    gene2   1300    30  14
    -----------
    :return: rpkm_dict, tpm_dict
    """
    count_table = pd.read_table(count_table_file, index_col=0)
    columns = count_table.columns

    gene_len = count_table[columns[0]]
    rpkm_dict = dict()
    tpm_dict = dict()
    for sample in columns[1:]:
        # Divide the read counts by the length of each gene in kilobases.
        # This gives you reads per kilobase (RPK)
        rpk = count_table[sample]/gene_len
        # get rpkm/fpkm
        total_counts = sum(count_table[sample])
        rpkm = rpk/total_counts*1000000*1000
        # get tpm
        norm_gene_len_total_counts = sum(rpk)
        tpm = rpk/norm_gene_len_total_counts*1000000
        # save
        rpkm_dict[sample] = rpkm
        tpm_dict[sample] = tpm
    # save results
    df_rpkm = pd.DataFrame(rpkm_dict, index=count_table.index)
    df_tpm = pd.DataFrame(tpm_dict, index=count_table.index)
    df_rpkm.to_csv(count_table_file+'.fpkm.xls', sep='\t')
    df_tpm.to_csv(count_table_file+'.tpm.xls', sep='\t')
    #
    return rpkm_dict, tpm_dict

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="count table file, example:\n"
                                                 "------------------------------------\n"
                                                 "gene_id\tgene_length\tsample1\tsample2\n"
                                                 "gene1\t1001\t29\t50\n"
                                                 "gene2\t1300\t30\t14\n"
                                                 "... ... ... ...\n"
                                                 "-----------------------------------\n")
    parser.add_argument('-count', type=str, required=True, metavar="count_table",
                        help="count table file with second column containing gene length")
    args = parser.parse_args()
    exp_calculator_with_count(args.count)


