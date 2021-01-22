import os
import GEOparse
import pandas as pd


def get_expr_table(gse_id_or_path, outdir='./'):
    if os.path.exists(gse_id_or_path):
        gse = GEOparse.parse_GSE(gse_id_or_path, {})
    else:
        gse = GEOparse.get_GEO(geo=gse_id_or_path, destdir=outdir, include_data=True, annotate_gpl=True)
        print(f'success to download {gse_id_or_path}')

    # warn info
    print(f'there are {len(gse.gsms)} samples in {gse.name}')
    print(f'there are multiple platforms {gse.gpls.keys()} in {gse.name}')

    # get sample data
    gpl_tables = dict()
    sample_info_lst = []
    for sample, sample_obj in gse.gsms.items():
        platform = sample_obj.metadata['platform_id'][0]
        characteristics = sample_obj.metadata.get('characteristics_ch1', 'None')
        sample_title = sample_obj.metadata['title'][0]
        sample_type = sample_obj.metadata['type'][0]
        sample_desc = sample_obj.metadata.get('description', ['None'])[0]
        sample_process = sample_obj.metadata.get('data_processing', ['None'])[0]
        sample_info_lst.append([sample, gse.name, platform, sample_type,
                                characteristics, sample_title, sample_desc, sample_process])
        first_col = sample_obj.table.columns[0]
        table = sample_obj.table.set_index(first_col)
        table.columns = [sample]
        if platform not in gpl_tables:
            gpl_tables[platform] = table
        else:
            gpl_tables[platform] = gpl_tables[platform].join(table)

    with open(f'{gse.name}.sample_info.txt', 'w') as f:
        f.write(f'gsm_id\tgse_id\tplatform\ttype\tcharacteristics\ttitle\tdescription\tdata_processing\n')
        for each in sample_info_lst:
            f.write('\t'.join(str(x) for x in each)+'\n')

    # annotate data
    for gpl_id, data in gpl_tables.items():
        annot_table = gse.gpls[gpl_id].table.set_index('ID')
        data = data.join(annot_table)
        data.to_csv(f'{gse.name}.{gpl_id}.expr.table.txt', sep='\t')




if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())



