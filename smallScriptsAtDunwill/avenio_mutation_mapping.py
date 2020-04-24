import os
import re
import logging
import pysam
import statistics
import numpy as np
import pandas as pd
from pysam import VariantFile, VariantHeader
import subprocess


def cross_map(vcf, chain=None, ref=None):
    chain = chain if chain else '/nfs2/database/CoordinateConvertion/hg38ToHg19.over.chain.gz'
    ref = ref if ref else '/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta'
    outname = 'hg19.' + os.path.basename(vcf)
    outdir = os.path.dirname(vcf)
    out = os.path.join(outdir, outname)
    if not os.path.exists(out):
        cmd = f'CrossMap.py vcf {chain} {vcf} {ref} {out}'
        subprocess.check_call(cmd, shell=True)
    else:
        print('mapped result already found !')
    return out


def extract_hot(vcfs:list, avenio_table, hotspot, chain=None, ref=None, out='in_hotspot.xlsx'):
    chain = chain if chain else '/nfs2/database/CoordinateConvertion/hg38ToHg19.over.chain.gz'
    ref = ref if ref else '/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta'
    mapping = []
    hg19vcfs = []
    for vcf in vcfs:
        if os.path.getsize(vcf) <= 2:
            print(vcf, 'is empty and skipped!')
            continue
        print('processing', vcf)
        hg38_id_lst = []
        with VariantFile(vcf) as f:
            for record in f:
                sample = list(f.header.samples)[0]
                af = round(record.info['AF'][0], 4)
                mut_id = f'{sample}:{record.chrom}:{record.pos}:{af:.2%}'
                hg38_id_lst.append(mut_id)

        hg19_id_lst = []
        vcf = cross_map(vcf, chain, ref)
        hg19vcfs.append(vcf)
        with VariantFile(vcf) as f:
            for record in f:
                sample = list(f.header.samples)[0]
                af = round(record.info['AF'][0], 4)
                mut_id = f'{sample}:{record.chrom}:{record.pos}:{af:.2%}'
                hg19_id_lst.append(mut_id)
        assert len(hg19_id_lst) == len(hg38_id_lst)
        for one, two in zip(hg38_id_lst, hg19_id_lst):
            mapping.append([one, two])
    # check
    for each in mapping:
        if mapping.count(each) > 1:
            print(f'Warning: mut_id {each} is not uniq')

    converter = dict(mapping)
    with open('converter.txt', 'w') as f:
        for k, v in converter.items():
            f.write(f'{k}\t{v}\n')

    # extract hotspot
    cmd = 'python /data/users/dqgu/PycharmProjects/biodev/smallScriptsAtDunwill/validation.py '
    cmd += 'batch_extract_hotspot '
    cmd += '-vcfs '
    for vcf in hg19vcfs:
        cmd += f'{vcf} '
    cmd += f'-hotspot {hotspot} '
    cmd += '-id_mode chr:start:ref:alt --af_in_info --dp_in_info -sample_index 0'
    subprocess.check_call(cmd, shell=True)
    hot_df = pd.read_csv('all.detected.hotspot.xls', header=0, sep='\t')
    hot_df['hg19siteID'] = hot_df['sample'] + ':' + hot_df['site'] + ':' + hot_df['AF(%)']
    hot_df = hot_df.set_index('hg19siteID')
    hot_df = hot_df[['mutation', 'ID']]
    print(f'found {hot_df.shape[0]} hotspot !')

    #
    avenio = pd.read_csv(avenio_table, header=0, sep='\t')
    formated_af = []
    for af in avenio['Allele Fraction']:
        if type(af) == float:
            formated_af.append(f'{af:.2%}')
        else:
            if af.endswith('%'):
                af = float(af[:-1])/100
            else:
                try:
                    af = float(af)
                except:
                    print(f'warn: failed to convert {af} to float')
                    af = 0.0
            formated_af.append(f'{af:.2%}')
    avenio['Allele Fraction'] = formated_af
    hg38_ids = avenio['Sample ID']+':'+avenio['Genomic Position']+':'+avenio['Allele Fraction']
    avenio['hg19siteID'] = [converter[x] if x in converter else x for x in hg38_ids]
    avenio = avenio.set_index('hg19siteID')
    result = hot_df.join(avenio)
    result.to_excel(out)


def match_avenio_hot(hot_table, avenio_table):
    avenio = pd.read_csv(avenio_table, header=0, sep='\t')
    hot = pd.read_csv(hot_table, header=0, sep='\t')
    hit_inds = []
    for ind, row in avenio.iterrows():
        gene = row['Gene']
        phgvs = row['Amino Acid Change']
        for i, r in hot.iterrows():
            try:
                if gene == r['Gene'] and r['pHgvs'].replace('(', '').replace(')', '') in phgvs:
                    if ind in hit_inds:
                        continue
                    else:
                        hit_inds.append(ind)
            except:
                print(r['Gene'], r['pHgvs'], phgvs)
    result = avenio.loc[hit_inds]
    result.to_excel('in_avenio_hot.xlsx')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['cross_map', 'extract_hot', 'match_avenio_hot'])
