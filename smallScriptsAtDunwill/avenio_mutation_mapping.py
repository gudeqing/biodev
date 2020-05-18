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


def aa_three2one(phgvs):
    if type(phgvs) != str:
        return phgvs
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    aad = {x.lower().capitalize(): y for x, y in d.items()}
    for k, v in aad.items():
        phgvs = phgvs.replace(k, v)
    return phgvs


def annotation(bed, comm_trans_file='/nfs2/database/CommonRefSeqTranscripts/canonical_ensgene_hg19_v90.txt', gtf="/nfs2/database/gencode_v29/gencode.v29.primary_assembly.annotation.gtf"):
    annot_out = f'{bed}.gtfanno'
    cmd = f"bedtools intersect -a {gtf}  -b {bed} -wo > {annot_out}"
    print(cmd)
    os.system(cmd)
    annot = dict()
    comm_trans = [x.strip().split()[1].split('.')[0] for x in open(comm_trans_file)]
    with open(annot_out) as f:
        "{gene:{t1:exon, t2:exon, ...}, gene2:{...}}"
        for line in f:
            lst = line.strip().split('\t')
            bkp = lst[9]+':'+lst[11]
            annot.setdefault(bkp, dict())
            col9dict = dict(x.replace('"', '').split() for x in lst[8].rstrip(';').split(';'))
            if 'gene_name' in col9dict:
                gene = col9dict['gene_name']
                annot[bkp].setdefault(gene, dict())

            if 'transcript_id' in col9dict:
                transcript = col9dict['transcript_id']
                annot[bkp][gene].setdefault(transcript, 'intronic')

            if 'exon_number' in col9dict:
                exon = 'exon' + col9dict['exon_number']
                annot[bkp][gene][transcript] = exon
    print(annot)
    bkps = []
    for line in open(bed):
        lst = line.strip().split()
        bkp = lst[0] + ':' + lst[2]
        bkps.append(bkp)

    result = []
    for bkp in bkps:
        if bkp in annot:
            gte = []
            for gene in annot[bkp]:
                hit = [x for x in annot[bkp][gene].keys() if x.split('.')[0] in comm_trans]
                if not hit:
                    info = '|'.join(gene+':'+t+':'+e for t, e in annot[bkp][gene].items())
                else:
                    info = gene+':'+hit[0]+':'+annot[bkp][gene][hit[0]]
                gte.append(info)
            result.append(gte)
        else:
            print(f'{bkp} is not annotated !')
    return result


def process_avenio_result(result, sample_info, snvindel_target, fusion_target, cnv_target, comm_trans=None, gtf=None):
    """
    加工一下罗氏的结果表：
    1. 加入罗氏的样本id
    2. 提取常用转录本集对应的cphgvs和phgvs和exon_number
    3. phgvs 转换为单字母表
    4. 断点转录本注释，使用bedtools从gtf中提取转录本和外显子号信息，另外可以通过http://www.genome.ucsc.edu/cgi-bin/hgGateway验证结果。
    :param result:avenio的提取结果表
    :param sample_info:两列，无header，第一列是result中包含的样本id，第二列是新提供的新样本id，最后的结果将使用新的id
    :param snvindel_target: 两列，有header，第一列是新样本id,二列是genomic position，后续将合并这两列作为id进行结果提取，获得target结果。
    :param fusion_target: 两列，有header，第一列是新样本id,二列是genomic position，后续将合并这两列作为id进行结果提取，获得target结果。
    :param cnv_target: 两列，有header，第一列是新样本id,二列是genomic position，后续将合并这两列作为id进行结果提取，获得target结果。
    :param comm_trans: Ensembl常用转录本信息，第一列为基因名，第二列为常用转录本id: /nfs2/database/CommonRefSeqTranscripts/canonical_ensgene_hg19_v90.txt
    :param gtf: gtf: /nfs2/database/gencode_v29/gencode.v29.primary_assembly.annotation.gtf
    :return: 输出文件有
        raw1.result.xls，对result处理的初步结果，处理包括1，2，3过程
        raw1.fusion.xls，依据上表提取出fusion的结果后进行断点注释
        raw2.result.xlsx，从raw1.result.xls或raw1.fusion.xls提取特定列（详见代码）后的结果
        target.result.xlsx 从raw2.result.xlsx中根据*target参数提取结果
    """
    table = pd.read_csv(result, header=0, sep='\t')
    # 规范化AF
    formated_af = []
    for af in table['Allele Fraction']:
        if type(af) == float:
            formated_af.append(f'{af:.2%}')
        else:
            if af.endswith('%'):
                af = float(af[:-1]) / 100
            else:
                try:
                    af = float(af)
                except:
                    print(f'warn: failed to convert {af} to float')
                    af = 0.0
            formated_af.append(f'{af:.2%}')
    table['Allele Fraction'] = formated_af
    #
    sample_info_dict = dict([x.strip().split('\t')[:2] for x in open(sample_info)])
    target_samples = [x in sample_info_dict for x in table['Sample ID']]
    table = table.iloc[target_samples]
    table['Sample'] = [sample_info_dict[x] for x in table['Sample ID']]
    # table.set_index('NewID', inplace=True)
    #
    comm_trans_set = {x.strip().split()[1].split('.')[0] for x in open(comm_trans)}
    cols = ['Transcript', 'Coding Change', 'Amino Acid Change', 'Exon Number']
    comm_t_id = []
    for transcripts in table['Transcript']:
        matched = False
        if type(transcripts) == str:
            for ind, trans in enumerate(transcripts.split(";")):
                if trans.split('.')[0] in comm_trans_set:
                    comm_t_id.append(ind)
                    matched = True
                    break
        if not matched:
            comm_t_id.append(0)
    get = lambda x, ind: x.split(';')[ind] if type(x)==str else x
    table['Transcript'] = list(map(get, table['Transcript'], comm_t_id))
    table['Coding Change'] = list(map(get, table['Coding Change'], comm_t_id))
    table['Amino Acid Change'] = list(map(get, table['Amino Acid Change'], comm_t_id))
    # aa 3letter -> 1letter
    table['Amino Acid Change'] = list(map(aa_three2one, table['Amino Acid Change']))
    table['Exon Number'] = list(map(get, table['Exon Number'], comm_t_id))
    table['Total Exons'] = list(map(get, table['Total Exons'], comm_t_id))
    table['Variant Description'] = list(map(get, table['Variant Description'], comm_t_id))
    table.to_csv('raw1.result.xls', sep='\t')
    # snv indel
    target_cols = [
        'Sample',
        'Input DNA Mass (ng)',
        'Isolated DNA Mass (ng)',
        'Mutation Class',
        'Gene',
        'Transcript',
        'Exon Number',
        'Coding Change',
        'Amino Acid Change',
        'Variant Description',
        'Allele Fraction',
        'No. Mutant Molecules per mL',
        'Genomic Position',
        'Variant Depth',
        'Unique Depth',
        'dbSNP ID',
        'COSMIC ID'
    ]

    snv_indel = table.loc[(table['Mutation Class']=='SNV') | (table['Mutation Class'] == 'INDEL'), target_cols]

    # cnv
    target_cols = [
        'Sample',
        'Input DNA Mass (ng)',
        'Isolated DNA Mass (ng)',
        'Mutation Class',
        'Gene',
        'Genomic Position',
        'CNV Score'
    ]
    cnv = table.loc[table['Mutation Class']=='CNV', target_cols]

    # 注释融合断点
    fusion = table.loc[table['Mutation Class'] == 'FUSION']
    bkp1 = fusion['Fusion Breakpoint 1'].str.split(':', expand=True)
    bkp1.columns = ['chr', 'end']
    bkp1['start'] = bkp1['end'].astype(int) - 1
    bkp1[['chr', 'start', 'end']].to_csv('bkp1.bed', sep='\t', header=False, index=False)
    annot1 = annotation('bkp1.bed', comm_trans, gtf)
    fusion['Fusion Partner Transcript 1'] = [x[0] for x in annot1]

    bkp2 = fusion['Fusion Breakpoint 2'].str.split(':', expand=True)
    bkp2.columns = ['chr', 'end']
    bkp2['start'] = bkp2['end'].astype(int) - 1
    bkp2[['chr', 'start', 'end']].to_csv('bkp2.bed', sep='\t', header=False, index=False)
    annot2 = annotation('bkp2.bed', comm_trans, gtf)
    fusion['Fusion Partner Transcript 2'] = [x[0] for x in annot2]
    fusion.to_csv('raw1.fusion.xls', sep='\t')
    target_cols = [
        'Sample',
        'Input DNA Mass (ng)',
        'Isolated DNA Mass (ng)',
        'Mutation Class',
        'Gene',
        'Genomic Position',
        'Fusion Partner Gene 1',
        'Fusion Breakpoint 1',
        'Fusion Partner Transcript 1',
        'Fusion Partner Gene 2',
        'Fusion Breakpoint 2',
        'Fusion Partner Transcript 2',
        'Paired Reads Spanning Fusion Breakpoint'
    ]
    fusion = fusion[target_cols]
    with pd.ExcelWriter('raw2.result.xlsx') as writer:
        snv_indel.to_excel(writer, sheet_name='SNV&INDEL', index=False)
        fusion.to_excel(writer, sheet_name='Fusion', index=False)
        cnv.to_excel(writer, sheet_name='CNV', index=False)

    # 提取目标数据
    snvindel_target = pd.read_csv(snvindel_target, sep='\t', header=0)
    snvindel_targets = snvindel_target.iloc[:,0] + snvindel_target.iloc[:,1]
    fusion_target = pd.read_csv(fusion_target, sep='\t', header=0)
    fusion_targets = fusion_target.iloc[:, 0] + fusion_target.iloc[:, 1]
    cnv_target = pd.read_csv(cnv_target, sep='\t', header=0)
    cnv_targets = cnv_target.iloc[:, 0] + cnv_target.iloc[:, 1]
    snv_indel.index = snv_indel['Sample'] + snv_indel['Genomic Position']
    print('DupIndex:', snv_indel.loc[snv_indel.index.duplicated(keep=False)])
    fusion.index = fusion['Sample'] + fusion['Genomic Position']
    cnv.index = cnv['Sample'] + cnv['Genomic Position']
    # print(set(snvindel_targets) - set(snv_indel.index))
    # print(set(snvindel_targets) & set(snv_indel.index))
    # print(snv_indel.head())
    with pd.ExcelWriter('target.result.xlsx') as writer:
        in_dup_ind = set(snvindel_targets) & set(snv_indel.index[snv_indel.index.duplicated()])
        if in_dup_ind:
            print(f"Caution: dup index of {in_dup_ind}")
        snv_indel.loc[snvindel_targets].to_excel(writer, sheet_name='SNV&INDEL', index=False)
        fusion.loc[fusion_targets].to_excel(writer, sheet_name='Fusion', index=False)
        cnv.loc[cnv_targets].to_excel(writer, sheet_name='CNV', index=False)


def extract_hot(vcfs:list, avenio_table, hotspot, chain=None, ref=None, out='in_hotspot.xlsx'):
    """
    由于现有hotspot是hg19版的，而avenio结果是hg38版的，为确定avenio结果表中的突变是否属于hotspot，故要用crossMap进行坐标转换vcf。
    根据chr,start,ref,alt确定是否为同义突变
    :param vcfs: avenio的vcf
    :param avenio_table: avenio结果表
    :param hotspot: vcf版
    :param chain: 用于坐标转换
    :param ref: 用于坐标转换
    :param out:
    :return:
    """
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
    xcmds.xcmds(locals(), include=['cross_map', 'extract_hot', 'match_avenio_hot', 'process_avenio_result', 'annotation'])
