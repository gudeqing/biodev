from pysam import VariantFile
import os
"""
这个脚本提供对vcf进行预处理，包括AF信息提取到INFO列，同时使用bcftools 进行 norm
"""


def reformat_vcf(vcf_file, out, reference, tumor_sample=None):
    """
    [split and left-normalization] -> [re-calculate AF] -> [bgzip and tabix]
    :param vcf_file:
    :param out:
    :param reference:
    :param tumor_sample: tumor sample name or position, if the second sample is tumor, then it is 1 else 0
    :return:
    """
    os.system(f'bcftools norm -f {reference} -m -both {vcf_file} -o tmp_.vcf')
    with VariantFile('tmp_.vcf') as fr:
        header = fr.header
        header.info.add('TDP', number=1, type='Integer', description='Tumor sample depth')
        header.info.add('NDP', number=1, type='Integer', description='Normal sample depth')
        header.info.add('TAF', number=1, type='Float', description='Tumor sample AF')
        header.info.add('NAF', number=1, type='Float', description='Normal sample AF')
        samples = list(header.samples)

        if tumor_sample is not None:
            if tumor_sample not in [0, 1, '0', '1']:
                if tumor_sample in samples:
                    tumor_idx = samples.index(tumor_sample)
                    normal_idx = 1 - tumor_idx
                else:
                    raise Exception(f'{tumor_sample} is not in samples {samples} recorded in vcf')
            else:
                tumor_idx = int(tumor_sample)
                normal_idx = 1 - tumor_idx
        else:
            tumor_idx = guess_tumor_idx(vcf_file)
            normal_idx = 1 - tumor_idx

        with VariantFile(out, 'w', header=header) as fw:
            for record in fr:
                record.info['TDP'] = record.samples[tumor_idx]['DP']
                record.info['NDP'] = record.samples[normal_idx]['DP']
                # re-calculate AF since caller like sentieon may report AF that is not consistent with AD info
                record.info['TAF'] = round(record.samples[tumor_idx]['AD'][1]/record.samples[tumor_idx]['DP'], 3)
                if record.samples[normal_idx]['DP'] != 0:
                    record.info['NAF'] = round(record.samples[normal_idx]['AD'][1]/record.samples[normal_idx]['DP'], 3)
                else:
                    record.info['NAF'] = 0
                fw.write(record)

    os.remove('tmp_.vcf')
    os.system(f'bgzip {out}')
    os.system(f'tabix {out}.gz')


def guess_tumor_idx(vcf_file):
    tumor_is_first = 0
    tumor_is_second = 0

    with VariantFile(vcf_file) as fr:
        samples = list(fr.header.samples)
        formats = list(fr.header.formats)
        if 'AF' not in formats:
            raise Exception('No AF in format info to detect tumor sample')
        for record in fr:
            if record.samples[0]['AF'][0] > record.samples[1]['AF'][0]:
                tumor_is_first += 1
            else:
                tumor_is_second += 1
    tumor_idx = tumor_is_second >= tumor_is_first
    print(f'we guess tumor sample is {samples[tumor_idx]} ')
    return tumor_idx


def prepare_for_purity_imputation(vcf, ploidy=2):
    """
    https://github.com/KhiabanianLab/All-FIT/blob/master/test/input/sampleFile1.xls
    :return:
    """
    tumor_idx = guess_tumor_idx(vcf)
    info = []
    with VariantFile(vcf) as fr:
        samples = list(fr.header.samples)
        tumor_name = samples[tumor_idx]
        for record in fr:
            var_id = f'{record.contig}:{record.pos}:{record.ref}:{record.alts[0]}'
            depth = record.samples[tumor_idx]['DP']
            freq = record.samples[tumor_idx]['AF'][0]*100
            info.append([var_id, freq, depth, ploidy])

    with open(f'{tumor_name}.var_freq_dp_ploidy.txt', 'w') as f:
        f.write('ID\tAllele_Freq\tDepth\tPloidy\n')
        for var_id, freq, depth, ploidy in info:
            f.write(f'{var_id}\t{freq}\t{depth}\t{ploidy}\n')


def simplify_report(out_prefix, tsv_files='*_tumor/*.pcgr_acmg.grch37.pass.tsv.gz', target_genes=None):
    import glob
    import pandas as pd

    files = glob.glob(tsv_files)

    tables = []
    target_genes = target_genes
    target_genes = ['ENSG00000100462', 'ENSG00000124523', 'ENSG00000142082', 'ENSG00000198890', 'ENSG00000133703']

    for i in files:
        a = pd.read_table(i, comment='#')
        a['sample'] = i.split('/')[0]
        if target_genes:
            b = a[[x in target_genes for x in a['Gene']]]
            tables.append(b)
        else:
            tables.append(a)

    bt = pd.concat(tables)

    tcols = ['sample', 'CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'STRAND', 'ENSEMBL_GENE_ID', 'SYMBOL',
             'ENSEMBL_TRANSCRIPT_ID', 'REFSEQ_MRNA', 'CANONICAL', 'UNIPROT_ACC', 'Consequence', 'VARIANT_CLASS', 'TDP',
             'TAF', 'HGVSc', 'HGVSp', 'HGVSp_short', 'Existing_variation', 'TCGA_DRIVER', 'TCGA_FREQUENCY',
             'EAS_AF_GNOMAD', 'CLINVAR_CLASSIFICATION', 'MUTATION_HOTSPOT', 'ONCOGENE', 'PUTATIVE_DRIVER_MUTATION',
             'STR']

    with pd.ExcelWriter(f'{out_prefix}.xlsx') as writer:
        bt[tcols].to_excel(writer, sheet_name='summary', index=False)
        header_desc = {
            "sample": "样本名称",
            "CHROM": "染色体",
            "POS": "突变的起始坐标",
            "REF": "参考序列",
            "ALT": "突变后的序列",
            "FILTER": "PASS表示通过TNFilter过滤条件",
            "STR": "突变是否在短串联重复区域",
            "TDP": "肿瘤样本测序深度",
            "TAF": "突变频率",
            "CLINVAR_CLASSIFICATION": "clinvar数据库的突变分级",
            "TCGA_FREQUENCY": "突变在TCGA数据库中不同癌症的统计信息，格式：癌种|百分比|有突变的case数量|总cases数量",
            "PUTATIVE_DRIVER_MUTATION": "依据TCGA数据库中9423个肿瘤外显子数据分析预测的驱动基因，格式：symbol:hgvsp:ensembl_transcript_ids:discovery_approaches",
            "MUTATION_HOTSPOT": "已知肿瘤热点突变，依据cancerhotspots.org_v2数据库，格式：Gene|Codon|Q-value",
            "CANONICAL": "是否被VEP标记为经典转录本",
            "Consequence": "Impact modifier for the consequence type",
            "EAS_AF_GNOMAD": "东亚人群突变频率",
            "ENSEMBL_GENE_ID": "Ensembl基因ID",
            "ENSEMBL_TRANSCRIPT_ID": "Ensembl转录本ID",
            "Existing_variation": "已报道的突变ID，包括NCBI和COSMIC数据库的ID",
            "HGVSc": "基于核酸序列的HGVS格式描述",
            "HGVSp": "基于氨基酸序列的HGVS格式描述",
            "HGVSp_short": "基于氨基酸序列的HGVS格式描述，用单字母表示氨基酸",
            "ONCOGENE": "是否被标记为致癌基因（from Network of Cancer Genes (NCG) & the CancerMine text-mining resource）",
            "REFSEQ_MRNA": "RefSeq转录本id",
            "STRAND": "表示基因在DNA的正链或者负链",
            "SYMBOL": "基因名称",
            "TCGA_DRIVER": "是否为癌症驱动基因（from Pan-Cancer analysis）",
            "UNIPROT_ACC": "UniprotKB accession identifiers",
            "VARIANT_CLASS": "突变类型（VEP格式）"
        }

        desc = pd.DataFrame([header_desc])
        desc[[x for x in tcols if x in desc.columns]].T.to_excel(writer, sheet_name='header_desc', header=False)

    bt[tcols].to_csv(f'{out_prefix}.txt', sep='\t', index=False)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())

