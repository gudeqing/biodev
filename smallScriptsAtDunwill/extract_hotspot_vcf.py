from pysam import VariantFile, VariantHeader
import pandas as pd


def parse_hotspot(hotspot):
    with open(hotspot) as f:
        line_dict = dict()
        for line in f:
            lst = line.strip().split('\t')
            if len(lst) != 7:
                print('skip bad line: ', line)
                continue
            line_dict['chr'] = lst[0]
            line_dict['start'] = int(lst[1])
            line_dict['end'] = int(lst[2])
            line_dict['coding_change'] = lst[3]
            line_dict['aa_change'] = lst[4]
            line_dict['transcript'] = lst[5]
            line_dict['gene_name'] = lst[6]
            yield line_dict


def get_hotspot(hotspot):
    """
    构建使用transcript_id和coding_change作为键的字典。
    由于annovar给出的注释是和已有的hotspot不一样的:
    c.7230_7286del:p.2410_2429del -> hotspot中del后面会跟着删除的序列
    c.872_873insCAGC:p.G291fs -> 这和hotspot中的表示是一致的
    :param hotspot:
    :return:
    """
    result = dict()
    # MUTYH:NM_001350650:exon13:c.932A>C:p.Q311P
    for line in parse_hotspot(hotspot):
        transcript = line['transcript'].split('.')[0]
        coding_change = line['coding_change']
        if 'del' in coding_change:
            # 更改表示方式，使其和annovar的注释一致
            coding_change = line['coding_change'].split('del')[0] + 'del'
        key = ':'.join((
            # line['gene_name'],
            transcript,
            coding_change
        ))
        result[key] = line

        if 'dup' in coding_change:
            # 转换成ins, 增加一个备选方式, 以防查询时由于表示方式不一致而漏掉
            coding_change = coding_change.replace('dup', 'ins')
            key = ':'.join((
                # line['gene_name'],
                transcript,
                coding_change
            ))
            result[key] = line

    return result


def format_header():
    header_info = [
        '##fileformat=VCFv4.2',
        '##assembly=hg19',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##INFO=<ID=AAChange_refGene,Number=.,Type=String,Description="AAChange_refGene annotation">',
        # '##FORMAT=<ID=None,Number=R,Type=Integer,Description="None">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    ]
    header = VariantHeader()
    for line in header_info:
        header.add_line(line)
    return header


def simplify_annovar_vcf(vcf, out_vcf):
    """
    从vcf提取包含AAChange_refGene的记录
    :param vcf:
    :param out_vcf:
    :return:
    """
    with VariantFile(vcf) as fr:
        # remove header
        for key in fr.header.info.keys():
            if key != 'AAChange_refGene':
                fr.header.info.remove_header(key)
        # 仅保留INFO中AAChange_refGene信息
        with VariantFile(out_vcf, 'w', header=fr.header) as fw:
            for record in fr:
                if record.info['AAChange_refGene'] != '.':
                    for key in record.info.keys():
                        if key != 'AAChange_refGene':
                            record.info.pop(key)
                    fw.write(record)
                else:
                    pass


def extract_hotspot(vcf, hotspot, exclude_hotspot=None):
    """
    从annovar注释中提取热点突变, 必须要求有'AAChange_refGene'注释
    :param vcf:
    :param hotspot:
    :param exclude_hotspot:
    :return:
    """
    hots = get_hotspot(hotspot)
    if exclude_hotspot:
        exclude_hots = get_hotspot(exclude_hotspot)
    else:
        exclude_hots = dict()
    # print('There are {} SNP/Ins/Del/Dup in our hotspot database'.format(len(hots)))

    with VariantFile(vcf) as fr:
        # open(vcf+'.hotspot.txt', 'w') as fw
        sample = list(fr.header.samples)[1]
        result = {sample:dict()}
        for record in fr:
            if not 'AAChange_refGene' in record.info:
                continue
            for each in record.info['AAChange_refGene']:
                # NOTCH2:NM_024408:exon34:c.6795T>G:p.N2265K
                if each == '.':
                    break
                tmp = each.split(':')
                if len(tmp) <= 3:
                    pass
                else:
                    sample = list(record.samples)[1]
                    key = ':'.join((tmp[1], tmp[3]))
                    if key in exclude_hots:
                        print(each, f'of {sample} is in low coverage region, and will be excluded')
                    if key in hots and key not in exclude_hots:
                        af = round(record.samples[sample]['AF'][0], 4)
                        print(f'{sample}'+':'+each+':'+f'AF={af}')
                        # fw.write(f'{sample}\t{each}\t{af}\n')
                        result[sample][each] = af
    return result


def batch_extract_hotspot(vcfs:list, hotspot, sample_info=None, exclude_hotspot=None, cmp_vcfs:list=None):
    """
    输入的vcf必须是经过annovar注释的，必须要求有'AAChange_refGene'注释，
    这里假设一个突变的唯一性可以由transcript_id和coding_change共同决定。
    :param vcfs: vcf 路径, 可以提供多个
    :param hotspot: 热点突变文件路径
    :param sample_info: 样本信息文件，第一列必须是样本id，第一行是header
    :param exclude_hotspot: 要排除的热点突变路径，默认不提供
    :param cmp_vcfs: vcf路径，与vcfs提供的数量应该一致，该参数是为了比较两种分析的结果是否一致而设计的。默认不提供
    :return: 生成all.detected.hotspot.xls文件，'sample\tmutation\tAF1\tAF2\tconsistency\n'
    """
    result = dict()
    for vcf in vcfs:
        result.update(extract_hotspot(vcf, hotspot, exclude_hotspot))

    out_file = 'all.detected.hotspot.xls'
    with open(out_file, 'w') as fw:
        fw.write('sample\tmutation\tAF\n')
        for sample, var_dict in result.items():
            if not var_dict:
                print(f'WARN: No hotspot mutation detected in {sample}')
            for mutation, af in var_dict.items():
                fw.write(f'{sample}\t{mutation}\t{af}\n')

    if cmp_vcfs:
        if len(cmp_vcfs) != len(vcfs):
            raise Exception('两组vcf数量不相等, 如何比较?')
        cmp_result = dict()
        for vcf in cmp_vcfs:
            cmp_result.update(extract_hotspot(vcf, hotspot, exclude_hotspot))
    else:
        return

    with open(out_file, 'w') as fw:
        fw.write('sample\tmutation\tAF1\tAF2\tconsistent\n')
        for sample, var_dict in result.items():
            mutations = list(var_dict.keys())
            if sample in cmp_result:
                mutations += list(cmp_result[sample].keys())
                var_dict2 = cmp_result[sample]
            else:
                raise Exception(f'{sample} 只在一组vcf中出现')
            if not var_dict and (not var_dict2):
                print(f'WARN: No hotspot mutation detected in {sample}')
            for mutation in mutations:
                if mutation in var_dict:
                    af = var_dict[mutation]
                else:
                    af = 0
                if mutation in var_dict2:
                    af2 = var_dict2[mutation]
                else:
                    af2 = 0
                if mutation in var_dict and (mutation in var_dict2):
                    consistent = 'yes'
                else:
                    consistent = 'no'
                fw.write(f'{sample}\t{mutation}\t{af}\t{af2}\t{consistent}\n')
    if sample_info:
        sample_info_df = pd.read_csv(sample_info, header=0, index_col=0, sep=None, engine='python')
        mutation_df = pd.read_csv(out_file, header=0, index_col=0, sep=None, engine='python')
        final_df = mutation_df.join(sample_info_df)
        final_df.to_excel(out_file[:-3]+'xlsx')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['batch_extract_hotspot'])


