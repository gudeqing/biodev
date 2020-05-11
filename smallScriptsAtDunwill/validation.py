import os
import re
import logging
import pysam
import statistics
import numpy as np
import pandas as pd
from pysam import VariantFile, VariantHeader
from statsmodels.stats.proportion import proportion_confint as pconf

"""
性能验证思路：
（0）如何确定两个突变一致？ 我目前的建议是使用'transcript:chgvs'作为唯一var_id，var_id相同则认为突变相同。
（1）获得样本信息，获得样本分组信息，比如重复性分组，LOD分组，input分组，这些组的样本共享已知突变
（2）获得已知突变信息，要精确到每个样本或分组，最终可以制成{已知突变字典}：{sample： {var_id: [mutation, AF]} }； 另外要确定背景突变位点数，用于计算真阴性等统计。
（3）提取每个样本的突变信息。通常从vcf出发，如对vcf用相同软件做注释，然后根据注释结果提取每个样本的突变结果，形成一个汇总表。
    我的做法或建议如下：
    for each sample：
        a. 如果vcf来自不同软件的结果，可以考虑先用‘bcftools norm’对vcf进行统一的标准化
        b. 对标准化的vcf进行annovar注释
        c. 对注释后的vcf进行简化，如提取出经典转录本的注释
        d. 制作{检测突变字典}， 形如 {sample： {var_id: [mutation, AF]} }
        e. 可选，用hotspot或其他信息进行突变筛选，得到最终的{检测突变字典}
    else:
        合并所有样本的突变字典-->{检测突变字典}，使用pandas生成汇总表，得到的文件中，你会发现每一行代表一个样本的一个突变。
（4）统计：
    a. 按照样本逐一循环，比较每个样本的{已知突变字典}和{检测突变字典}，得到每个样本的tp,fp等统计信息
    b. 汇总所有样本的tp,fp等统计信息, 得到总体tp,fp等统计信息。 最后的汇总表中每一行代表一个样本，最后一行是汇总结果。
    c. lod和input和replicate设计的统计形式一样，可统称为重复性分析，即针对一组样本进行统计：
        for samples in each_sample_group:
            for mutation in known_mutation_dict of current group:
                len(<? mutation in sample_detected_mutation_dict> for detected_mutation_dict of AllSamples)/len(samples)
                --> 10/11, 11个样本中有10个样本检测到该突变
            else:
                得到当前组的已知突变检测统计结果
"""

# fasta_file = '/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta'
# fasta_file = '/nfs2/database/gencode_v29/GRCh38.primary_assembly.genome.fa'
# "/nfs2/software/STAR-Fusion-v1.6.0/refgenome/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa"

def parse_coding_change(cc):
    m_type = None
    tp1 = re.compile(r'[cg].[*]?[0-9+-]+([ATCG])>([ATCG])')
    tp2 = re.compile(r'[cg].[0-9]+del([ATCG]*)')
    # tp3 = re.compile(r'[cg].([0-9]+)_([0-9]+)([ATCG]+)')
    tp3 = re.compile(r'[cg].([0-9]+)_([0-9]+)[ATCG]*[>]?([ATCG]+)')
    tp4 = re.compile(r'[cg].([0-9]+)_([0-9]+)ins([ATCG]+)')
    tp5 = re.compile(r'[cg].([0-9]+)[_]?([0-9]*)del([ATCG]*)ins([ATCG]+)')
    tp6 = re.compile(r'[cg].([0-9]+)_([0-9]+)del([ATCG]*[0-9]*)')
    tp7 = re.compile(r'[cg].([0-9]+)_([0-9]+)dup([ATCG]+)')
    tp8 = re.compile(r'[cg].[0-9]+dup([ATCG]?)')
    tp9 = re.compile(r'[cg].[0-9]+ins([ATCG])')
    matched_ind = None
    matched_obj = None
    for ind, pattern in enumerate([tp1, tp2, tp3, tp4, tp5, tp6, tp7, tp8, tp9], start=1):
        matched_obj = pattern.fullmatch(cc)
        if matched_obj:
            matched_ind = ind
            break

    if not matched_ind:
        # raise Exception(f"{cc} cannot be parsed!")
        return None, f"{cc} cannot be parsed!"
    return matched_ind, matched_obj


def reverse_complement(seq):
    """
    :param seq: 输入序列
    :return: 返回reversed的序列
    """
    seq = seq.upper()
    complement = dict(zip(list("ATCG"), list("TAGC")))
    return ''.join(complement[base] if base in complement else base for base in seq[::-1])


def translate_to_vcf_5cols(chr_name, start, coding_change, genome, strand='+', mut_id=None):
    """
    :param chr_name: chromosome name, such as chr7
    :param start: 1-based start coordinate, 这个坐标是突变发生时的第一个碱基，是满足vcf格式要求的pos
    :param coding_change: coding change, such as 'c.2156_2157AT', 因为这个信息是基于转录本序列的，为此我们需要知道基因所在链信息
    :param genome: genome fasta file or pysam object of genome
    :param strand: '-' or '+' indicate which strand transcribes the gene
    """
    start = int(start)
    if type(genome) == str:
        genome = pysam.FastaFile(genome)
    ind, matched = parse_coding_change(coding_change)
    if mut_id is None:
        mut_id = f'{chr_name}:{start}:{coding_change}'
    if not ind:
        return f'{mut_id} cannot be parsed!'

    if ind == 1:
        # snv
        ref = genome.fetch(chr_name, start-1, start).upper()
        base, alt = matched.groups()
        if strand == '-' and (not coding_change.startswith('g.')):
            base = reverse_complement(base)
            alt = reverse_complement(alt)
        if ref != base:
            return f"{mut_id}中的碱基和参考基因组ref不一致:{base} vs {ref}"
        return chr_name, start, mut_id, ref, alt

    elif ind == 2:
        # single del
        alt = matched.groups()[0]
        if strand == '-' and (not coding_change.startswith('g.')):
            alt = reverse_complement(alt)
        ref = genome.fetch(chr_name, start-1, start+1).upper()
        if alt and alt != ref[1]:
            return f"{mut_id}中的碱基alt和参考基因组ref不一致:{alt} vs {ref[1]}"
        return chr_name, start, mut_id, ref, ref[0]

    elif ind == 3:
        # substitution, 本质上是一种等长替换
        s, t, alt = matched.groups()
        if strand == '-' and (not coding_change.startswith('g.')):
            alt = reverse_complement(alt)
        if int(t) - int(s) + 1 != len(alt):
            raise Exception(f'从{mut_id}获得的替换长度不一致')
        ref = genome.fetch(chr_name, start-1, start-1+len(alt)).upper()
        return chr_name, start, mut_id, ref, alt

    elif ind == 4:
        # insertion
        s, t, alt = matched.groups()
        if strand == '-' and (not coding_change.startswith('g.')):
            alt = reverse_complement(alt)
        if int(t) - int(s) != 1:
            raise Exception(f'{mut_id}中cc坐标似乎不对，两个坐标相减应该=1')
        ref = genome.fetch(chr_name, start-1, start).upper()
        return chr_name, start, mut_id, ref, ref+alt

    elif ind == 5:
        # del and insertion, 本质上是一种替换
        s, t, deletion, alt = matched.groups()
        if not t:
            t = s
        if strand == '-' and (not coding_change.startswith('g.')):
            alt = reverse_complement(alt)
            deletion = reverse_complement(deletion)
        del_len = int(t) - int(s) + 1
        if deletion and len(deletion) != del_len:
            raise Exception(f"{mut_id}的deletion长度不一致")
        ref = genome.fetch(chr_name, start-1, start-1+del_len).upper()
        if deletion and ref != deletion:
            return f"{mut_id}中的deletion碱基和参考基因组ref不一致:{deletion} vs {ref}"
        return chr_name, start, mut_id, ref, alt 

    elif ind == 6:
        # deletion(>=2)
        s, t, alt = matched.groups()
        if alt.isnumeric():
            alt = ''
        if strand == '-' and (not coding_change.startswith('g.')):
            alt = reverse_complement(alt)
        del_len = int(t) - int(s) + 1
        if alt and del_len != len(alt):
            raise Exception(f"{mut_id}中deletion长度不一致")
        ref = genome.fetch(chr_name, start-1, start+del_len).upper()
        if alt and ref[1:] != alt:
            # print(f'{mut_id}获得deletion和参考基因组的不一致')
            # 如果给的不是起始坐标，而是最后的坐标，可以通过下面矫正回来
            # 例子如：chr10   43609939    p.L629_D631>H           RET NM_020975   c.1886_1891delTGTGCG
            start = start - del_len
            ref = genome.fetch(chr_name, start-1, start+del_len).upper()
            if ref[1:] != alt:
                return f"{mut_id}中的deletion碱基和参考基因组ref仍然不一致:{alt} vs {ref[1:]}"

        return chr_name, start, mut_id, ref, ref[0]

    elif ind == 7:
        # dup(>2)
        s, t, alt = matched.groups()
        if strand == '-' and (not coding_change.startswith('g.')):
            alt = reverse_complement(alt)
        dup_len = int(t) - int(s) + 1
        if len(alt) != dup_len:
            raise Exception(f'{mut_id}中dup长度信息不一致')
        ref = genome.fetch(chr_name, start-1, start-1+len(alt)).upper()
        if ref != alt:
            # 当成终止坐标试试
            start = start - dup_len + 1
            ref = genome.fetch(chr_name, start-1, start-1+len(alt)).upper()
            if ref != alt:
                return f"{mut_id}中的dup碱基和参考基因组的不一致:{alt} vs {ref}"
        return chr_name, start, mut_id, ref, ref+alt

    elif ind == 8:
        # single dup
        alt = matched.groups()[0]
        if alt and strand == '-' and (not coding_change.startswith('g.')):
            alt = reverse_complement(alt)
        ref = genome.fetch(chr_name, start-1, start).upper()
        if alt and alt != ref:
            return f"{mut_id}中的碱基和参考基因组不一致:{alt} vs {ref}"
        return chr_name, start, mut_id, ref, ref+ref

    elif ind == 9:
        # single insertion
        alt = matched.groups()[0]
        if strand == '-' and (not coding_change.startswith('g.')):
            alt = reverse_complement(alt)
        ref = genome.fetch(chr_name, start-1, start).upper()
        return chr_name, start, mut_id, ref, ref+alt
    else:
        raise Exception("unknown matching type error")


def set_logger(name='log.info', logger_id='x'):
    logger = logging.getLogger(logger_id)
    fh = logging.FileHandler(name, mode='w')
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    return logger


def vcf_header():
    header_info = [
        '##fileformat=VCFv4.2',
        '##assembly=hg19',
        "##contig=<ID=chr1,length=249250621>",
        "##contig=<ID=chr2,length=243199373>",
        "##contig=<ID=chr3,length=198022430>",
        "##contig=<ID=chr4,length=191154276>",
        "##contig=<ID=chr5,length=180915260>",
        "##contig=<ID=chr6,length=171115067>",
        "##contig=<ID=chr7,length=159138663>",
        "##contig=<ID=chrX,length=155270560>",
        "##contig=<ID=chr8,length=146364022>",
        "##contig=<ID=chr9,length=141213431>",
        "##contig=<ID=chr10,length=135534747>",
        "##contig=<ID=chr11,length=135006516>",
        "##contig=<ID=chr12,length=133851895>",
        "##contig=<ID=chr13,length=115169878>",
        "##contig=<ID=chr14,length=107349540>",
        "##contig=<ID=chr15,length=102531392>",
        "##contig=<ID=chr16,length=90354753>",
        "##contig=<ID=chr17,length=81195210>",
        "##contig=<ID=chr18,length=78077248>",
        "##contig=<ID=chr20,length=63025520>",
        "##contig=<ID=chr19,length=59128983>",
        "##contig=<ID=chr22,length=51304566>",
        "##contig=<ID=chr21,length=48129895>",
        "##contig=<ID=chrMT,length=16569>",
        "##contig=<ID=chrY,length=59373566>",
        "##contig=<ID=chrX,length=156040895>",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##INFO=<ID=Source,Number=.,Type=String,Description="source of mutation">',
        '##INFO=<ID=AAChange_refGene,Number=.,Type=String,Description="AAChange_refGene annotation">',
        # '##FORMAT=<ID=None,Number=R,Type=Integer,Description="None">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    ]
    return header_info


def trans(raw_mutation, out='result.vcf', gene_strand='/nfs2/database/gencode_v29/gene_strand.pair', genome='/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta'):
    """
    当初写这个脚本的本意是为了根据chr,start,chgvs的信息反推出ref，alt等信息，形成vcf.
    如此，当从数据库中下载的突变信息不够完整的时候就可以重新得到vcf进行注释. 当然去mutalyzer转换也是一种选择。
    :param raw_mutation: 文件路径,文件需要四列(chr_name\tstart_position\tgene_name\tcoding_change),第四列的信息可以是g.xxx|c.xxx
    :param gene_strand: 文件路径，两列(gene\tstrand),记录基因所在链的信息
    :param genome: 基因组fasta路径
    最后生成vcf
    """
    strand_dict = dict(x.strip().split()[:2] for x in open(gene_strand))
    logger = set_logger('trans.log')
    genome = pysam.FastaFile(genome)
    col5_list = []
    for raw_line in open(raw_mutation):
        line = dict(zip(['chr', 'start', 'gene_name', 'coding_change'], raw_line.strip().split()))
        # 预处理，染色体名称修正或排除，coding_change筛选或修正，gene_name更正
        if line['chr'].isnumeric():
            line['chr'] = 'chr'+line['chr']
        if line['chr'] in ['X', 'Y', 'MT']:
            line['chr'] = 'chr' + line['chr']

        if ',' in line['coding_change']:
            # 为了处理civic的数据而作的特殊处理
            tmplst = line['coding_change'].split(',')
            valid_cc = []
            for each in tmplst:
                tmp = each
                if ':' in each:
                    tmp = each.split(':')[1]
                if tmp.startswith('c.') or tmp.startswith('g.'):
                    valid_cc.append(tmp)
            if valid_cc:
                if len(valid_cc) > 0:
                    line['coding_change'] =  sorted(valid_cc)[0]
            else:
                raise Exception(f"{raw_line} -> has no valid coding_change expression!")
        else:
            if ':' in line['coding_change']:
                line['coding_change'] = line['coding_change'].split(':')[1]

        if line['coding_change'].startswith('c.1-1_'):
            # 由于要用到坐标计算长度信息，这里必须要修改一下，规范的表达不可以‘c.0'开的头
            line['coding_change'] = 'c.0' + line['coding_change'][5:]

        if 'MT' in line['chr']:
            continue

        gene = line['gene_name']
        if gene in strand_dict:
            strand = strand_dict[gene]
        elif gene.split('_')[0] in strand_dict:
            strand = strand_dict[gene.split('_')[0]]
        else:
            print(f'{gene} is not found in {gene_strand}, and we assume the gene is on + strand')
            strand = '+'
       
        # 开始转换
        start = int(float(line['start']))
        for st in [start, start-1, start+1]:
            col5 = translate_to_vcf_5cols(line['chr'], st, line['coding_change'], genome, strand=strand)
            if type(col5) == str:
                continue
            else:
                break

        if type(col5) == str:
            logger.info(col5)
        else:
            if col5 not in col5_list:
                col5_list.append(col5)
            else:
                print(f'{line} --> has a previous found variant!')
    # 写入vcf
    col5_list = sorted(col5_list, key=lambda x:(x[0], x[1]))
    vcf = set_logger(out, logger_id='result')
    for line in vcf_header():
        vcf.info(line)
    for each in col5_list:
        each = list(each) + ['.']*3
        vcf.info('\t'.join((str(x) for x in each)))


def cat_vcf(vcfs:list, out='merged.vcf'):
    """
    合并vcf并依据1-2-4-5列信息去重，且去掉3, 6，7，8列信息, 当初写这个脚本是为了合并trans的结果
    :param vcfs: vcf 路径，可以是多个，如果是一个，只仅仅是达到了去重的功能
    """
    vcf_info = set()
    for vcf in vcfs:
        with open(vcf) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                info = line.strip().split()[:5]
                # 清空第三列的id信息，保证利用1，2，4，5列信息去重
                # if len(vcfs) > 1:
                info[2] = '.'
                chr_name = info[0]
                if chr_name.isnumeric():
                    info[0] = 'chr'+chr_name
                if chr_name in ['X', 'Y', 'MT']:
                    info[0] = 'chr' + chr_name
                vcf_info.add(tuple(info))
    
    vcf = set_logger(out, logger_id='merged')
    print(f'find {len(vcf_info)} unique variants')

    # for line in vcf_header():
        # vcf.info(line)
    with open(vcfs[0]) as f:
        for line in f:
            if line.startswith('#'):
                vcf.info(line.strip())

    for each in sorted(vcf_info, key=lambda x:(x[0],x[1])):
        each = list(each) + ['.', 'PASS', '.']
        vcf.info('\t'.join((str(x) for x in each)))


def dedup_vcf(vcf, out='dedup.vcf'):
    """
    当用bcftools norm完vcf后，会发现起始有重复的突变，这时可以
    根据chr+start+ref+alt信息去重，仅保留第一个出现的记录
    :param vcf:
    :param out:
    :return:
    """
    with open(vcf) as f, open(out, 'w') as f2, open('old2new.pair', 'w') as f3:
        id_set = set()
        for line in f:
            if line.startswith('#'):
                f2.write(line)
            else:
                uniq_id = line.strip().split()[:5]
                old_id = uniq_id[2]
                uniq_id[2] = '.'
                uniq_id = ':'.join(uniq_id)
                f3.write(f'{old_id}\t{uniq_id}\n')
                if uniq_id in id_set:
                    print('discard duplicated mutation:', line)
                    continue
                else:
                    id_set.add(uniq_id)
                    f2.write(line)


def hotspot2msk(gene_phgvs, out='msk.list.xls'):
    """
    根据phgvs中ref氨基酸信息作为id进行合并，如p.E454K -> p.E454, p.E454_E457delinsG -> p.E454_E457
    :param gene_phgvs: 文件，两列，无header，第一列时gene，第二列是phgvs
    :param out:
    :return:
    """
    with open(gene_phgvs) as f:
        result = dict()
        for line in f:
            lst = line.strip().split()
            if len(lst) < 2:
                print('Skip line:', line)
            gene, phgvs = lst
            result.setdefault(gene, dict())
            try:
                if '_' not in phgvs:
                    first_aa = re.match('p.[A-Z][0-9]+', phgvs).group()
                else:
                    first_aa = re.match('p.[A-Z][0-9]+_[A-Z][0-9]+', phgvs).group()
                result[gene].setdefault(first_aa, list())
                result[gene][first_aa].append(phgvs)
            except Exception as e:
                print(line)
                print(e)

    with open(out, 'w') as f:
        for gene, first_aa in result.items():
            tmp_list = list(first_aa.keys())
            tmp_list = [(x, re.match('p.([A-Z])([0-9]+)', x).groups()) for x in tmp_list]
            # 按位置和氨基酸简写排序
            tmp_list = sorted(tmp_list, key=lambda x:(int(x[1][1]), x[1][0]))
            print(tmp_list)
            lst = [gene, ','.join(x[0] for x in tmp_list)]
            f.write('\t'.join(lst)+'\n')

    with open(out+'.source.xls', 'w') as f:
        for gene, first_aa in result.items():
            for faa, mutations in first_aa.items():
                lst =[gene, faa, ','.join(mutations)]
                f.write('\t'.join(lst)+'\n')


def parse_mycancergenome_var(mutation_file, out='mycancergeneome.var.txt'):
    """
    处理从mycancergenome获得的突变信息，转换成trans需要的四列文件
    """
    with open(mutation_file) as f, open(out, 'w') as fw:
        _ = f.readline()
        for line in f:
            gene, info = line.split('\t', 2)[:2]
            for each in info.split('__'):
                if len(each.split(':')) <= 1:
                    raise Exception(f'skip invalid: {each}')
                else:
                    chr_name, gc = each.split(':')
                start = re.match('g.([0-9]+)', gc).groups()[0]
                tmp = [chr_name, start, gene, gc]
                fw.write('\t'.join(tmp)+'\n')


def parse_hotspot(hotspot):
    """解析之前张慧给的hotspot文件，因为重新整理hotspot为vcf，不再使用"""
    with open(hotspot) as f:
        line_dict = dict()
        for line in f:
            lst = line.strip().split('\t')
            if len(lst) != 7:
                print('skip bad line: ', line)
                continue
            line_dict['chr'] = lst[0]
            # line_dict['start'] = int(lst[1])
            line_dict['start'] = lst[1]
            # line_dict['end'] = int(lst[2])
            line_dict['end'] = lst[2]
            line_dict['coding_change'] = lst[3]
            line_dict['aa_change'] = lst[4]
            line_dict['transcript'] = lst[5]
            line_dict['gene_name'] = lst[6]
            yield line_dict


def rm_dup_hotspot(hotspot, often_used, out='new_hotspot.bed'):
    """
    应该是再也用不到了，当初是用了对张慧提供的hotspot进行去重
    1. 尽量用坐标和核酸变化作为key，其他值作为value
    2. 如果key对应仅有一个value，则保留，不能使用的是否为常用转录本
    3. 如果key对应有多个values，如果使用的转录本均不是常用转录本，那么仅保留第一个value
    4. 如果key对应有多个values，仅保留包含有常用转录本记录的那一条
    """
    often_dict = dict(x.strip().split('\t') for x in open(often_used))
    transcripts = set(often_dict.values())
    transcripts = {x.split('.')[0] for x in transcripts}
    raw_hots = dict()
    for line_dict in parse_hotspot(hotspot):
        line = [str(x) for x in line_dict.values()]
        key = line[:3]
        coding_change = line[3]
        tmp = coding_change.split('|')
        tmp_len = [len(x) for x in tmp]
        coding_change = sorted(zip(tmp, tmp_len), key=lambda x:x[1])[-1][0]
        line[3] = coding_change
        nt_change = re.findall('[ATCG]+.*', coding_change)
        if not nt_change:
            # 发现存在几条没有核酸序列的记录，如c.1648_1674del27，但这个是可以的
            pass
        else:
            nt_change = nt_change[0]
            key.append(nt_change)
        key = tuple(key)
        if key not in raw_hots:
            raw_hots[key] = [line[3:]]
        else:
            raw_hots[key].append(line[3:])
    print('Raw mutation number', len(raw_hots))
    out = open(out, 'w')
    for key in raw_hots:
        values = raw_hots[key]
        if len(values) == 1:
            mutation = list(key)[:3]+values[0]
            if values[0][2].split('.')[0] not in transcripts:
                if values[0][3] in often_dict:
                    print(f'没有用常用转录本{often_dict[values[0][3]]}, 但仍然保留，因为只有一记录', mutation)
                else:
                    print('我们没有记录这个基因的常用转录本', values[0][3], '但仍然保留，因为只有一条记录',mutation)
            out.write('\t'.join(mutation)+'\n')
        else:
            ts = {x[2].split('.')[0] for x in values}
            
            if not (ts & transcripts):
                for each in values:
                    mutation = list(key)[:3]+each
                    out.write('\t'.join(mutation)+'\n')
                    print('下面的突变均没有记录的常用转录本,仅保留第一个', key, values)
                    break
            else:
                used_ts = set()
                for each in values:
                    mutation = list(key)[:3]+each
                    if each[2].split('.')[0] not in transcripts:
                        print(each[2].split('.')[0], '不是常用转录本,但我们找到包含常用转录本的其他记录,因此跳过该突变:', mutation)
                    else:
                        if each[2].split('.')[0] in used_ts:
                            continue
                        used_ts.add(each[2].split('.')[0])
                        mutation = list(key)[:3]+each
                        print('保留包含常用转录本的突变:', mutation)
                        out.write('\t'.join(mutation)+'\n')
                #print(used_ts)
    out.close()


def get_hotspot_old(hotspot):
    """
    应该也是不再需要了，当初编写该脚本是为了处理张慧这边提供的hotspot信息
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
        tmp = coding_change.split('|')
        tmp_len = [len(x) for x in tmp]
        coding_change = sorted(zip(tmp, tmp_len), key=lambda x:x[1])[-1][0]
        if ('del' in coding_change) and ('ins' not in coding_change):
            # 更改表示方式，使其和annovar的注释一致
            if '_' in coding_change:
                coding_change = coding_change.split('del')[0] + 'del'
            else:
                # coding_change只包含一个核酸del时，要保留核酸字母
                pass
        elif ('ins' in coding_change) and ('del' not in coding_change):
            # 转换成ins, 增加一个备选方式, 以防查询时由于表示方式不一致而漏掉
            coding_change = coding_change.replace('dup', 'ins')
        elif ('ins' in coding_change) and ('del' in coding_change):
            #print('Skip special mutation with both ins and del ! 因为annovar的注释好像不存在这种情况！')
            continue

        key = ':'.join((
            # line['gene_name'],
            transcript,
            coding_change
        ))

        result[key] = line

    return result


def get_hotspot(vcf, id_mode='transcript:chgvs'):
    """
    :param vcf: annovar注释过的且info列只保留AAChange_refGene信息且仅保留经典转录本注释的vcf文件。
    :param id_mode: 指定mutation的唯一id格式, 默认为'transcript:chgvs'.
        支持的字段有 ['chr', 'start', 'id', 'ref', 'alt', 'gene', 'transcript', 'exon', 'chgvs', 'phgvs']
    :return: 返回一个字典，以字典的形式存储突变, 方便查询
    """
    print('你指定的mutaion唯一id格式为:', id_mode)
    result = dict()
    with open(vcf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            lst = line.strip().split()
            hgvs =  lst[7].split('=')[1]
            info_dict = dict(zip(['gene', 'transcript', 'exon', 'chgvs', 'phgvs'], hgvs.split(':')))
            info_dict.update(dict(zip(['chr', 'start', 'id', 'ref', 'alt'], lst[:5])))
            key = ':'.join(info_dict[x] for x in id_mode.split(':'))
            if key in result:
                # raise Exception(f"在 {hotspot_vcf} 中, 根据当前vard_id={key}作为突变id时，突变不唯一!")
                print(f"在 {vcf} 中, 根据{key}作为突变id时，突变不唯一!")
            # result[key] = line
            result[key] = hgvs
    return result


def simplify_annovar_vcf(vcf, out_prefix, often_trans=None, filter=None):
    """
    从vcf提取包含AAChange_refGene的记录，如果提供经典转录本信息，则仅保留经典转录本的注释。
    另外如果还提供filter参数, 则依据filter里的每一行条件去匹配注释信息，默认仅保留匹配到的第一个条目。
    先filter，然后进行经典转录本筛选。
    :param vcf:
    :param out_prefix:
    :param filter: 文件路径，文件的每行包括匹配条件，到时结果仅显示匹配到的条目。
     示例如下，注意包含染色体色信息和基因信息可以提高匹配速度:
     (1) g=ABL1;t=NM_005157;e=exon4;c=c.763G>A;p=p.E255K
     (2) 支持的key有: [r, s, g, t, e, c, p] , 他们的含义分别是 [chr,start,gene,transcript,exon,chgvs,phgvs]
    :return:
    """
    if often_trans:
        often_dict = dict(x.strip().split('\t')[:2] for x in open(often_trans) if len(x.strip().split())>1)
    else:
        often_dict = dict()

    filter_lst = []
    filter_genes = set()
    filter_chrs = set()
    filter_starts = set()
    if filter: 
        with open(filter) as f:
            for line in f:
                line_dict = dict(x.strip().split('=') for x in line.split(';'))
                invalid_key = set(line_dict.keys()) - set('rsgtecp')
                if invalid_key:
                    raise(f'Invalid keys: {invalid_key}')
                gene = line_dict.get('g')
                if gene:
                    filter_genes.add(gene)
                chr_name = line_dict.get('r')
                if chr_name:
                    filter_chrs.add(chr_name)
                start = line_dict.get('s')
                if start:
                    filter_starts.add(start)
                transcript = line_dict.get('t')
                if transcript:
                    # 去掉版本号信息，因为annovar注释没有版本号信息
                    line_dict['t'] = line_dict['t'].split('.', 1)[0]
                if line_dict:
                    filter_lst.append(line_dict)

    # raw_filter_lst = [x for x in filter_lst]

    with VariantFile(vcf) as fr:
        # remove header
        for key in fr.header.info.keys():
            if key != 'AAChange_refGene':
                fr.header.info.remove_header(key)

        # 仅保留INFO中AAChange_refGene信息
        with VariantFile(out_prefix+'.vcf', 'w', header=fr.header) as fw:
            for record in fr:
                if record.info['AAChange_refGene'][0] != '.':
                    # for key in record.info.keys():
                    #     if key != 'AAChange_refGene':
                    #         record.info.pop(key)
                    info = record.info['AAChange_refGene']
                    record.info.clear()
                    record.info['AAChange_refGene'] = info

                    # 实施过滤步骤
                    if filter and (not filter_lst):
                        print('成功匹到所有目标突变！')
                        break

                    if filter and filter_lst:
                        chrom = record.chrom
                        pos = record.pos
                        gene = record.info['AAChange_refGene'][0].split(':', 1)[0]
                        # 提前过滤，加速
                        if filter_genes and (gene not in filter_genes):
                            continue
                        if filter_chrs and (chrom not in filter_chrs):
                            continue
                        if filter_starts and (pos not in filter_starts):
                            continue 

                        keep = []
                        matched_ind = set()
                        for each in record.info['AAChange_refGene']:
                            annot = set(each.split(':') + [chrom, pos])
                            found = False
                            for ind, filter_dict in enumerate(filter_lst):
                                # 使用集合相减的方式判断是否匹配
                                if not (set(filter_dict.values()) - annot):
                                    found = True
                                    matched_ind.add(ind)
                                    if each not in keep:
                                        keep.append(each)
                                    else:
                                        print(f'{filter_dict}匹配到的信息之前已经被匹配到过')
                                    
                            if found and ('t' in filter_dict):
                                # 转录本在当前的info中是独一无二的
                                break

                        if not keep:
                            # 如果匹配不到，则跳过当前record
                            continue
                        else:
                            for ind in matched_ind:
                                filter_lst.pop(ind)
                            record.info['AAChange_refGene'] = keep

                    # 筛选出经典转录本的注释
                    if often_dict:
                        for each in record.info['AAChange_refGene']:
                            tmp = each.split(':')
                            trans = often_dict.get(tmp[0])
                            if len(tmp[0]) <= 1:
                                print('Invalid AAChange_refGene:', each)
                            if trans:
                                if (trans.split('.')[0] == tmp[1].split('.')[0]):
                                    record.info['AAChange_refGene'] = [each]
                                    break
                            else:
                                print(f'{tmp[0]} is not recorded in {often_trans}!')
                    fw.write(record)
                else:
                    pass

    if filter:
        if filter_lst:
            with open(f'{out_prefix}.unmatched.list', 'w') as f:
                    _ = [f.write('\t'.join(x.values())+'\n') for x in filter_lst]
        else:
            print('No thing matched!')


def extract_hotspot_from_vcf(vcf, hotspot, exclude_hotspot=None, id_mode='transcript:chgvs', sample_index=1,
                             af_in_info=False, dp_in_info=False):
    """
    仅适用于hotspot和分析结果都用annovar注释的情况，且只关心有AAChange_refGene注释的突变。
    目的：hotspot是用annovar注释的vcf，如果分析结果也用annovar注释, 可以根据AAChange_refGene信息比较两个突变是否为同一个突变
    从而可以从分析结果中抽提出hotspot中的突变，用于后续的性能验证统计。
    注意：如果要根据start|ref|alt等信息判定突变是否相同，输入的所有vcf应该进行'bcftools norm'标准化
    :param vcf: 分析结果vcf，如果id_mode中包含['gene', 'transcript', 'exon', 'chgvs', 'phgvs'],则要求info列有AAChange_refGene信息
    :param hotspot: vcf, annovar注释过的且info列只保留AAChange_refGene信息且仅保留经典转录本的注释。
    :param exclude_hotspot: 格式同hotspot
    :return:
    """
    hots = get_hotspot(hotspot, id_mode) if type(hotspot) != dict else hotspot
    # print(f'There are {len(hots)} mutation in our hotspot database')
    if exclude_hotspot:
        exclude_hots = get_hotspot(exclude_hotspot, id_mode) if type(exclude_hotspot) != dict else exclude_hotspot
        # print(f'There are {len(exclude_hots)} mutations to be excluded from hotspot')
    else:
        exclude_hots = dict()
    id_mode_lst = id_mode.split(':')
    require_ac = {'gene', 'transcript', 'exon', 'chgvs', 'phgvs'} & set(id_mode_lst)
    with VariantFile(vcf) as f:
        sample = list(f.header.samples)[sample_index]
        result = {sample:dict()}
        for record in f:
            # if not list(record.filter)[0] == 'PASS':
            #     continue
            site_info = [record.chrom, str(record.pos), record.id, record.ref, record.alts[0]]
            if 'COSMIC_ID' in record.info:
                if site_info[2]:
                    site_info[2] += ',' + ','.join(record.info['COSMIC_ID'])
                else:
                    site_info[2] = ','.join(record.info['COSMIC_ID'])
            site_dict = dict(zip(['chr', 'start', 'id', 'ref', 'alt'], site_info))
            if af_in_info:
                af = round(record.info['AF'][0], 4)
            else:
                af = round(record.samples[sample]['AF'][0], 4)

            if dp_in_info:
                depth = record.info['DP']
            else:
                depth = round(record.samples[sample]['DP'], 4)

            if require_ac:
                if not 'AAChange_refGene' in record.info or record.info['AAChange_refGene'][0] == '.':
                    continue
                for hgvs in record.info['AAChange_refGene']:
                    # NOTCH2:NM_024408:exon34:c.6795T>G:p.N2265K
                    # sample = list(record.samples)[1]
                    hgvs_dict = dict(zip(['gene', 'transcript', 'exon', 'chgvs', 'phgvs'], hgvs.split(':')))
                    info_dict = dict(site_dict, **hgvs_dict)
                    key = None
                    try:
                        key = ':'.join(info_dict[x] for x in id_mode_lst)
                    except Exception as e:
                        print(hgvs)
                        print(e)

                    if key in exclude_hots:
                        print(hgvs, f'of {sample} is in low coverage region, and will be excluded')

                    elif key in hots:
                        # fw.write(f'{sample}\t{each}\t{af}\n')
                        if key not in result[sample]:
                            print(f'{sample}' + ':' + hots[key] + ':' + f'AF={af}' + f'DP={depth}')
                            result[sample][key] = [hots[key], af, depth, site_dict['id'], ':'.join(site_info[:2])]
                        else:
                            print(f'{key} is duplicated!')
            else:
                key = ':'.join(site_dict[x] for x in id_mode_lst)
                if key in exclude_hots:
                    print(key, f'of {sample} is in low coverage region, and will be excluded')
                elif key in hots:
                    if key not in result[sample]:
                        # print(f'{sample}' + ':' + hots[key] + ':' + f'AF={af}' + ':' + f'DP={depth}')
                        result[sample][key] = [hots[key], af, depth, site_dict['id'], ':'.join(site_info[:2])]
                    else:
                        print(f'{key} is duplicated!')
    return result


def batch_extract_hotspot(vcfs:list, hotspot, sample_info=None, id_mode='transcript:chgvs',
                          exclude_hotspot=None, cmp_vcfs:list=None, col_names:list=None,
                          prefix="all.detected.hotspot", sample_index=1, cmp_pair=None,
                          af_in_info=False, dp_in_info=False):
    """
    :param vcfs: vcf 路径, 可以提供多个,每个样本对应一个vcf，样本名将从vcf中提取获得.
        分析结果vcf，如果id_mode中包含['gene', 'transcript', 'exon', 'chgvs', 'phgvs'],则要求info列有AAChange_refGene信息.
    :param hotspot: 热点突变文件路径, 使用simplify_annovar_vcf简化的vcf文件
    :param sample_info: 样本信息文件，第一列必须是样本id，第一行是header, 即提供样本注释信息, 注释信息列数不限
    :param id_mode: 指定mutation的唯一id格式, 默认为'transcript:chgvs'.
        支持的字段有 ['chr', 'start', 'id', 'ref', 'alt', 'gene', 'transcript', 'exon', 'chgvs', 'phgvs']
    :param exclude_hotspot: 要排除的热点突变路径，默认不提供
    :param cmp_vcfs: vcf路径，与vcfs提供的数量应该一致，该参数是为了比较两种分析的结果是否一致而设计的, 他们涉及的样本必须一致。默认不提供。
    :param col_names: 项目名称，体现在分析结果的header中，如果提供cmp_vcfs，则需提供两个名称，空格分开即可
    :param prefix: 输出结果文件的前缀，后续会自动加上.xls
    :param sample_index: 指示vcf中第几个样本是肿瘤样本，即目标样本
    :param af_in_info: 指示af信息是否在INFO列，如果提供该参数，则从INFO列中AF提取信息, 该参数的设置是发现罗氏的vcf出现这种情况
    :param dp_in_info: 指示dp信息是否在INFO列，如果提供该参数，则从INFO列中AF提取信息, 该参数的设置是发现罗氏的vcf出现这种情况
    :param cmp_pair: 比较信息，即第一组vcf与第二组vcf的配对信息，第一列为第一组vcf中的样本名，第二列为第二组vcf中的样本名
    :return: 默认生成all.detected.hotspot.xls文件，'sample\tmutation\tAF1\tAF2\tconsistency\n'
    """
    hots = get_hotspot(hotspot, id_mode)
    print(f'There are {len(hots)} mutation in our hotspot database')
    if exclude_hotspot:
        exclude_hots = get_hotspot(exclude_hotspot, id_mode)
        print(f'There are {len(exclude_hots)} mutations to be excluded from hotspot')
    else:
        exclude_hots = dict()

    result = dict()
    for vcf in vcfs:
        if os.path.getsize(vcf) <= 2:
            print('Empty', vcf)
            continue
        var_dict = extract_hotspot_from_vcf(vcf, hots, exclude_hots, id_mode, sample_index, af_in_info, dp_in_info)
        sample = list(var_dict.keys())[0]
        if sample in result:
            result[sample].update(var_dict[sample])
        else:
            result.update(var_dict)

    out_file = prefix + '.xls'
    if not cmp_vcfs:
        with open(out_file, 'w') as fw:
            if col_names:
                name = '_' + col_names[0]
            else:
                name = ''
            fw.write(f'sample\tmutation\tAF{name}(%)\tDepth{name}\tID\tsite\n')
            for sample, var_dict in result.items():
                if not var_dict:
                    print(f'WARN: No hotspot mutation detected in {sample}')
                    fw.write(f'{sample}\tNone\t0.00%\t0\tNone\tNone\n')
                for m_id, detail in var_dict.items():
                    fw.write(f'{sample}\t{detail[0]}\t{detail[1]:.2%}\t{detail[2]}\t{detail[3]}\t{detail[4]}\n')

        if sample_info:
            sample_info_df = pd.read_csv(sample_info, header=0, index_col=0, sep=None, engine='python')
            mutation_df = pd.read_csv(out_file, header=0, index_col=0, sep=None, engine='python')
            mutation_df = mutation_df.round(5)
            final_df = mutation_df.join(sample_info_df)
            final_df.index.name = 'sample'
            final_df.to_excel(out_file[:-3]+'xlsx')
            final_df.to_csv(out_file, sep='\t')
        return result

    else:
        cmp_result = dict()
        for vcf in cmp_vcfs:
            if os.path.getsize(vcf) <= 2:
                print('Empty', vcf)
                continue
            var_dict = extract_hotspot_from_vcf(vcf, hots, exclude_hots, id_mode, sample_index, af_in_info, dp_in_info)
            sample = list(var_dict.keys())[0]
            if sample in cmp_result:
                cmp_result[sample].update(var_dict[sample])
            else:
                cmp_result.update(var_dict)
        # 形成比较信息
        if cmp_pair is None:
            intersection = result.keys() & cmp_result.keys()
            if not intersection:
                raise Exception('两组vcf没有共同的样本，你也没有提供cmp_pair信息，无法开展比较！')
            else:
                print(f'共有{len(intersection)}个交集样本, 将对这些样本进行比较!')
                pair_dict = dict(zip(intersection, intersection))
        else:
            pair_dict = dict(x.strip().split()[:2] for x in open(cmp_pair))

    if col_names:
        name, name2 = '_' + col_names[0], '_' + col_names[1]
        out_file = prefix+f'.{col_names[0]}_vs_{col_names[1]}' + '.xls'
    else:
        name, name2 = '1', '2'
    with open(out_file, 'w') as fw:
        if cmp_pair:
            fw.write(f'ID{name}\tID{name2}\tmutation\tAF{name}(%)\tAF{name2}(%)\tconsistent\tDepth{name}\tDepth{name2}\n')
        else:
            fw.write(f'ID\tmutation\tAF{name}(%)\tAF{name2}(%)\tconsistent\tDepth{name}\tDepth{name2}\n')
        for sample in pair_dict:
            var_dict = result[sample]
            mutations = list(var_dict.keys())
            if pair_dict[sample] in cmp_result:
                mutations += list(cmp_result[pair_dict[sample]].keys())
                var_dict2 = cmp_result[pair_dict[sample]]
            else:
                print(list(cmp_result.keys()))
                raise Exception(f'{pair_dict[sample]} 不在cmp_vcfs中出现, 请核查')
            if not var_dict and (not var_dict2):
                print(f'WARN: No hotspot mutation detected in {sample} and {pair_dict[sample]}')
                if cmp_pair:
                    fw.write(f'{sample}\t{pair_dict[sample]}\tNone\tNone\tNone\tyes\n')
                else:
                    fw.write(f'{sample}\tNone\tNone\tNone\tyes\n')
            for mutation in set(mutations):
                if mutation in var_dict:
                    name = var_dict[mutation][0]
                    af = var_dict[mutation][1]
                    dp = var_dict[mutation][2]
                else:
                    af = 0
                    dp = 'None'
                if mutation in var_dict2:
                    name = var_dict2[mutation][0]
                    af2 = var_dict2[mutation][1]
                    dp2 = var_dict2[mutation][2]
                else:
                    af2 = 0
                    dp2 = 'None'
                if mutation in var_dict and (mutation in var_dict2):
                    consistent = 'yes'
                else:
                    consistent = 'no'
                if not cmp_pair:
                    fw.write(f'{sample}\t{name}\t{af:.2%}\t{af2:.2%}\t{consistent}\t{dp}\t{dp2}\n')
                else:
                    fw.write(f'{sample}\t{pair_dict[sample]}\t{name}\t{af:.2%}\t{af2:.2%}\t{consistent}\t{dp}\t{dp2}\n')
    if sample_info:
        sample_info_df = pd.read_csv(sample_info, header=0, index_col=0, sep=None, engine='python')
        mutation_df = pd.read_csv(out_file, header=0, index_col=0, sep='\t')
        # mutation_df = mutation_df.round(5)
        final_df = mutation_df.join(sample_info_df)
        final_df.index.name = 'sample'
        final_df.to_excel(out_file[:-3]+'xlsx')
        final_df.to_csv(out_file, sep='\t')


def model_statistic(true_positive, false_positive, false_negative, true_negative, confint=False):
    total = true_positive + false_positive + false_negative + true_negative
    overall_accuracy = round((true_positive+true_negative)/total*100, 2)
    if true_positive + false_positive == 0:
        positive_accuracy = 100.0
    else:
        positive_accuracy = round(true_positive/(true_positive+false_positive)*100, 2)

    negative_accuracy = round(true_negative/(true_negative+false_negative)*100, 2)

    if true_positive + false_negative == 0:
        sensitive = 100.0
    else:
        sensitive = round(true_positive/(true_positive+false_negative)*100, 2)

    specificity = round(true_negative/(true_negative+false_positive)*100, 2)
    
    percent = overall_accuracy, positive_accuracy, negative_accuracy, sensitive, specificity
    if confint:
        oa_conf = pconf(true_positive+true_negative, total, method='wilson')
        pa_conf = pconf(true_positive, true_positive+false_positive, method='wilson')
        na_conf = pconf(true_negative, true_negative+false_negative, method='wilson')
        se_conf = pconf(true_positive, true_positive+false_negative, method='wilson')
        sp_conf = pconf(true_negative, true_negative+false_positive, method='wilson')
        confs = oa_conf, pa_conf, na_conf, se_conf, sp_conf
        return percent, confs
    else:
        return percent 


def extract_hotspot_from_avenio_result(infile, hotspot, out, var_id_mode='gene:chgvs', af_is_percent=False):
    # hotdict = dict()
    # with open(hotspot) as f:
    #     header = f.readline().strip().split()
    #     header = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'gene', 'type', 'chgvs', 'phgvs', 'transcript']
    #     for line in f:
    #         line = dict(zip(header, line.strip().split()))
    #         var_id = ':'.join(line[x] for x in var_id_mode.split(':'))
    #         if var_id in hotdict:
    #             raise Exception(f"在{hotspot}中, 根据当前vard_id={var_id}作为突变id时，突变不唯一!")
    #         hotdict[var_id] = ':'.join([line['gene'], line['transcript'], line['chgvs'], line['phgvs']])
    hotdict = get_hotspot(hotspot, id_mode=var_id_mode)

    raw_table = pd.read_csv(infile, sep='\t', header=0)
    targets = ['Sample ID', 'Gene', 'Transcript', 'Coding Change', 'Amino Acid Change', 'Allele Fraction']
    mutations = raw_table.loc[:, targets]
    mutations.columns = ['Sample', 'Gene', 'Transcript', 'c.HGVS', 'p.HGVS', 'AF']
    af_list = []
    for each in mutations['AF']:
        if each:
            each = str(each)
            if '%' in each:
                af_list.append(str(float(each[:-1])/100))
            else:
                if af_is_percent:
                    af_list.append(str(float(each)*0.01))
                else:
                    af_list.append(str(float(each)))
        else:
            af_list.append('0')
    mutations['AF'] = af_list
    out_rows = []
    with open(out, 'w') as f:
        f.write('sample\tmutation\tAF\n')
        var_id_set = set()
        for ind, row in mutations.iterrows():
            if row[2] is np.nan or row[3] is np.nan:
                continue
            is_hot = False
            if row[4] is np.nan:
                row[4] = 'N/A'

            for t, c, p in zip(row[2].split(';'), row[3].split(';'), row[4].split(';')):
                line = dict(zip(['gene', 'transcript', 'chgvs', 'phgvs'], [row[1], t, c, p]))
                var_id = ':'.join(line[x] for x in var_id_mode.split(':'))
                if var_id in var_id_set:
                    raise Exception(f"在{infile}中, 根据当前vard_id={var_id}作为突变id时，突变不唯一!")
                else:
                    var_id_set.add(row[0]+var_id)
                if var_id in hotdict:
                    is_hot = True
                    # variant = ':'.join([row[1], t, 'exon?', c, p])
                    variant = hotdict[var_id]
                    new_line = f'{row[0]}\t{variant}\t{row[5]}\n'
                    if new_line not in out_rows:
                        f.write(new_line)
                        out_rows.append(new_line)
                    else:
                        print('duplicated record:', new_line.strip())
                    # 匹配到一个则跳过剩下的
                    break
            # if not is_hot:
            #     print('not hot:', row)
    return out
 

def parse_formated_mutation(detected, af_cutoff, var_id_mode='transcript:chgvs'):
    """
    :param detected: 文件路径，至少三列信息，示例如下, 如果AF已经是百分比, 需要在header中添加百分比标致:
        sample  mutation        AF
        EPS19F1X1L7     CASP8:NM_001228:exon9:c.904G>C:p.D302H  0.07
    :af_cutoff: af 阈值
    """
    with open(detected, 'r') as f:
        header = f.readline().strip().split()
        if '%' in header[2]:
            af_is_percent = True
        else:
            af_is_percent = False

        dec_dict = dict()
        full_dec_dict = dict()
        for line in f:
            lst = line.strip().split('\t')[:3]
            sample, af = lst[0], lst[2]
            if af.endswith('%'):
                af = float(af[:-1])/100
            else:
                if af_is_percent:
                    af = float(af)/100
                else:
                    af = float(af)
            af = round(af, 5)
            if lst[1] != 'None':
                gene, transcript, exon, chgvs, phgvs = lst[1].split(':')
                transcript = transcript.split('.')[0]
                tmp = dict(zip(
                    ['gene', 'transcript', 'exon', 'chgvs', 'phgvs'],
                    [gene, transcript, exon, chgvs, phgvs]
                ))
                var_id = ':'.join(tmp[x] for x in var_id_mode.split(':'))
                dec_dict.setdefault(sample, dict())
                lst[1] = ':'.join([gene, transcript, chgvs, phgvs])
                if af >= af_cutoff:
                    dec_dict[sample][var_id] = (lst[1], f'{af:.2%}')
                if var_id in dec_dict[sample] and af < af_cutoff:
                    print(f'最后一个突变信息{var_id}覆盖掉前面的突变信息')
                    dec_dict[sample].pop(var_id)
                full_dec_dict.setdefault(sample, dict())
                full_dec_dict[sample][var_id] = (lst[1], f'{af:.2%}')
            else:
                dec_dict.setdefault(sample, dict())
                full_dec_dict.setdefault(sample, dict())
    return dec_dict, full_dec_dict


def overall_stat(detected, known, var_num:int, sample_info, date_col='PCR1完成时间', operator_col='PCR1操作者',
                 replicate_design=None, group_sample=None, lod_group:list=None, report_false_positive=False,
                 detected_af_cutoff=0.0, known_af_cutoff=0.0, lod_cutoff=0.0, lod_deviation=0.0,
                 prefix='final_stat', var_id_mode='transcript:chgvs', include_lod_for_accuracy=False):
    """
    :param detected: 格式举例, 该文件可由batch_extract_hotspot产生的文件的前3列得到
        sample  mutation        AF
        EPS19F1X1L7     CASP8:NM_001228:exon9:c.904G>C:p.D302H  0.07
    :param known: 格式举例, 格式要求必须和detected一样, 必须包含所有需要分析的样本:
        group  mutation        ExpectedAF
        group_01     CASP8:NM_001228:exon9:c.904G>C:p.D302H  0.07
        group_02     None(表示阴性样本)                       0
    :param var_num: 背景总数，用于计算真阴性
    :param sample_info: 文件路径. 需要header. 必须参数
        第一列为样本id，和detected文件里的样本id一致; 其他列可有可无，
        第二列可以是样本组名，和已知突变known文件里的样本id相对应, 用于指示属于相同组的样本具有相同的已知突变.
        这个信息也可以通过group_sample参数提供,而且group_sample优先.
    :param date_col: 指定sample_info中哪一列记录时pcr完成时间信息，用于重复性统计
    :param operator_col: 指定sample_info中哪一列记录操作人员信息，用于重复性 统计
    :param group_sample: 文件路径. 默认无. 该文件用于指示哪些样本共享相同的已知突变
        第一列是样本组名，和已知突变known文件里的样本id相对应, 如果一组有多个样本，可以用';'分割, 也可以采用多行表示分组信息
        第二列为样本id，和detected文件里的样本id一致;
    :param replicate_design: 文件路径, 两列, 第一列是组名; 第二列是样本id, 可以用';'分割的多个样本, 也可以采用多行表示分组信息。
        增加该参数的目的是以防一个样本同时参与LOD和input等复杂设计。 如果不提供该参数，则用样本分组信息替代。需要header
        group       sample_id
        group_lod   EPS19F1X1L7
        group_lod   EPS19F1X1L8
        group_lod   EPS19F1X1L9
    :param lod_group: 组名,指定哪一组是LOD设计样本, 可以赋予多个值, 空格隔开即可，组名必须是relicate_design或sample_info中的
    :param report_false_positive: 如果设置该参数，那么在重复性/LOD统计时，也会对不在已知突变列表的突变进行统计。
    :param detected_af_cutoff: 实际检测到的突变的af阈值, 用于准确性, 敏感度统计所需，该值不影响LOD分析，默认0.0001
    :param known_af_cutoff: 已知突变的af阈值，当af阈值>=设定的值时才认为是有效已知突变，默认0.02
    :param lod_cutoff: 对于lod的分组样本，对已知突变进行>=lod_cutoff过滤
    :param lod_deviation: 百分比，配合lod_cutoff使用，对于lod的分组样本，对检测到的突变进行>=lod_cutoff*(1- lod_deviation)过滤.
        如果设置为1，则相当于不对检测结果进行过滤
    :param prefix: 准确性统计表的前缀，输出结果会在文件名后面加xls
    :param var_id_mode :指定mutation的唯一id格式, 默认为'transcript:chgvs'.
    :param include_lod_for_accuracy: 如果提供该参数, 统计准确性时需要把LOD设计样本包含进来，默认为False，即统计时排除lod设计的样本
    :return: 输出多个文件，其中{prefix}.xlsx文件为最主要的结果，其中包含多个sheet，几乎囊括所有分析结果.
    """
    outdir = os.path.dirname(prefix)
    lod_groups = lod_group if lod_group is not None else []
    dec_dict, full_dec_dict = parse_formated_mutation(detected, detected_af_cutoff, var_id_mode=var_id_mode)
    known_dict, full_known_dict = parse_formated_mutation(known, known_af_cutoff, var_id_mode=var_id_mode)

    # 根据样本信息，对已知突变信息进行拓展，使其中的键表示样本而不再是group
    if group_sample:
        mutation_group = dict()
        with open(group_sample, 'r') as f:
            for line in f:
                share, samples = line.strip().split()[:2]
                sample_lst = samples.split(';')
                for sample in sample_lst:
                    mutation_group.setdefault(share, list())
                    mutation_group[share].append(sample)
    else:
        col2 = pd.read_csv(sample_info, sep=None, header=0, engine='python').iloc[:, [1, 0]]
        col2.columns = ['mgroup', 'sample']
        mutation_group = col2.groupby('mgroup')['sample'].apply(list).to_dict()

    new_known_dict = dict()
    for g, d in known_dict.items():
        for sample in mutation_group[g]:
            new_known_dict[sample] = d
    known_dict = new_known_dict

    new_full_known_dict = dict()
    for g, d in full_known_dict.items():
        for sample in mutation_group[g]:
            new_full_known_dict[sample] = d
    full_known_dict = new_full_known_dict

     # 更新dec_dict, 把没有检测结果的样本包含进来
    no_detect_samples = set(known_dict.keys()) - set(dec_dict.keys())
    for each in no_detect_samples:
        print(f'{each} has no detected result')
        dec_dict[each] = dict()
        full_dec_dict[each] = dict()
    
    # 对于那些由于频率低于cutoff而过滤掉的已知突变，即使检测出来了，也要当作没有检测出来
    new_dec_dict = dict()
    possible_LOD = []
    for sample, var_dict in dec_dict.items():
        if sample not in full_known_dict:
            raise Exception(f"查不到{sample}的已知突变信息，你是否忘了将其列在样本信息表中？")
        new_dec_dict[sample] = dict()
        for var_id, detail in var_dict.items():
            if var_id not in known_dict[sample] and (var_id in full_known_dict[sample]):
                print(sample, detail, '因为低频被从已知突变中过滤，因此即使检测出也将忽略掉')
                possible_LOD.append([sample, detail[0], detail[1], full_known_dict[sample][var_id][1]])
            else:
                new_dec_dict[sample][var_id] = detail
    dec_dict = new_dec_dict
    if possible_LOD:
        with open(os.path.join(outdir, 'Low_expected_AF_but_detected.xls'), 'w') as f:
            possible_LOD.sort(key=lambda x:x[3])
            f.write('sample\tmutation\tdetected\texpected\n')
            for sample, mut, af, af2 in possible_LOD:
                f.write(f'{sample}\t{mut}\t{af}\t{af2}\n')

    # 输出新的检测结果和已知结果
    with open(os.path.join(outdir, 'detected.muation.xls'), 'w') as f:
        f.write('sample\tmutation\tAF(%)\n')
        for sample, var_dict in dec_dict.items():
            for _, (var_detail, af) in var_dict.items():
                f.write(f'{sample}\t{var_detail}\t{af}\n')

    with open(os.path.join(outdir, 'known.muation.xls'), 'w') as f:
        f.write('sample\tmutation\tAF(%)\n')
        for sample, var_dict in known_dict.items():
            for _, (var_detail, af) in var_dict.items():
                f.write(f'{sample}\t{var_detail}\t{af}\n')
    
    # 解析重复设计信息
    if replicate_design is None:
        replicate_dict = mutation_group
    else:
        with open(replicate_design, 'r') as f:
            replicate_dict = dict()
            for line in f:
                share, samples = line.strip().split()[:2]
                sample_lst = samples.split(';')
                for sample in sample_lst:
                    replicate_dict.setdefault(share, list())
                    replicate_dict[share].append(sample)
    replicate_dict = {x:y for x, y in replicate_dict.items() if len(y) > 1}
    
    lod_samples = []
    for lod_group in lod_groups:
        if lod_group in replicate_dict:
            lod_samples += replicate_dict[lod_group]
        else:
            raise Exception(f'{lod_group} 不在重复设计信息当中')

    summary = list()
    concordance = list()
    accuracy = list()
    summary_header = [
        'Sample',
        '真阳性', '假阳性', '真阴性', '假阴性',
        '总符合率', '阳性符合率', '阴性符合率',
        '敏感性', '特异性',
        '真阳性位点', '真阳性频率',
        '假阳性位点', '假阳性频率',
        '假阴性位点假', '阴性频率'
    ]
    concordance_header = ['Sample',	'Mutation', 'Expected AF(%)', 'Detected AF(%)', 'Concordance']
    accuracy_header = ['Sample', 'TP', 'FP', 'TN', 'FN']

    for sample, var_dict in dec_dict.items():
        if not include_lod_for_accuracy and (sample in lod_samples):
            print(f"LOD设计样本{sample}不参与准确性灵敏度等统计")
            continue
        for mid in (set(var_dict.keys()) | set(known_dict[sample].keys())):
            # mutation = known_dict[sample][mid][0] if mid in known_dict[sample] else var_dict[mid][0]
            mutation = var_dict[mid][0] if mid in var_dict else known_dict[sample][mid][0]
            expected_af = known_dict[sample][mid][1] if mid in known_dict[sample] else 'Unexpected'
            detected_af = var_dict[mid][1] if mid in var_dict else 'Not detected'
            consistent = 'Yes' if (mid in var_dict) and (mid in known_dict[sample]) else 'No'
            concordance.append([sample, mutation, expected_af, detected_af, consistent])

        tp_ids = set(var_dict.keys()) & set(known_dict[sample].keys())
        true_positive = len(tp_ids)
        fp_ids = set(var_dict.keys()) - set(known_dict[sample].keys())
        false_positive = len(fp_ids)
        fn_ids = set(known_dict[sample].keys()) - tp_ids
        false_negative = len(fn_ids)
        true_negative = var_num - len(known_dict[sample]) - false_positive
        overall_accuracy, positive_accuracy, negative_accuracy, sensitive, specificity = model_statistic(
            true_positive, false_positive, false_negative, true_negative)
        tp_mutations = [var_dict[x][0] for x in tp_ids]
        tp_mutations_af = [var_dict[x][1] for x in tp_ids]
        fp_mutations = [var_dict[x][0] for x in fp_ids]
        fp_mutations_af = [var_dict[x][1] for x in fp_ids]
        fn_mutations = [known_dict[sample][x][0] for x in fn_ids]
        fn_mutations_af = []
        for mid in fn_ids:
            if mid in full_dec_dict[sample]:
                fn_mutations_af.append(full_dec_dict[sample][mid][1])
            else:
                fn_mutations_af.append(0)
        # fn_mutations_af = [known_dict[sample][x][1] for x in fn_ids]
        lst = [
            sample,
            true_positive, false_positive, true_negative, false_negative,
            overall_accuracy, positive_accuracy, negative_accuracy,
            sensitive, specificity,
            tp_mutations, tp_mutations_af,
            fp_mutations, fp_mutations_af,
            fn_mutations, fn_mutations_af
        ]
        summary.append(lst)
        accuracy.append([sample, true_positive, false_positive, true_negative, false_negative])

    out = prefix + '.xls'
    with pd.ExcelWriter(out[:-3]+'xlsx') as writer:
        result = pd.DataFrame(summary, columns=summary_header)
        result = result.set_index('Sample')
        true_positive = result['真阳性'].sum()
        false_positive = result['假阳性'].sum()
        true_negative = result['真阴性'].sum()
        false_negative = result['假阴性'].sum()
        percent, confs = model_statistic(true_positive, false_positive, false_negative, true_negative, confint=True)
        percent_confs = [f'{x}%'+'(['+f'{y[0]:.2%}'+','+f'{y[1]:.2%}'+'])' for x,y in zip(percent, confs)]
        result.loc['Total'] = [true_positive, false_positive, true_negative, false_negative] + percent_confs + ['[]']*6
        # 添加样本信息
        sample_df = pd.read_csv(sample_info, sep=None, header=0, index_col=0, engine='python')
        sample_df.index.name = 'sample'
        result = sample_df.join(result, how='right')
        # 输出结果
        result.to_csv(out, sep='\t')
        result.to_excel(writer, sheet_name='Summary')

        concordance_df = pd.DataFrame(concordance, columns=concordance_header)
        mut_cols = concordance_df['Mutation'].str.split(':', expand=True)
        mut_cols.columns = ['Gene', 'Transcript', 'cHgvs', 'pHgvs']
        mut_cols = concordance_df.loc[:, ['Sample']].join(mut_cols)
        concordance_df = mut_cols.join(concordance_df.iloc[:, 2:])
        concordance_df = concordance_df.set_index(['Sample', 'Gene', 'Transcript', 'cHgvs', 'pHgvs'])
        concordance_df.to_excel(writer, sheet_name='Concordance')

        accuracy_df = pd.DataFrame(accuracy, columns=accuracy_header)
        accuracy_df.set_index('Sample', inplace=True)
        accuracy_df.loc['total'] = accuracy_df.sum()
        accuracy_df.to_excel(writer, sheet_name='Accuracy')

        lod_detect_dict = dict()
        for sample, detail in full_dec_dict.items():
            lod_detect_dict[sample] = {
                x:y for x, y in detail.items() if float(y[1][:-1])*0.01 >= lod_cutoff*(1-lod_deviation)
            }

        lod_known_dict = dict()
        for sample, detail in full_known_dict.items():
            lod_known_dict[sample] = {
                x: y for x, y in detail.items() if float(y[1][:-1])*0.01 >= lod_cutoff
            }

        # replicate_stat
        for group, samples in replicate_dict.items():
            if group not in lod_groups:
                summary_df, concordance_df = replicate_stat(samples, group, dec_dict, known_dict, report_false_positive)
                day_operator_df = sample_df.loc[:, [date_col, operator_col]]
                day_operator_df.columns = ['Day', 'Operator']
                concordance_df = concordance_df.set_index('Sample').join(day_operator_df)
                concordance_df.index.name = 'Sample'
                concordance_df.reset_index(inplace=True)
                concordance_df.sort_values(by=['Mutation', 'Expected AF(%)', 'Day', 'Operator', 'Sample'], inplace=True)
                concordance_df.set_index(['Mutation', 'Expected AF(%)', 'Day', 'Operator', 'Sample'], inplace=True)
                concordance_df.to_csv(os.path.join(outdir, f'{group}.replicate.xls'), sep='\t')
                summary_df.to_csv(os.path.join(outdir, f'{group}.replicate.summary.xls'), sep='\t', index=False)
                # to sheet
                concordance_df.reset_index(inplace=True)
                mut_cols = concordance_df['Mutation'].str.split(':', expand=True)
                mut_cols.columns = ['Gene', 'Transcript', 'cHgvs', 'pHgvs']
                concordance_df = mut_cols.join(concordance_df.iloc[:, 1:])
                concordance_df.set_index(
                    ['Gene', 'Transcript', 'cHgvs', 'pHgvs', 'Expected AF(%)', 'Day', 'Operator', 'Sample'],
                    inplace=True
                )

                concordance_df.to_excel(writer, sheet_name=f'{group}.Rep')

                mut_cols = summary_df['Mutation'].str.split(':', expand=True)
                mut_cols.columns = ['Gene', 'Transcript', 'cHgvs', 'pHgvs']
                mut_cols = summary_df.loc[:, ['Sample']].join(mut_cols)
                summary_df = mut_cols.join(summary_df.iloc[:, 2:])

                summary_df.to_excel(writer, sheet_name=f'{group}.RepSum', index=False)
            else:
                summary_df, concordance_df = replicate_stat(
                    samples, group, lod_detect_dict, lod_known_dict, report_false_positive
                )
                order = ['Mutation', 'Expected AF(%)', 'Sample', 'Detected AF(%)', 'Concordance']
                concordance_df = concordance_df.loc[:, order]
                concordance_df.sort_values(by=['Mutation', 'Expected AF(%)', 'Sample'], inplace=True)
                concordance_df.set_index(['Mutation', 'Expected AF(%)', 'Sample'], inplace=True)
                concordance_df.to_csv(os.path.join(outdir, f'{group}.LOD.xls'), sep='\t')
                summary_df.to_csv(os.path.join(outdir, f'{group}.LOD.summary.xls'), sep='\t', index=False)
                # to sheet
                concordance_df.reset_index(inplace=True)
                mut_cols = concordance_df['Mutation'].str.split(':', expand=True)
                mut_cols.columns = ['Gene', 'Transcript', 'cHgvs', 'pHgvs']
                concordance_df = mut_cols.join(concordance_df.iloc[:, 1:])
                concordance_df.set_index(
                    ['Gene', 'Transcript', 'cHgvs', 'pHgvs', 'Expected AF(%)', 'Sample'],
                    inplace=True
                )

                concordance_df.to_excel(writer, sheet_name=f'{group}.LOD')

                mut_cols = summary_df['Mutation'].str.split(':', expand=True)
                mut_cols.columns = ['Gene', 'Transcript', 'cHgvs', 'pHgvs']
                mut_cols = summary_df.loc[:, ['Sample']].join(mut_cols)
                summary_df = mut_cols.join(summary_df.iloc[:, 2:])
                summary_df.to_excel(writer, sheet_name=f'{group}.LODSum', index=False)

    # new_lod_stat
    for group, samples in replicate_dict.items():
        if group in lod_groups:
            out = os.path.join(outdir, group+'.Gradient.LOD.xls')
            new_lod_stat(samples, lod_detect_dict, lod_known_dict, out=out, gradient=(0.01, 0.02, 0.03, 0.04, 0.05))

    samples = tuple(full_known_dict.keys())
    out = os.path.join(outdir, 'Overall.Gradient.LOD.xls')
    new_lod_stat(samples, full_dec_dict, full_known_dict, out=out, gradient=(0.01, 0.02, 0.03, 0.04, 0.05))


def replicate_stat(samples, group_name, detected_dict, known_dict, report_false_positive=False):
    """
    以 突变为单位的 一致性/检出率 统计
    根据report_false_positive参数决定是考察已知突变还是考察所有检测到的突变
    假设这里的samples都具有相同的已知突变
    """
    if report_false_positive:
        detected_mid = set(x for s, d in detected_dict.items() for x in d if s in samples)
        known_mid = set(known_dict[samples[0]].keys())
        false_positive_mid = detected_mid - known_mid
        target_mutations = detected_mid | known_mid
    else:
        target_mutations = known_dict[samples[0]].keys()
        false_positive_mid = set()

    concordance = []
    concordance_header = ['Sample',	'Mutation', 'Expected AF(%)', 'Detected AF(%)', 'Concordance']

    sample_num = len(samples)
    summary = []
    summary_header = [
        'Sample', 'Mutation', 'MAF_Range', 'MAF_Mean', 'MAF_median',
        'MAF(SD)', 'MAF(%CV)', 'Positive_Calls', 'Positive_Calls_Rate'
    ]

    for mid in target_mutations:
        positive = 0
        af_lst = []
        mutation = None
        for sample in samples:
            if sample not in detected_dict:
                print(sample, 'not in detected_dict')
                continue
            if (mid not in detected_dict[sample]) and (mid not in known_dict[sample]):
                continue
            # concordance
            # mutation = known_dict[sample][mid][0] if mid in known_dict[sample] else detected_dict[sample][mid][0]
            mutation = detected_dict[sample][mid][0] if mid in detected_dict[sample] else known_dict[sample][mid][0]
            expected_af = known_dict[sample][mid][1] if mid in known_dict[sample] else 'Unexpected'
            detected_af = detected_dict[sample][mid][1] if mid in detected_dict[sample] else 'Not detected'
            consistent = 'Yes' if (mid in detected_dict[sample]) and (mid in known_dict[sample]) else 'No'
            concordance.append([sample, mutation, expected_af, detected_af, consistent])
            # positive_stat
            if mid in detected_dict[sample]:
                positive += 1
                af_lst.append(detected_dict[sample][mid][1])

        if not af_lst:
            af_lst = ['0', '0']
        if len(af_lst) == 1:
            # 一个数据无法计算标准差
            af_lst.append(af_lst[0])
        af_lst = [float(x.strip('%')) * 0.01 if x.endswith('%') else float(x) for x in af_lst]
        af_lst = [round(x, 5) for x in af_lst]
        # af_range = str(min(af_lst)) + '-' + str(max(af_lst))
        af_range = f'{min(af_lst):.2%}' + '-' + f'{max(af_lst):.2%}'
        af_mean = statistics.mean(af_lst)
        af_median = statistics.median(af_lst)
        af_std = statistics.stdev(af_lst)
        af_cv = af_std / af_mean if af_mean else 0
        ratio = f'{positive}|{sample_num}'
        confint = pconf(positive, sample_num, method='wilson')
        rate = positive / sample_num
        call = f'{rate:.2%}([{confint[0]:.2%}, {confint[1]:.2%}])'
        group = group_name +'|Unexpected' if mid in false_positive_mid else group_name
        summary.append([
            group, mutation, af_range, f'{af_mean:.2%}', f'{af_median:.2%}',
            f'{af_std:.2%}', f'{af_cv:.2%}', ratio, call
        ])

    summary_df = pd.DataFrame(summary, columns=summary_header)
    concordance_df = pd.DataFrame(concordance, columns=concordance_header)
    return summary_df, concordance_df


def new_lod_stat(samples, detected_dict, known_dict, out=None, gradient=(0.01, 0.02, 0.03, 0.04, 0.05)):
    # 对突变按照理论频率进行梯度分组如>1,2,3,4,5
    grad_dict = dict()
    gradient = [f'{x:.2%}' for x in gradient]
    for sample, var_dict in known_dict.items():
        if sample in samples:
            for k, v in var_dict.items():
                if v[1] in gradient:
                    grad_dict.setdefault(v[1], list())
                    grad_dict[v[1]].append([sample, k])
    # 根据gradient排序
    grad_dict = {x:grad_dict[x] for x in gradient if x in grad_dict}

    result = []
    for grad, mutations in grad_dict.items():
        expected_num = len(mutations)
        positive = 0
        mutation_names = list()
        af_lst = []
        for sample, mutation in mutations:
            if sample not in detected_dict:
                print(sample, 'not in detected_dict')
                continue
            if mutation in detected_dict[sample]:
                positive += 1
                mutation_name = detected_dict[sample][mutation][0]
                mutation_names.append(mutation_name)
                af_lst.append(detected_dict[sample][mutation][1])

        if not af_lst:
            af_lst = ['0', '0']
        if len(af_lst) == 1:
            # 一个数据无法计算标准差
            af_lst.append(af_lst[0])
        af_lst = [float(x.strip('%')) * 0.01 if x.endswith('%') else float(x) for x in af_lst]
        af_lst = [round(x, 5) for x in af_lst]
        # af_range = str(min(af_lst)) + '-' + str(max(af_lst))
        af_range = f'{min(af_lst):.2%}' + '-' + f'{max(af_lst):.2%}'
        af_mean = statistics.mean(af_lst)
        af_median = statistics.median(af_lst)
        af_std = statistics.stdev(af_lst)
        af_cv = af_std / af_mean if af_mean else 0
        ratio = f'{positive}|{expected_num}'
        confint = pconf(positive, expected_num, method='wilson')
        rate = positive / expected_num
        call = f'{rate:.2%}([{confint[0]:.2%}, {confint[1]:.2%}])'
        result.append([
            str(grad), ';'.join(set(mutation_names)), af_range,
            f'{af_mean:.2%}', f'{af_median:.2%}', f'{af_std:.2%}', f'{af_cv:.2%}',
            ratio, call])

    # summary
    header = [
        'ExpectedAF', 'Mutations', 'MAF_Range',
        'MAF_Mean', 'MAF_median',
        'MAF(SD)', 'MAF(%CV)', 'Positive_Calls', 'Positive_Calls_Rate'
    ]
    pd.DataFrame(result, columns=header).to_csv(out, sep='\t', index=False)
    pd.DataFrame(result, columns=header).to_excel(out + 'x')


def cosmic2vcf(mutation_txt, genome, out='cosmic.vcf', skip_chrom:tuple=('chrMT')):
    genome = pysam.FastaFile(genome)
    vcf = set_logger(out, logger_id='merged')
    for line in vcf_header():
        vcf.info(line)

    with open(mutation_txt) as f:
        for line in f:
            if line.startswith('#'):
                continue

            chr_id, start, end, ref, alt, idx = line.strip().split()
            
            if alt == '-' and ref == '-':
                continue
            if alt == '0' and ref == '0':
                continue

            start = int(start)
            end = int(end)
            if chr_id.isnumeric() or (chr_id in ['X', 'Y', 'M', 'MT']):
                chr_id = 'chr' + chr_id
                if chr_id in skip_chrom:
                    # chr_id = 'chrM'
                    print(f'skip {chr_id}')
                    continue

            if ref == '-':
                ref = genome.fetch(chr_id, start-1, start)
                alt = ref + alt
            elif alt == '-':
                if len(ref) != end - start + 1:
                    raise Exception(f'ref长度和坐标暗示的不一致：{line}')
                start = start-1
                alt = genome.fetch(chr_id, start-1, start)
                ref = alt+ref
            else:
                pass
            vcf.info(f'{chr_id}\t{start}\t{idx}\t{ref}\t{alt}\t.\tPASS\t.')


def create_variant_db(vcfs:list, sources:list, canonical_refseq=None, host='0.0.0.0', port=27017, db_name='variants', coll_name='hgvs'):
    import pymongo
    client = pymongo.MongoClient(host=host, port=port) 
    db = client[db_name]
    coll = db[coll_name]

    if len(vcfs) > 1 and len(sources) == 1:
        sources = sources*len(vcfs)
    elif len(vcfs) == len(sources):
        pass
    else:
        raise Exception("vcfs and sources are not properly matched")

    if canonical_refseq is not None:
        g2t = dict([x.strip().split()[:2] for x in open(canonical_refseq) if len(x.strip().split()) > 1])
    else:
        g2t = dict()

    for vcf, source in zip(vcfs, sources):
        with VariantFile(vcf) as f:
            main_keys = ['source', 'chrom', 'start', 'var_id', 'gene', 'hgvs', 'ref', 'alt', 'canonical']
            hgvs_keys = ['gene', 'transcript', 'exon', 'chgvs', 'phgvs']
            for record in f:
                if record.info['AAChange_refGene'][0] != '.':
                    chrom = record.chrom
                    start = record.pos
                    ref = record.ref
                    alt = record.alts[0]
                    var_id = record.id
                    canonical = 'None'
                    hgvs_values = []
                    for hgvs in record.info['AAChange_refGene']:
                        detail = hgvs.split(':')
                        gene = detail[0]
                        if len(detail) == 5:
                            g, t, e, c, p = detail
                            hgvs_values = detail
                        elif len(detail) == 4:
                            g, t, e, c = detail
                            hgvs_values = [g, t, e, c, 'None']
                        else:
                            print(f'{detail} is out of expectation, the stored info maybe Wrong!')
                            hgvs_values = detail

                        if gene in g2t:
                            ct = g2t[gene]
                            if ct.startswith(detail[1]):
                                canonical = hgvs

                        hgvs_dict = dict(zip(hgvs_keys, hgvs_values))
                        main_values = [source, chrom, start, var_id, gene, hgvs_dict, ref, alt, canonical]
                        coll.insert_one(dict(zip(main_keys, main_values)))
                else:
                    print('skip line with empty "AAChange_refGene"')


def query_variant_db(query, out='query_result.vcf', host='0.0.0.0', port=27017, db_name='variants', coll_name='hgvs'):
    """
    支持的header有 ['gene', 'transcript', 'exon', 'chgvs', 'phgvs']
    """
    import pymongo
    client = pymongo.MongoClient(host=host, port=port) 
    db = client[db_name]
    coll = db[coll_name]

    with open(query) as f, open(out, 'w') as fw:
        for line in vcf_header():
            fw.write(line+'\n')

        header = f.readline().strip().split()
        header = ['hgvs.'+x if x in ('transcript', 'exon', 'chgvs', 'phgvs') else x for x in header]
        new_rows = []
        for line in f:
            condition = dict(zip(header, line.strip().split()))
            docs = list(coll.find(condition))
            if not docs:
                print('failed query:', line.strip())
                continue
            # sources = [doc['source'] for x in docs]
            # mutation_set = {(doc['chrom'], doc['start'], doc['ref'], doc['alt']) for x in docs}
            filtered_docs = [doc for doc in docs if doc["canonical"] != "None"]
            if filtered_docs:
                docs = filtered_docs
            # if len(docs) > 1:
                # print(f'{len(docs)} records for {line.strip()}')
            for doc in docs:
                lst = [doc[x] for x in ['chrom', 'start', 'var_id', 'ref', 'alt']]
                lst.append('.')
                lst.append("PASS")
                query = f'Query={condition}'
                if doc["canonical"] == 'None':
                    # lst.append(f'{query};Source={doc["source"]};AAChange_refGene={doc["hgvs"]}')
                    hgvs = ':'.join(doc['hgvs'].values())
                    lst.append(f'AAChange_refGene={doc["hgvs"]}')
                else:
                    # lst.append(f'{query};Source={doc["source"]};AAChange_refGene={doc["canonical"]}')
                    lst.append(f'AAChange_refGene={doc["canonical"]}')
                new_line = '\t'.join((str(x) for x in lst))+'\n'
                if new_line not in new_rows:
                    fw.write(new_line)
                    new_rows.append(new_line)


def del_db(host='0.0.0.0', port=27017, db_name='variants', coll_name='hgvs', condition=None):
    import pymongo
    client = pymongo.MongoClient(host=host, port=port) 
    if db_name not in client.list_database_names():
        exit(f'databse {db_name} is not found')
    db = client[db_name]

    colls = db.list_collection_names()
    print('Find collections:', colls)
    if coll_name not in colls:
        exit(f'{coll_name} is not found')
    coll = db[coll_name]

    if condition is None:
        print(f'You are dropping the whole collection {coll_name}')
        password = input('Password:')
        if password != 'yes':
            exit('Password is Invalid')
        else:
            coll.drop()
    else:
        query = dict(x.split('=') for x in condition.split(','))
        print(f'You are deleting docs by querying {query}')
        password = input('Password:')
        if password != 'yes':
            exit('Password is Invalid')
        else:
            if len(list(coll.find(query))) > 0:
                coll.delete_many(query)
            else:
                print('Nothing matched and no deletion')


def change_transcript_version(mutalyzer_coor_query_result, out='new.result'):
    with open(mutalyzer_coor_query_result) as f, open(out, 'w') as f2:
        _ = f.readline()
        for line in f:
            lst = line.strip().split('\t')
            if lst[1]:
                err_info = lst[1]
                transcripts = re.findall(r'NM_[0-9]+\.[0-9]+', err_info)
                transcripts.sort(key=lambda x:int(x.split('.')[1]))
                f2.write(lst[0]+'\t'+transcripts[0] + ':' + lst[0].split(':')[1] +'\n')
            else:
                f2.write(lst[0]+'\t'+lst[0]+'\n')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(
        locals(), 
        include=[
            'trans', 
            'cat_vcf',
            'dedup_vcf',
            'translate_to_vcf_5cols',
            'parse_mycancergenome_var', 
            'batch_extract_hotspot',
            'hotspot2msk',
            'overall_stat',
            'simplify_annovar_vcf',
            'cosmic2vcf',
            'create_variant_db',
            'query_variant_db',
            'extract_hotspot_from_avenio_result',
            'del_db',
            'change_transcript_version',
        ]
    )

# import matplotlib
# matplotlib.use('agg')
# import seaborn as sns
# import matplotlib.pyplot as plt
# sns.set(style="white")
# import numpy as np
# import pandas as pd

# col_name = ['sphere1', 'sphere2']
# col_name2 = ['non-sphere1', 'non-sphere2']

# a = np.random.randn(120, 2) + 5
# b = np.random.randn(120, 2) + 2
# a = pd.DataFrame(a, columns=col_name)
# b = pd.DataFrame(b, columns=col_name2)
# c1 = a.join(b)
# a = np.random.randn(80, 2) + 2
# b = np.random.randn(80, 2) + 5
# a = pd.DataFrame(a, columns=col_name)
# b = pd.DataFrame(b, columns=col_name2)
# c2 = a.join(b)
# c = pd.concat([c1, c2])
# g = sns.clustermap(c, z_score=0, cmap='RdBu', yticklabels=False)
# plt.savefig('clustermap.png', dpi=200, bbox_inches='tight')
