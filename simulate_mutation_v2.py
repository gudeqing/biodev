# coding=utf-8
from __future__ import print_function
import pysam
import random
import re
import os
import logging
import time
import gzip
import json
import cPickle as pickle
__author__ = 'gdq'
"""
bam是0-based坐标，pysam也是用0-based方法读取信息
0. fetch命令可以读取某个区域比对上的所有reads，随后可以决定要随机修改的reads集合
1. 用fetch读取bam文件，结合输入的突变位点信息，获得需要修改的reads集合
3. 修改原reads并生成fastq文件，有两种思路：
    a. 按顺序读取bam文件，并一行行写入新的bam文件当中去，在写之前进行判断并修改，
       需要过滤掉secondary-alignment，后直接获得修改后的bam文件。
       最后把bam文件转换为fastq, 这中间需要对bam进行排序（相对第二种，无需输入原fastq，但速度相对慢）
    b. 直接遍历原fastq文件，修改相应reads，生成新的fastq文件（速度较快，需输入原fastq）
4. 弊端：暂时不能模拟大片段的indel，即indel超过read长度;仅针对单端测序数据模拟
"""


def set_logger(name='log.info', logger_id='x'):
    logger = logging.getLogger(logger_id)
    fh = logging.FileHandler(name, mode='w+')
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    return logger


def replace(string, old, new, start=0):
    """
    从指定位置开始，替换字符串string中字符串
    :param string:
    :param old: string after 'start' position to be replaced
    :param new: replacement
    :param start: start index of 'old'
    :return: string after replacing
    """
    # 'xxx:old:yyy' --> 'xxx:new_string:yyyy'
    new_str = string[:start] + new + string[start+len(old):]
    return new_str


def target_mutation_info(in_file, max_indel=120, af=None):
    with open(in_file) as fr:
        for line in fr:
            line = line.strip()
            if not line:
                continue
            if af:
                chr_, pos, x, ref, alt = line.split('\t')[:5]
            else:
                chr_, pos, x, ref, alt, af = line.split('\t')[:6]
            if len(alt) > max_indel or len(ref) > max_indel:
                print("indel length is too large: {}. And we will discard this one.".format(line))
                continue
            yield chr_, int(pos), ref, alt, af


def pickle_read_pos(bam_file):
    # 获得所有paired reads的比对的起始位置和终止位置
    read_pos_dict_file = os.path.join(os.path.dirname(bam_file), 'read_start_end_pos.pkl')
    if not os.path.exists(read_pos_dict_file):
        bam_obj = pysam.AlignmentFile(bam_file, "rb")
        read_pos_dict = dict()
        for read in bam_obj.fetch():
            if (not read.is_secondary) and read.is_paired and (not read.is_unmapped):
                read_name = read.query_name + '_r1' if read.is_read1 else read.query_name + '_r2'
                read_pos_dict[read_name] = (read.reference_start, read.reference_end)
        bam_obj.close()
        with open(read_pos_dict_file, 'wb') as f:
            pickle.dump(read_pos_dict, f)

    else:
        with open(read_pos_dict_file, 'rb') as f:
            read_pos_dict = pickle.load(f)
    return read_pos_dict


def get_candidate_info(bam_file, chr_, ref_start, ref, alt, depth_threshold=300, expected_af=0.1, logger=None):
    """
    read的分类说明：
    # dup-pair-read: read 和 read' 对应的r1和r2的比对起始位置相同,且他们的mate的比对终止位置也相同， 模拟时read和read'需同时修改
    # pair-read： 若r1 and r2 都包含突变位点则属于pair, 模拟时r1和r2需同时修改
    # unique-pair： pair-read 减去 dup-pair-read
    # single：r1和r2中只有一条read比对当前要模拟的区域，由于此时若要同时获得两条read的比对位置，势必要通读bam文件，耗时较长。
    # 因此下面用likely的策略来对read进行分类，以供后续突变模拟的选择使用。
    # likely unique-single： likely指仅用当前read的比对起止位置作为标记来判断是否duplicated
    # likely dup-single： likely指仅用当前read的比对起止位置作为标记来判断是否duplicated
    # bad-read: (not read.is_proper_pair) | read.is_secondary | cigar-string不全为M | 当前突变位点处已经存在突变

    候选read的筛选条件：
    1. 突变起始位置要落在read比对上的区间
    2. read能比对上，且cigar-string不能为空, 且是is-proper—pair，且不能是is-secondary
    3. read1和read2在当前模拟的突变位点处不能存在突变；cigar-string只能全是M
    4. 如read末端仅包含部分ref（突变前的参考序列）序列，则该read及其mate也不能参与模拟，且有日志记录
    5. 如果某组dup read 数量已经达到了目标模拟reads数的85%， 那么不应该把这组dup read作为模拟候选，且有日志记录
    6. bad-read-number的比例占coverage的0.1以上时，不对当前提供的位点进行突变模拟，且日志有记录
    7. coverage低于提供的阈值depth_threshold值时，不对当前提供的位点进行突变模拟，且日志有记录

    :param bam_file: bam 文件的路径
    :param chr_: 染色体名称
    :param ref_start: 突变起始位置
    :param ref: 突变前的参考序列
    :param alt: 突变后的序列
    :param depth_threshold: 测序深度阈值
    :param expected_af: 期望模拟的突变频率
    :param logger: 日志文件对象，默认从头创建日志
    :return:
    """
    if logger is None:
        logger = set_logger('get_candidates.log', 'get_candidates')
    pattern = re.compile('\d+M$')
    locate_start = ref_start - 1
    locate_end = locate_start + len(ref)
    logger.info("Start to screen valid reads for mutation simulation")
    bam_obj = pysam.AlignmentFile(bam_file, "rb")
    read_pos_dict = dict()
    pair_read_dict = dict()
    coverage = 1
    candidates = set()
    bad_pair_set = set()  # 存储那些至少有一端read存在变异的pair，他们需要从candidates中被排除
    bad_read_number = 0
    for ind, read in enumerate(bam_obj.fetch(chr_, locate_start, locate_end)):
        if locate_start not in read.get_reference_positions():
            continue
        if read.cigarstring is None or read.is_unmapped:
            continue
        coverage += 1
        seq = read.seq
        start = ref_start - 1 - read.pos  # ref序列在当前read上的起始位置
        read_name = read.query_name + '_r1' if read.is_read1 else read.query_name + '_r2'
        if (not read.mate_is_unmapped) and (not read.is_secondary):
            read_pos_dict[read_name] = (read.reference_start, read.reference_end)
            # update pair dict
            pair_read_dict.setdefault(read.query_name, set())
            pair_read_dict[read.query_name].add(ind)

        if (not read.is_proper_pair) or read.is_secondary:
            bad_read_number += 1
            continue
        if not pattern.match(read.cigarstring):
            bad_read_number += 1
            logger.info("Cigar: {}显示当前read: {}存在indel，该read及其mate将不参与模拟".format(read.cigarstring, read_name))
            bad_pair_set.add(read.query_name)
            continue
        if not seq[start:].startswith(ref):
            bad_read_number += 1
            bad_pair_set.add(read.query_name)
            if ref.startswith(seq[start:]):
                logger.info("当前read的末端seq[{}:]刚好仅为ref的一部分,该read及其mate将不参与模拟".format(start))
                logger.info('   即{}中出现{}:{},ref序列{}包含read末端{},意味着ref序列未测完整'.format(
                    read_name, read.reference_name, ref_start, ref, seq[start:start + len(ref)]
                ))
            else:
                logger.info('{}被修改前已经和ref冲突{}:{},{}>{},而将模拟的突变为:{}>{},该read及其mate将不参与模拟'.format(
                    read_name, read.reference_name, ref_start, ref, seq[start:start + len(ref)], ref, alt
                ))
            continue
        candidates.add(ind)
    bam_obj.close()
    logger.info("coverage is {}".format(coverage))

    # remover bad pair from candidates
    for each in bad_pair_set:
        candidates -= pair_read_dict[each]
    logger.info("bad read ratio is {}/{}={}".format(bad_read_number, coverage, bad_read_number/float(coverage)))

    # classify read
    pair_dict = dict()
    single_dict = dict()
    for name, read_set in pair_read_dict.items():
        if len(read_set) == 2:
            pair_key = str(read_pos_dict[name+'_r1'][0]) + '_' + str(read_pos_dict[name+'_r2'][1])
            pair_dict.setdefault(pair_key, set())
            pair_dict[pair_key].add(tuple(read_set))
        else:
            single_key = read_pos_dict[name+'_r1'] if name+'_r1' in read_pos_dict else read_pos_dict[name+'_r2']
            single_dict.setdefault(single_key, set())
            single_dict[single_key].update(read_set)
    dup_pair_dict = {x: y for x, y in pair_dict.items() if len(y) >= 2}
    if dup_pair_dict:
        dup_pair_total = sum(len(dup_pair_dict[x])*2 for x in dup_pair_dict)
        logger.info("mean dup-pair-group size is {}".format(float(dup_pair_total)/len(dup_pair_dict)))
        logger.info("total {}/{}={} dup-pair reads".format(dup_pair_total, coverage, dup_pair_total/float(coverage)))
    unique_pair_dict = {x: y for x, y in pair_dict.items() if len(y) == 1}
    if unique_pair_dict:
        unique_pair_total = sum(len(unique_pair_dict[x])*2 for x in unique_pair_dict)
        logger.info("total {}/{}={} unique-pair reads".format(unique_pair_total, coverage,
                                                              unique_pair_total/float(coverage)))
    unique_single_dict = {x: y for x, y in single_dict.items() if len(y) == 1}
    if unique_single_dict:
        unique_single_total = sum(len(unique_single_dict[x]) for x in unique_single_dict)
        logger.info("下面`likely`指仅用当前read的比对起止位置作为标记来判断是否duplicated ")
        logger.info("total {}/{}={} likely unique-single reads".format(unique_single_total, coverage,
                                                                       unique_single_total/float(coverage)))
    dup_single_dict = {x: y for x, y in single_dict.items() if len(y) >= 2}
    if dup_single_dict:
        dup_single_total = sum(len(dup_single_dict[x]) for x in dup_single_dict)
        logger.info("mean dup-single-group size is {}".format(float(dup_single_total) / len(dup_single_dict)))
        logger.info("total {}/{}={} likely dup-single reads".format(dup_single_total, coverage, dup_single_total/float(coverage)))

    # coverage filtering
    target_changed_read_num = int(round(coverage * expected_af))
    for each in dup_pair_dict:
        # 如果某组dup read 数量已经达到了目标模拟reads数的85%， 那么不应该把这组dup read作为模拟候选
        if len(dup_pair_dict[each])*2 > target_changed_read_num*0.85:
            logger.info("At {}:{}, dup_pair_group {} has {} reads. But target simulation number is {}."
                        "They will be excluded from candidates".format(
                chr_, ref_start, each, len(dup_pair_dict[each])*2, target_changed_read_num)
            )
            for reads in dup_pair_dict[each]:
                candidates -= set(reads)

    for each in dup_single_dict:
        # 如果某组dup read 数量已经达到了目标模拟reads数的85%， 那么不应该把这组dup read作为模拟候选
        if len(dup_single_dict[each]) > target_changed_read_num*0.85:
            logger.info("At {}:{}, dup_single_group {} has {} reads. But target simulation number is {}."
                        "They will be excluded from candidates".format(
                chr_, ref_start, each, len(dup_single_dict[each]), target_changed_read_num)
            )
            candidates -= dup_single_dict[each]

    if coverage < depth_threshold:
        candidates = set()
        logger.warning("{}:{},{}>{} is not simulated for low coverage {}".format(chr_, ref_start, ref, alt, coverage))
    elif len(candidates) < target_changed_read_num - 1:
        candidates = set()
        logger.warning("{}:{},{}>{} is not simulated for not enough valid candidates to change, {}/{}. ".format(
            chr_, ref_start, ref, alt, len(candidates), target_changed_read_num
        ))
    if bad_read_number > coverage*0.1:
        candidates = set()
        logger.warning("{}:{},{}>{} is not simulated for high bad read ratio {}/{}={}".format(
            chr_, ref_start, ref, alt, bad_read_number, coverage, bad_read_number/float(coverage)
        ))

    logger.info("Success to screen {} valid reads from {} ones for simulation".format(len(candidates), coverage))
    return candidates, coverage, dup_pair_dict, dup_single_dict, unique_pair_dict, unique_single_dict


def sieve_candidates(candidates, coverage, dup_pair_dict, dup_single_dict, unique_pair_dict, unique_single_dict, af, logger):
    """
    随机挑选要修改的read过程：
    1. 大体思路：以比对的坐标信息为key，相应的reads为values； 随机抽取这些key，从而随机获得相应需修改的reads信息
    2. 随机的实现： random.sample, 无需设置seed，每次都是随机实验
    3. 假设候选总数为N，随机挑出N-1个，此时获得的是无序列表S
    4. 将S的第一个元素替换为dup-pair中的一条，该条随机挑选获得，为了保证至少选中一个dup-pair。
    5. 总共需要挑选的read数量为：int(round(coverage * af))
    5. 对S进行循环取值，并累加统计挑选出来的reads总数, 且设每轮循环时还需挑选的read数量为Diff。
    6. 如果当前循环取出的read数目超过Diff+3，那么放弃本轮循环挑选出来的reads，目的是为了控制最后的AF偏差
    7. 循环终止的条件：最后挑选出来的候选read数目大于等于目标数量
    """
    target_changed_read_num = int(round(coverage * af))
    all_locations = dup_pair_dict.keys() + dup_single_dict.keys() + unique_pair_dict.keys() + unique_single_dict.keys()
    sampling = random.sample(all_locations, len(all_locations)-1)
    if len(dup_pair_dict.keys()) >= 1:
        sampling[0] = random.sample(dup_pair_dict.keys(), 1)[0]
    final_candidates = set()
    selected_dup_pair = set()
    selected_dup_single =set()
    selected_unique_pair = set()
    selected_unique_single = set()
    for each in sampling:
        if each in dup_pair_dict:
            add_num = len(dup_pair_dict[each])*2
            diff = target_changed_read_num - len(final_candidates)
            if add_num < diff + 3:
                for reads in dup_pair_dict[each]:
                    final_candidates.update(set(reads))
                    selected_dup_pair.update(set(reads))
        elif each in dup_single_dict:
            add_num = len(dup_single_dict[each])
            diff = target_changed_read_num - len(final_candidates)
            if add_num < diff + 3:
                final_candidates.update(dup_single_dict[each])
                selected_dup_single.update(dup_single_dict[each])
        elif each in unique_pair_dict:
            add_num = len(unique_pair_dict[each])*2
            diff = target_changed_read_num - len(final_candidates)
            if add_num < diff + 3:
                for reads in unique_pair_dict[each]:
                    final_candidates.update(set(reads))
                    selected_unique_pair.update(set(reads))
        elif each in unique_single_dict:
            add_num = len(unique_single_dict[each])
            diff = target_changed_read_num - len(final_candidates)
            if add_num < diff + 3:
                final_candidates.update(unique_single_dict[each])
                selected_unique_single.update(unique_single_dict[each])
        final_candidates &= candidates
        over_select = len(final_candidates) - target_changed_read_num
        if over_select >= 0:
            break
    assert abs(len(final_candidates)-target_changed_read_num) <= 3
    logger.info("{} reads classified as `dup-paired` were selected for simulation".format(
        len(selected_dup_pair & candidates)))
    logger.info("{} reads classified as `likely dup-single` were selected for simulation".format(
        len(selected_dup_single & candidates)))
    logger.info("{} reads classified as `unique-paired` were selected for simulation".format(
        len(selected_unique_pair & candidates)))
    logger.info("{} reads classified as `likely unique-single` were selected for simulation".format(
        len(selected_unique_single & candidates)))
    return final_candidates


def change_read(read, ref, alt, start, max_read_len=151):
    """
    修改reads的过程：
    1. 用alt替换ref得到新的read
    2. 如果新的read长度超过max_read_len, 则需要从左端或右端去除部分原来的序列
    """
    qual, seq, read_name = read.qual, read.seq, read.query_name
    read.seq = replace(seq, ref, alt, start=start)
    old_qual = qual[start:start + len(ref)]
    new_qual = ''.join([qual[start]] * len(alt))
    read.qual = replace(qual, old_qual, new_qual, start=start)
    if len(read.seq) > max_read_len:  # 插入突变导致后read长度不能太长
        diff_len = len(alt) - len(ref)
        left_ori_num = start - 1
        right_ori_num = len(read.seq) - start - len(alt) + 1
        seq, qual = read.seq, read.qual
        if right_ori_num >= left_ori_num and right_ori_num > diff_len:
            # 去除右端若干碱基以保证长度
            read.seq = seq[:max_read_len]
            read.qual = qual[:max_read_len]
        elif left_ori_num >= right_ori_num and left_ori_num > diff_len:
            # 去除左端若干碱基以保证长度
            read.seq = seq[diff_len:]
            read.qual = qual[diff_len:]
        else:
            read.seq = seq[start:max_read_len + start]
            read.qual = qual[start:max_read_len + start]
    return read


def simulate_reads(bam_file, mutation_file, af=0.1, depth_threshold=300, expect_mutation_number=100,
                   max_read_len=151, output_dir=os.getcwd()):
    logger = set_logger(output_dir+'/simulate_mutation.log', 'simulate_mutation')
    read_change_logger = set_logger(output_dir+"/track_read_mutation.log")
    mutation_reads = dict()
    mutation_num = 0
    # step one: get mutated reads dict
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(output_dir+'/Introduced_mutation.vcf', 'w') as fw:
        for chr_, pos, ref, alt, af in target_mutation_info(mutation_file, af=af):
            locate_start = pos - 1
            locate_end = locate_start + len(ref)
            candidates_ind, coverage, dup_pair_dict, dup_single_dict, unique_pair_dict, unique_single_dict = \
                get_candidate_info(bam_file, chr_, pos, ref, alt,
                                   depth_threshold=depth_threshold,
                                   expected_af=af, logger=logger
            )
            if not candidates_ind:
                continue
            target_changed_read_num = int(round(coverage * af))
            logger.info('Begin to simulate mutation {}:{},{}>{}.'.format(chr_, pos, ref, alt))
            final_candidates = sieve_candidates(candidates_ind, coverage,
                                                dup_pair_dict, dup_single_dict,
                                                unique_pair_dict, unique_single_dict,
                                                af, logger)

            # mutate read
            changed_read_num = 0
            fetch_obj = bam.fetch(chr_, locate_start, locate_end)
            for ind, read in enumerate(fetch_obj):
                if ind in final_candidates:
                    start = pos - 1 - read.pos
                    qual, seq, read_name = read.qual, read.seq, read.query_name
                    read_name += '_r1' if read.is_read1 else '_r2'
                    # read_name += '_secondary' if read.is_secondary else ''
                    tmp_mark = ''.join([' ']*start + ['|']*len(ref))
                    read_change_logger.info("{}: {}:{},{}>{} start={}".format(read_name, chr_, pos, ref, alt, start))
                    read = change_read(read, ref, alt, start, max_read_len)
                    read_change_logger.info("{}\n{}\n{}".format(seq, tmp_mark, read.seq))
                    mutation_reads[read_name] = (read.seq, read.qual, read.is_reverse)
                    changed_read_num += 1
            # print(changed_read_num, target_changed_read_num, af)
            assert abs(changed_read_num - target_changed_read_num) <= 3
            logger.info('Success to introduce {}/{}={} mutated reads for mutation {}:{},{}>{}.'.format(
                changed_read_num, coverage, changed_read_num/float(coverage), chr_, pos, ref, alt,
            ))
            # write out introduced mutation information
            fw.write("{}\t{}\t.\t{}\t{}\t{}\n".format(chr_, pos, ref, alt, changed_read_num/float(coverage)))
            mutation_num += 1
            if mutation_num == expect_mutation_number:
                break
        # end of loop
        if mutation_num < expect_mutation_number:
            logger.warning("Cannot introduce {} mutations, but only {} ones".format(expect_mutation_number,
                                                                                    mutation_num))
        logger.info("success to introduced {} mutation".format(mutation_num))
        logger.info("total need modified reads number: {}".format(len(mutation_reads)))
    return mutation_reads


def generate_bam_fastq(bam_file, mutation_reads, out_fastq_prefix="simulated"):
    """基于修改的bam排序后生成模拟的fastq，相对耗时，如果没有原fastq文件时使用"""
    # step two: generate modified bam file
    out_bam = bam_file[:-3] + 'faked.bam'
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with pysam.AlignmentFile(out_bam, "wb", template=bam) as new_bam:
            for read in bam.fetch():
                if read.is_secondary:  # bam转换fastq时，碰到secondary reads会发生警告,fastq也不允许一条序列出现在多个位置
                    continue
                read_name = read.query_name
                if read.is_read1:
                    read_name += '_r1'
                else:
                    read_name += '_r2'
                if read_name in mutation_reads:
                    seq, qual, reverse = mutation_reads.pop(read_name)
                    read.seq = seq
                    read.qual = qual
                    read.cigarstring = '{}M'.format(len(read.seq))  # 强制修改，不论正确性
                new_bam.write(read)
    # step three: convert bam to fastq
    sort_bam_name = out_bam + '.readname.sorted.bam'
    os.system('samtools sort -n -@ 4 -o {out} {bam}'.format(out=sort_bam_name, bam=out_bam))
    os.system('bedtools bamtofastq -fq {pre}_R1_001.fastq -fq2 {pre}_R2_001.fastq -i {bam} '.format(
        pre=out_fastq_prefix, bam=sort_bam_name))


def reverse_complement(seq):
    """
    :param seq: 输入序列
    :return: 返回reversed的序列
    """
    seq = seq.upper()
    complement = dict(zip(list("ATCG"), list("TAGC")))
    return ''.join(complement[base] if base in complement else base for base in seq[::-1])


def generate_fastq(mutation_reads, fq1=None, fq2=None, out_fastq_prefix="simulated", compress_out=False):
    """
    生成fastq：
    通读原fastq文件，根据read name找到要修改的read，并用新read替换，替换前，如果发现是reverse比对结果，需进行序列反转后再替换。
    """
    print("total need modified reads number: {}".format(len(mutation_reads)))
    fq_list = list()
    if fq1:
        fq1_obj = gzip.open(fq1, 'rb') if fq1.endswith('.gz') else open(fq1, 'r')
        fq1_out = out_fastq_prefix+'_R1_001.fastq'
        # fq1_out_obj = gzip.open(fq1_out + '.gz', 'wb') if compress_out else open(fq1_out, 'w')
        fq1_out_obj = open(fq1_out, 'w')
        fq_list.append((fq1_obj, fq1_out_obj, '_r1'))
    if fq2:
        fq2_obj = gzip.open(fq2, 'rb') if fq2.endswith('.gz') else open(fq2, 'r')
        fq2_out = out_fastq_prefix+'_R2_001.fastq'
        # fq2_out_obj = gzip.open(fq2_out + '.gz', 'wb') if compress_out else open(fq2_out, 'w')
        fq2_out_obj = open(fq2_out, 'w')
        fq_list.append((fq2_obj, fq2_out_obj, '_r2'))
    for fq_obj, fq_out_obj, read_type in fq_list:
        m = 0
        while 1:
            line = fq_obj.readline()
            if not line:
                break
            if line.startswith('@'):
                read_info = mutation_reads.get(line[1:].split()[0] + read_type)
                if read_info:
                    m += 1
                    fq_out_obj.write(line)
                    seq, qual, reverse = read_info
                    fq_obj.readline()
                    if not reverse:
                        fq_out_obj.write(seq+'\n')
                    else:
                        fq_out_obj.write(reverse_complement(seq) + '\n')
                    fq_out_obj.write(fq_obj.readline())
                    fq_obj.readline()
                    fq_out_obj.write(qual+'\n')
                else:
                    fq_out_obj.write(line)
                    fq_out_obj.write(fq_obj.readline())
                    fq_out_obj.write(fq_obj.readline())
                    fq_out_obj.write(fq_obj.readline())
        fq_obj.close()
        fq_out_obj.close()
        print("read{} modified number:{}".format(read_type, m))
    if fq1 and compress_out:
        os.system('gzip {} &'.format(fq1_out))
    if fq2 and compress_out:
        os.system('gzip {} '.format(fq2_out))


def simulate_pipeline(bam_file, mutation_file, af=0.1, depth_threshold=800, expect_mutation_number=100, output_dir=None,
                      max_read_len=151, fq1=None, fq2=None, out_fastq_prefix='simulated', compress_out=True):
    """
    read的分类说明：
    # dup-pair-read: read 和 read' 对应的r1和r2的比对起始位置相同,且他们的mate的比对终止位置也相同， 模拟时read和read'需同时修改
    # pair-read： 若r1 and r2 都包含突变位点则属于pair, 模拟时r1和r2需同时修改
    # unique-pair： pair-read 减去 dup-pair-read
    # single：r1和r2中只有一条read比对当前要模拟的区域，由于此时若要同时获得两条read的比对位置，势必要通读bam文件，耗时较长。
    因此下面用likely的策略来对read进行分类，以供后续突变模拟的选择使用。
    # likely unique-single： likely指仅用当前read的比对起止位置作为标记来判断是否duplicated，
    # likely dup-single： likely指仅用当前read的比对起止位置作为标记来判断是否duplicated，模拟时需同时修改
    # bad-read: (not read.is_proper_pair) | read.is_secondary | cigar-string不全为M | 当前突变位点处已经存在突变

    候选read的筛选条件：
    1. 突变起始位置要落在read比对上的区间
    2. read能比对上，且cigar-string不能为空, 且是is-proper—pair，且不能是is-secondary
    3. read1和read2在当前模拟的突变位点处不能存在突变；cigar-string只能全是M
    4. 如read末端仅包含部分ref（突变前的参考序列）序列，则该read及其mate也不能参与模拟，且有日志记录
    5. 如果某组dup read 数量已经达到了目标模拟reads数的85%， 那么不应该把这组dup read作为模拟候选，且有日志记录
    6. bad-read-number的比例占coverage的0.1以上时，不对当前提供的位点进行突变模拟，且日志有记录
    7. coverage低于提供的阈值depth_threshold值时，不对当前提供的位点进行突变模拟，且日志有记录

    随机挑选要修改的read过程：
    1. 大体思路：以比对的坐标信息为key，相应的reads为values； 随机抽取这些key，从而随机获得相应需修改的reads信息
    2. 随机的实现：使用 random.sample, 无seed设置，因为每次都是随机实验
    3. 假设候选总数为N，随机挑出N-1个，此时获得的是无序列表S
    4. 将S的第一个元素替换为dup-pair中的一条，该条随机挑选获得，目的是为了保证至少选中一个dup-pair。
    5. 总共需要挑选的read数量为：int(round(coverage * af))
    5. 对S进行循环取值，并累加统计挑选出来的reads总数, 且设每轮循环时还需挑选的read数量为Diff。
    6. 如果当前循环取出的read数目超过Diff+3，那么放弃本轮循环挑选出来的reads，目的是为了控制最后的AF偏差
    7. 循环终止的条件：最后挑选出来的候选read数目大于等于目标数量

    修改reads的过程：
    1. 对所有挑选出来的read，用alt替换ref得到新的read
    2. 如果新的read长度超过max_read_len, 则需要从左端或右端去除部分原来的序列

    生成新fastq文件，两种思路：
    1. 通读原fastq文件，根据read name找到要修改的read，并用新read替换，替换前，如果发现是reverse比对结果，需进行序列反转后再替换。
    2.  直接修改原bam文件，生成新的bam文件，然后用samtools排序，最后用bedtools生成最fastq

    程序功能：基于bam文件模拟snv或indel的突变，支持snv和indel的同时模拟，indel长度不能比read还长，最终生成fastq文件

    :param bam_file: bam文件路径
    :param mutation_file: 目标模拟的突变信息文件，vcf格式，第6列为可选项，为每个位点提供af信息
    :param af: 指定模拟的突变频率，如已经在mutation_file的第6列提供，可以指定为0.
    :param depth_threshold: 突变位点最低覆盖度，低于该覆盖度的位点将不被模拟，因此mutation_file需提供足够多的模拟位点.
    :param expect_mutation_number: 期望从mutation_file中提取的模拟位点个数，多余的信息不被考虑.
    :param output_dir: fastq，log等的输出目录
    :param max_read_len: 由于插入突变会导致模拟的read增加长度，需限制最大read长度，也作为最大模拟indel的限制
    :param fq1: 原read1的fastq文件，可为压缩文件，默认：None， 为None时利用bam文件提取fastq，耗时较长
    :param fq2: 原read2的fastq文件，可为压缩文件，默认：None， 为None时利用bam文件提取fastq，耗时较长
    :param out_fastq_prefix: 输出的fastq文件名前缀，不可带路径信息
    :param compress_out: 是否输出压缩的fastq文件
    :return: 除了生成fastq文件外，还有4个记录文件，分别记录参数信息，引入的突变信息*vcf，模拟过程的log，模拟结果的log
    """
    arg_dict = dict(locals())
    if output_dir is None:
        output_dir = os.getcwd()
    if not os.path.exists(output_dir):
            os.system("mkdir -p {} ".format(output_dir))
    with open(output_dir+"/Argument_detail.json", 'w') as f:
        json.dump(arg_dict, f, indent=2, sort_keys=True)
    if af == 0:
        af = None
    new_read_dict = simulate_reads(bam_file, mutation_file, af=af, depth_threshold=depth_threshold,
                                   output_dir=output_dir, expect_mutation_number=expect_mutation_number,
                                   max_read_len=max_read_len)
    out_fastq_prefix = os.path.join(output_dir, out_fastq_prefix)
    if fq1 or fq2:
        generate_fastq(new_read_dict, fq1=fq1, fq2=fq2, out_fastq_prefix=out_fastq_prefix, compress_out=compress_out)
    else:
        generate_bam_fastq(bam_file, new_read_dict, out_fastq_prefix=out_fastq_prefix)


def introduce_command(func):
    import argparse
    import inspect
    import json
    parser = argparse.ArgumentParser(description=func.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    func_args = inspect.getargspec(func)
    arg_names = func_args.args
    arg_defaults = func_args.defaults
    arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
    for arg, value in zip(arg_names, arg_defaults):
        if value == 'None':
            parser.add_argument('-'+arg, required=True, metavar=arg)
        elif type(value) == bool:
            if value:
                parser.add_argument('--'+arg, action="store_false", help='default: True')
            else:
                parser.add_argument('--'+arg, action="store_true", help='default: False')
        elif value is None:
            parser.add_argument('-' + arg, default=value, metavar='Default:' + str(value), )
        else:
            parser.add_argument('-' + arg, default=value, type=type(value), metavar='Default:' + str(value), )
    if func_args.varargs is not None:
        print("warning: *varargs is not supported, and will be neglected! ")
    if func_args.keywords is not None:
        print("warning: **keywords args is not supported, and will be neglected! ")
    args = parser.parse_args().__dict__
    with open("Argument_detail_for_{}.json".format(os.path.basename(__file__).split(".")[0]), 'w') as f:
        json.dump(args, f, indent=2, sort_keys=True)
    start = time.time()
    func(**args)
    print("total time: {}s".format(time.time() - start))


if __name__ == '__main__':
    introduce_command(simulate_pipeline)

