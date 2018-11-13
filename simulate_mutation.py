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
__author__ = 'gdq'
"""
bam是0-based坐标，而pysam用0-based方法读取信息
0. fetch命令可以读取某个区域比对上的所有reads，随后可以决定要随机修改的reads集合
1. 用fetch读取bam文件，结合输入的突变位点信息，获得需要修改的reads集合
2. 按顺序读取bam文件，并一行行写入新的文件当中去，在写之前进行判断并修改，需要过滤掉secondary-alignment，后直接获得修改后的bam文件
3. 把bam文件转换为fastq, 这中间需要对bam进行排序
4. 弊端：暂时不能模拟大片段的indel，即indel超过read长度
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


def get_candidate_info(bam_obj, chr_, ref_start, ref, alt, depth_threshold=300, expected_af=0.1, logger=None):
    # dup-read: read 和 read' 对应的r1和r2的比对起始位置相同,且他们的mate的比对起始位置也相同， 模拟时read和read'需同时修改
    # pair-read： r1 and r2 都包含突变位点属于pair, 模拟时r1和r2需同时修改
    if logger is None:
        logger = set_logger('get_candidates.log', 'get_candidates')
    pattern = re.compile('\d+M$')
    locate_start = ref_start - 1
    locate_end = locate_start + len(ref)
    dup_read_dict = dict()
    pair_read_dict = dict()
    coverage = 0
    candidates = set()
    logger.info("Start to screen valid reads for mutation simulation")
    bad_pair_set = set()  # 存那些至少有一端read存在变异的pair，他们需要从candidates中被排除
    for ind, read in enumerate(bam_obj.fetch(chr_, locate_start, locate_end)):
        if locate_start not in read.get_reference_positions():
            continue
        if read.cigarstring is None or read.is_unmapped:
            continue
        coverage += 1
        seq = read.seq
        start = ref_start - 1 - read.pos  # ref序列在当前read上的起始位置
        read_name = read.query_name + '_r1' if read.is_read1 else read.query_name + '_r2'

        # update dup dict
        if not read.mate_is_unmapped:
            # mate = bam_obj.mate(read)
            # dup_key = [read.reference_start, read.reference_end, read.next_reference_start, mate.reference_end]
            # if read.is_duplicate is True:
            dup_key = (read.reference_start, read.reference_end, read.next_reference_start, read.cigarstring)
            dup_read_dict.setdefault(dup_key, set())
            dup_read_dict[dup_key].add(ind)
            # update pair dict
            pair_read_dict.setdefault(read.query_name, set())
            pair_read_dict[read.query_name].add(ind)

        # other filtering
        if not read.is_proper_pair or read.is_secondary:
            continue
        if not pattern.match(read.cigarstring):
            logger.info("Cigar: {}显示当前read: {}存在indel，该read及其mate将不参与模拟".format(read.cigarstring, read_name))
            bad_pair_set.add(read.query_name)
            continue
        if not seq[start:].startswith(ref):
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
    # remover bad pair from candidates
    for each in bad_pair_set:
        candidates -= pair_read_dict[each]

    # classify reads
    dup_read_dict = {x: y for x, y in dup_read_dict.items() if len(y) >= 2}
    pair_read_dict = {x: y for x, y in pair_read_dict.items() if len(y) >= 2}

    # coverage filtering
    target_changed_read_num = int(round(coverage * expected_af))
    for each in dup_read_dict:
        # 如果某组dup read 数量已经达到了目标模拟reads数的85%， 那么不应该把这组dup read作为模拟候选
        if len(dup_read_dict[each]) > target_changed_read_num*0.85:
            logger.info("At {}:{},dup_group {} has {} reads. But target simulate number is {} ->Don't change".format(
                chr_, ref_start, each, len(dup_read_dict[each]), target_changed_read_num))
            for ind in dup_read_dict[each]:
                if ind in candidates:
                    candidates.remove(ind)
    if not coverage >= depth_threshold:
        candidates = set()
        logger.warning("{}:{},{}>{} is not simulated for low coverage {}".format(chr_, ref_start, ref, alt, coverage))
    if len(candidates) < target_changed_read_num - 1:
        candidates = set()
        logger.warning("{}:{},{}>{} is not simulated for not enough valid candidates to change, {}/{}. ".format(
            chr_, ref_start, ref, alt, len(candidates), target_changed_read_num
        ))
    logger.info("Success to screen {} valid reads from {} ones for simulation".format(len(candidates), coverage))
    return candidates, coverage, dup_read_dict, pair_read_dict


def sieve_candidates(candidates, coverage, dup_read_dict, pair_read_dict, af, logger):
    """
    dup-read: read 和 read' 对应的r1和r2的比对起始位置相同,且他们的mate的比对起始位置也相同， 模拟时read和read'需同时修改
    pair-read： r1 and r2 都包含突变位点属于pair, 模拟时r1和r2需同时修改
    先根据dup—read-dict和pair-read-dict和candidates的交集获得可以修改的dup-pair-read集合
    根据（dup-pair-read/coverage）*（coverage*af），计算出模拟时应该修改多少条dup-pair-read
    抽样dup-pair-read集合，获得要模拟的dup-pair-read，如果数量足够（过多时，随机删除多余的），抽样结束；
    否则从剩下的candidates（排除所有dup-pair-read）继续抽样，直到抽到足够数量的，注意每抽到一个pair-read，则pair也要一起算入。
    """
    duplicate_group_size_list = []
    target_read_num = int(round(af * coverage))

    # get reverse dict of dup read dict
    reverse_dup_read_dict = dict()
    for key, value in dup_read_dict.items():
        duplicate_group_size_list.append(len(value))
        for each in value:
            reverse_dup_read_dict[each] = key
    logger.info("coverage is {}".format(coverage))
    logger.info('candidates number is {}'.format(len(candidates)))
    mean_dup_num = float(sum(duplicate_group_size_list)) / len(duplicate_group_size_list)
    logger.info("mean dup group size is {}".format(mean_dup_num))
    logger.info("total {}/{}={} dup reads".format(len(reverse_dup_read_dict), coverage,
                                                  len(reverse_dup_read_dict)/float(coverage)))
    # get reverse dict of pair read dict
    reverse_pair_read_dict = dict()
    for key, value in pair_read_dict.items():
        for each in value:
            reverse_pair_read_dict[each] = key
    logger.info("total {}/{}={} pair reads".format(len(reverse_pair_read_dict), coverage,
                                                   len(reverse_pair_read_dict)/float(coverage)))

    # select dup pair read for simulation
    dup_pair_read_set = set()
    for k, v in pair_read_dict.items():
        tmp = list(v)
        if tmp[0] in reverse_dup_read_dict and tmp[1] in reverse_dup_read_dict:
            dup_pair_read_set.update(v)
    # print('{} strict dup pair read'.format(len(dup_pair_read_set)))
    loose_dup_pair_read_set = set(reverse_dup_read_dict.keys()) & set(reverse_pair_read_dict.keys())
    logger.info("total {}/{}={} dup-pair reads".format(len(dup_pair_read_set), coverage,
                                                       len(dup_pair_read_set)/float(coverage)))
    candidates_num = float(len(candidates))
    dup_pair_candidates = dup_pair_read_set & candidates
    logger.info("total {} dup-pair valid candidates".format(len(dup_pair_candidates)))
    target_dup_pair_num = int(af*len(dup_pair_read_set))
    if target_dup_pair_num == 0 and len(dup_pair_candidates) >= 1:
        target_dup_pair_num = 1
    logger.info("target modified dup-pair reads number is {}".format(target_dup_pair_num))
    if target_dup_pair_num > len(dup_pair_candidates):
        target_dup_pair_num = len(dup_pair_candidates)
    random_select_dup_pair = random.sample(dup_pair_candidates, len(dup_pair_candidates)-1)
    final_select_dup_pair = set()
    for i, each in enumerate(random_select_dup_pair):
        # if len(pair_read_dict[reverse_pair_read_dict[each]] & candidates) % 2 != 0:
        #     print(reverse_pair_read_dict[each], ', its mate is filtered')
        final_select_dup_pair.update(pair_read_dict[reverse_pair_read_dict[each]] & candidates)
        final_select_dup_pair.update(dup_read_dict[reverse_dup_read_dict[each]] & candidates)
        if len(final_select_dup_pair & dup_pair_candidates) >= target_dup_pair_num + 4:
            final_select_dup_pair -= (pair_read_dict[reverse_pair_read_dict[each]] & candidates)
            final_select_dup_pair -= (dup_read_dict[reverse_dup_read_dict[each]] & candidates)
        elif len(final_select_dup_pair & dup_pair_candidates) >= target_dup_pair_num-1:
            logger.info('Get enough dup-pair candidates in advance !')
            break
        if len(final_select_dup_pair) >= target_read_num:
            logger.info('Get enough candidates while selecting dup-pair-reads !')
            break
    all_dup_singles = set(reverse_dup_read_dict.keys()) - loose_dup_pair_read_set
    selected_dup_singles = final_select_dup_pair - dup_pair_candidates
    target_dup_single_num = int(len(all_dup_singles)*af*0.9)
    still_need_dup_single_num = target_dup_single_num - len(selected_dup_singles)
    if still_need_dup_single_num > 0:
        tmp_candidates = (all_dup_singles - selected_dup_singles) & candidates
        if still_need_dup_single_num > len(tmp_candidates):
            still_need_dup_single_num = len(tmp_candidates)
        tmp_set = random.sample(tmp_candidates, still_need_dup_single_num)
        final_select_dup_pair.update(set(tmp_set))
    else:
        logger.warning("Have selected more dup single reads than expected: {} vs {}".format(
            len(selected_dup_singles), target_dup_single_num))

    # select other reads for simulation
    remained_num = target_read_num - len(final_select_dup_pair)
    if remained_num <= -2:
        logger.info("remove {} over-selected candidates".format(abs(remained_num)))
        for i in range(abs(remained_num)):
            final_select_dup_pair.remove(random.sample(final_select_dup_pair, 1)[0])  # 去掉多选的
    other_select_reads = set()
    if remained_num > 0:
        # new_candidates = candidates-final_select_dup_pair-loose_dup_pair_read_set-all_dup_singles
        new_candidates = candidates-set(reverse_dup_read_dict.keys())
        # print(len(new_candidates), remained_num)
        if len(new_candidates) < remained_num:
            logger.warning("Fewer candidates than needed!")
            remained_num = len(new_candidates)
        random_reads = random.sample(new_candidates, remained_num)
        for each in random_reads:
            if each in reverse_pair_read_dict:
                other_select_reads.update(pair_read_dict[reverse_pair_read_dict[each]] & candidates)
            elif each in reverse_dup_read_dict:
                other_select_reads.update(dup_read_dict[reverse_dup_read_dict[each]] & candidates)
            else:
                other_select_reads.add(each)
            if len(other_select_reads) >= remained_num:
                over_selected = len(other_select_reads) - remained_num
                if over_selected >= 2:
                    logger.info("remove {} over-selected read".format(over_selected))
                    for i in range(over_selected):
                        other_select_reads.remove(random.sample(other_select_reads, 1)[0])  # 去掉多选的
                break
    selected_candidates = final_select_dup_pair | other_select_reads
    # print(len(selected_candidates), target_read_num, af)
    assert abs(len(selected_candidates) - target_read_num) <= 3
    dup_num, pair_num, dup_pair_num, other = 0, 0, 0, 0
    dup_pair_group_num_of_pair = set()
    dup_pair_group_num_of_dup = set()
    only_pair_group_num = set()
    only_dup_group_num = set()
    for each in selected_candidates:
        if each in dup_pair_read_set:
            dup_pair_group_num_of_pair.add(reverse_pair_read_dict[each])
            dup_pair_group_num_of_dup.add(reverse_dup_read_dict[each])
            dup_pair_num += 1
        elif each in reverse_pair_read_dict:
            only_pair_group_num.add(reverse_pair_read_dict[each])
            pair_num += 1
        elif each in reverse_dup_read_dict:
            only_dup_group_num.add(reverse_dup_read_dict[each])
            dup_num += 1
        else:
            other += 1
    logger.info("{}({} groups) reads classified as `dup-read-notPaired` were selected for simulation".format(
        dup_num, len(only_dup_group_num)))
    logger.info("{}({} groups) reads classified as `pair-read-notDuplicated` were selected for simulation".format(
        pair_num, len(only_pair_group_num)))
    logger.info("{}({} dup groups,{} pair groups) reads classified as `dup-pair-read` were selected for simulation".format(
        dup_pair_num, len(dup_pair_group_num_of_dup), len(dup_pair_group_num_of_pair)))
    logger.info("{} reads classified as `other` were selected for simulation".format(other))
    # return selected_candidates, (dup_num, pair_num, dup_pair_num, other)
    return selected_candidates


def change_read(read, ref, alt, start, max_read_len=151):
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
            # 先去除左端若干，再去除右端若干
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
            candidates_ind, coverage, dup_read_dict, pair_read_dict = get_candidate_info(
                bam, chr_, pos, ref, alt, depth_threshold=depth_threshold, expected_af=af, logger=logger
            )
            if not candidates_ind:
                continue
            target_changed_read_num = int(round(coverage * af))
            logger.info('Begin to simulate mutation {}:{},{}>{}.'.format(chr_, pos, ref, alt))
            final_candidates = sieve_candidates(candidates_ind, coverage, dup_read_dict, pair_read_dict, af, logger)

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
    seq = seq.upper()
    complement = dict(zip(list("ATCG"), list("TAGC")))
    return ''.join(complement[base] if base in complement else base for base in seq[::-1])


def generate_fastq(mutation_reads, fq1=None, fq2=None, out_fastq_prefix="simulated", compress_out=False):
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


def simulate_pipeline(bam_file, mutation_file, af=0.1, depth_threshold=300, expect_mutation_number=100, output_dir=None,
                      max_read_len=151, fq1=None, fq2=None, out_fastq_prefix='simulated', compress_out=True):
    """
    dup-read: read 和 read' 对应的r1和r2的比对起始位置相同,且他们的mate的比对起始位置也相同， 模拟时read和read'需同时修改
    pair-read： r1 and r2 都包含突变位点, 称之为pair, 模拟时r1和r2需同时修改
    先根据dup—read-dict和pair-read-dict和candidates的交集获得可以修改的dup-pair-read集合
    根据（dup-pair-read/coverage）*（coverage*af），计算出模拟时应该修改多少条dup-pair-read
    抽样dup-pair-read集合，获得要模拟的dup-pair-read，如果数量足够（过多时，随机删除多余的），抽样结束；
    否则从剩下的candidates（排除所有dup-pair-read）继续抽样，直到抽到足够数量的模拟reads，注意每抽到一个pair-read，则pair也要一起算入。

    :param bam_file: bam文件路径
    :param mutation_file: 目标模拟的突变信息文件，vcf格式，第6列为可选项，为每个位点提供af信息
    :param af: 指定模拟的突变频率，如已经在mutation_file的第6列提供，可以指定为None. 默认:0.1
    :param depth_threshold: 突变位点最低覆盖度，低于该覆盖度的位点将不被模拟，因此mutation_file需提供足够多的模拟位点. 默认：300
    :param expect_mutation_number: 期望从mutation_file中提取的模拟位点个数，多余的信息不被考虑. 默认:100
    :param output_dir: fastq，log等的输出目录
    :param max_read_len: 由于插入突变会导致模拟的read增加长度，需限制最大read长度，另外也作为最大模拟indel的限制，默认：151bp
    :param fq1: 原read1的fastq文件，可为压缩文件，默认：None， 为None时利用bam文件提取fastq，耗时较长
    :param fq2: 原read2的fastq文件，可为压缩文件，默认：None， 为None时利用bam文件提取fastq，耗时较长
    :param out_fastq_prefix: 输出的fastq文件名前缀，不可带路径信息，默认：'simulated'
    :param compress_out: 是否输出压缩的fastq文件，默认：True
    :return: 无返回值
    """
    arg_dict = dict(locals())
    if output_dir is None:
        output_dir = os.getcwd()
    if not os.path.exists(output_dir):
            os.system("mkdir -p {} ".format(output_dir))
    with open(output_dir+"/Argument_detail.json", 'w') as f:
        json.dump(arg_dict, f, indent=2, sort_keys=True)
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
    with open("Argument_detail.json", 'w') as f:
        json.dump(args, f, indent=2, sort_keys=True)
    start = time.time()
    func(**args)
    print("total time: {}s".format(time.time() - start))


if __name__ == '__main__':
    introduce_command(simulate_pipeline)

    # import glob
    # start = time.time()
    # bam_file = glob.glob('../*.RAW.bam')[0]
    # mutation_file = 'filtered.vcf_0'
    # sample_name = 'snv0af0002'
    # simulate_pipeline(bam_file, mutation_file, af=0.002, out_fastq_prefix=sample_name)
    # print('simulation total time: {}s'.format(time.time() - start))
    #
    # with open('simulate_project.cfg', 'w') as fw:
    #     fw.write('<project>\n')
    #     fw.write('fastq_dir = ./\n')
    #     fw.write('sample_id = {name}\n'.format(name=sample_name))
    #     fw.write('sample_name = {name}\n'.format(name=sample_name))
    #     fw.write('sample_cfg = LK291-Oncoscreen.520.v3\n')
    #     fw.write('sample_type = FFPE\n')
    #     fw.write('experiment_id = simulate_project\n')
    #     fw.write('</project>\n')
    #
    # cmd = "python ~/../test.rd/app/rdscripts_release/pipeline/run_pipeline.py "
    # cmd += "-p simulate_project.cfg "
    # cmd += "-c ~/../test.rd/app/rdscripts_release/configure_file/cfg.template.rd "
    # cmd += "-o ./"
    # print(cmd)
    # os.system(cmd)
