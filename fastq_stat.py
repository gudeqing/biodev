# coding=utf-8
import os
import gzip
import re
import time
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pandas as pd
import argparse
from collections import Counter
# from concurrent.futures import ProcessPoolExecutor as Pool
from multiprocessing import Pool


def chunk_file(in_file, chunk_size=40000):
    seq_list = []
    i = 0
    with gzip.open(in_file, 'rb') as fr:
        for num, seq in enumerate(fr):
            if num % 4 != 1:
                continue
            seq_list.append(seq)
            i += 1
            if i == chunk_size:
                yield seq_list
                i = 0
                seq_list = []
        else:
            if seq_list:
                yield seq_list


def counter(line_list, exp):
    pattern = re.compile(exp)
    px_count_dict = Counter()
    for line in line_list:
        match = pattern.match(line.upper().strip())
        if match:
            for group_id, nt_str in enumerate(match.groups(), start=1):
                for pos, nt in enumerate(nt_str, start=1):
                    key = str(group_id)+'_'+str(pos)
                    px_count_dict.setdefault(key, Counter())
                    px_count_dict[key].setdefault(nt, 0)
                    px_count_dict[key][nt] += 1
    return px_count_dict


def multi_process(in_file, exp, p_num=3, chunk_size=40000):
    pool = Pool(p_num)
    results = list()
    for line_list in chunk_file(in_file, chunk_size):
        future = pool.apply_async(counter, (line_list, exp))
        results.append(future)
    pool.close()
    pool.join()
    # merge result
    px_count_dict = Counter()
    for future in results:
        # Counter的update方法可以累加统计结果
        result = future.get()
        if result:
            for each in result:
                if each in px_count_dict:
                    px_count_dict[each].update(result[each])
                else:
                    px_count_dict[each] = result[each]
    return px_count_dict


def output_result(px_count_dict, in_file):
    # process result
    count_dict_list = []
    count_ratio_dict_list = []
    order_key = [[x] + x.split('_') for x in px_count_dict.keys()]
    order_key = sorted(order_key, key=lambda x:(int(x[1]), int(x[2])))
    order_key = [x[0] for x in order_key]
    for pos in order_key:
        count_dict = px_count_dict[pos]
        count_dict_list.append(count_dict)
        total = float(sum(count_dict.values()))
        # print(count_dict)
        count_ratio_dict = Counter({x: y/total for x, y in count_dict.items()})
        count_ratio_dict_list.append(count_ratio_dict)
        print("position {} count:\n A={}, T={}, C={}, G={}".format(
            pos, count_dict['A'], count_dict['T'], count_dict['C'], count_dict['G'])
        )
        print("position {} count ratio:\n A={}, T={}, C={}, G={}".format(
            pos, count_ratio_dict['A'], count_ratio_dict['T'],
            count_ratio_dict['C'], count_ratio_dict['G'])
        )
    # plot and save
    p = pd.DataFrame(count_dict_list, index=order_key).fillna(0)
    p.plot(kind='bar', stacked=True)
    plt.xticks(fontsize=7, rotation=90)
    plt.savefig(os.path.basename(in_file)+'.baseCount.pdf')
    plt.close()
    p.to_csv(os.path.basename(in_file) + '.baseCount.txt', index=True, header=True, sep='\t')
    p = pd.DataFrame(count_ratio_dict_list, index=order_key).fillna(0)
    p.plot(kind='bar', stacked=True)
    plt.xticks(fontsize=7, rotation=90)
    plt.savefig(os.path.basename(in_file)+'.baseRatio.pdf')
    plt.close()
    p.to_csv(os.path.basename(in_file) + '.baseRatio.txt', index=True, header=True, sep='\t')
    print('Total time: {}s'.format(time.time() - start))


def get_read_len(seq_list):
    len_dict = Counter()
    for line in seq_list:
        key = len(line.strip())
        len_dict.setdefault(key, 0)
        len_dict[key] += 1
    return len_dict


def get_all_read_len(in_file, p_num=3, chunk_size=50000, plot_rank=None):
    pool = Pool(p_num)
    results = list()
    for line_list in chunk_file(in_file, chunk_size):
        future = pool.apply_async(get_read_len, (line_list, ))
        results.append(future)
    pool.close()
    pool.join()
    len_dict = Counter()
    for future in results:
        len_dict.update(future.get())
    data = pd.Series(len_dict)
    # plot after sort
    name = os.path.basename(in_file)
    data= data.sort_values(ascending=False)
    ratio = data/data.sum()
    out_data = pd.DataFrame(dict(number=data, ratio=ratio))
    out_data.index.name = 'length'
    out_data.to_csv('{}.len.txt'.format(name), header=True, index=True, sep='\t')
    if plot_rank is None:
        sum_ratio = 0
        for rank, i in enumerate(ratio):
            sum_ratio += i
            if sum_ratio >= 0.999:
                plot_rank = rank
                break
    # plot frequency bar
    data = data[:plot_rank+1]
    data = data.sort_index()
    plt.subplot(211)
    data.plot.bar()
    plt.xticks(fontsize=7, rotation=90)
    plt.grid(axis='y', color='gray', linestyle='--')
    # plot ratio bar
    plt.subplot(212)
    ratio = ratio[:plot_rank]
    ratio = ratio.sort_index()
    ratio.plot.bar()
    plt.xticks(fontsize=7, rotation=90)
    plt.xlabel('Read Length')
    plt.grid(axis='y', color='gray', linestyle='--')
    plt.savefig('{}.len.distribution.pdf'.format(name))
    plt.close()


if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subparser",
                                       description="base: 碱基组成比例统计；length: 序列长度分布统计")
    #  base stat
    subparser1 = subparsers.add_parser("base", description="碱基组成比例统计")
    subparser1.add_argument('-f', help='fastq.gz文件路径', required=True, nargs='+')
    subparser1.add_argument('-p', help='process number', type=int, default=5)
    subparser1.add_argument('-e', default=r'^[T]{3,15}([ATCG]{4}).*',
                        help='默认为"^[T]{3,15}([ATCG]{4}).*", '
                             '匹配正则表达式，该脚本将统计所有()号里匹配到的碱基成分.')
    # length stat
    subparser2 = subparsers.add_parser('length', description="序列长度分布统计")
    subparser2.add_argument('-f', help='fastq.gz文件路径', required=True, nargs='+')
    subparser2.add_argument('-p', help='process number', type=int, default=5)
    subparser2.add_argument('-n', help='画出排名前n的分布,默认画累计占比超过0.95的', type=int, default=None)

    #
    args = parser.parse_args()
    if args.subparser == 'base':
        in_file = args.f
        exp = args.e
        p_num = args.p
        for each in in_file:
            output_result(multi_process(each, exp, p_num, chunk_size=30000), each)
    elif args.subparser == 'length':
        in_file = args.f
        p_num = args.p
        for each in in_file:
            get_all_read_len(each, p_num, chunk_size=50000, plot_rank=args.n)
    else:
        print("Wrong sub-command: {}".format(args.subparser))

