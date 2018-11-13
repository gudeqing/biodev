# coding=utf-8
import os
import glob
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pandas as pd
__author__ = 'gdq'


def get_introduced_mutation(out_dir):
    vcf = glob.glob('{}/Introduced_mutation.vcf'.format(out_dir))[0]
    introduced = set()
    with open(vcf) as fr:
        for line in fr:
            chr_, pos, _, ref, alt = line.split('\t')[0:5]
            introduced.add('{}:{},{}>{}\n'.format(chr_, pos, ref, alt))
    return introduced


def get_original_mutation(pipeline_out, vardict=True):
    if vardict:
        report_files = glob.glob('{}/*/report_vardict/*_report.xls'.format(pipeline_out))
    else:
        report_files = glob.glob('{}/*/report/*_report.xls'.format(pipeline_out))
    for each in report_files:
        cmd = "less {} | grep -v '#' | cut -f8,9,13,14 | sed 1d | sed 's/\t/:/;s/\t/,/;s/\t/>/' > {}/report_detect.list"
        cmd = cmd.format(each, os.path.dirname(each))
        os.system(cmd)
    if vardict:
        out_path = "{}/vardict_original_mutation.list".format(pipeline_out)
        cmd = 'cat {}/*/report_vardict/report_detect.list | sort | uniq > {}'.format(pipeline_out, out_path)
    else:
        out_path = "{}/varscan_original_mutation.list".format(pipeline_out)
        cmd = 'cat {}/*/report/report_detect.list | sort | uniq > {}'.format(pipeline_out, out_path)
    print(cmd)
    os.system(cmd)
    return out_path


def summary(out_dir, original_existed_mutation, allowed_distance=100, panel_size=1789799, vardict=True):
    sample_result_dirs = glob.glob("{}/*af*seed*".format(out_dir))
    # print(sample_result_dirs)
    result_dict = dict()
    for sample in sample_result_dirs:
        # print("Stat for: {}".format(sample))
        if vardict:
            report_file = glob.glob('{}/report_vardict/*_report.xls'.format(sample))
        else:
            report_file = glob.glob('{}/report/*_report.xls'.format(sample))
        if not report_file:
            print("Find no *_report.xls in {}".format(sample))
            continue
        else:
            report_file = report_file[0]
        cmd = "less {} | grep -v '#' | cut -f8,9,13,14 | sed 1d | sed 's/\t/:/;s/\t/,/;s/\t/>/' > {}/report_detect.list"
        cmd = cmd.format(report_file, os.path.dirname(report_file))
        os.system(cmd)
        report_list = open('{}/report_detect.list'.format(os.path.dirname(report_file))).readlines()
        origin_list = open(original_existed_mutation).readlines()
        report = set(report_list) - set(origin_list)
        simulate_out_dir = out_dir.split('_runPipeline')[0]
        introduced_list = get_introduced_mutation(simulate_out_dir + '/' + os.path.basename(sample))
        introduced_total_len = 0
        for each in introduced_list:
            ref, alt = each.strip().split(',')[1].split(">")
            if len(ref) > len(alt):
                introduced_total_len += len(ref)-len(alt)
            elif len(ref) == len(alt):
                introduced_total_len += 1

        out1 =  os.path.join(os.path.dirname(report_file), 'report_detect.list')
        success = set(introduced_list) & set(report)
        out2 = os.path.join(os.path.dirname(report_file), 'success_detect.list')
        failed = set(introduced_list) - set(report)
        out3 = os.path.join(os.path.dirname(report_file), 'failed_detect.list')
        more = set(report) - set(introduced_list)
        out4 = os.path.join(os.path.dirname(report_file), 'more_detect.list')
        for out, content in zip([out1, out2, out3, out4], [report, success, failed, more]):
            with open(out, 'w') as fw:
                for each in sorted(content):
                    fw.write(each)

        # rescue some failed
        # 进行并非严格的挽救程序，比如检测到了，位置稍微偏了，而且不考虑突变的内容是否完全与引入的一致，也算是成功检测到
        rescued_fail = set()
        rescued_more = set()
        for each_fail in failed:
            chr_, pos = each_fail.split(',')[0].split(":")
            for every_fail in more:
                chr_2, pos2 = every_fail.split(',')[0].split(":")
                if chr_ == chr_2 and abs(int(pos) - int(pos2)) < allowed_distance:
                    rescued_fail.add(each_fail)
                    rescued_more.add(every_fail)
        if rescued_more:
            print('允许一些位置的偏差后，成功从假阳性集合中挽救回{}个实际为真阳性的突变: {}'.format(len(rescued_more), sample))
        final_more = more - rescued_more
        final_failed = failed - rescued_fail
        final_success = success | rescued_fail

        # find those filtered by polish
        # 进行并非严格的挽救程序，有些引入的突变由于在"repeat区域"或"af稍微低了些"或"failed polish" 而没有报告出来，
        # 不考虑突变的内容是否完全与引入的一致，也算是成功检测到
        if vardict:
            report_raw = glob.glob('{}/report_vardict/*.variations.RAW.xls'.format(sample))
        else:
            report_raw = glob.glob('{}/report/*.variations.RAW.xls'.format(sample))
        failed_polish = set()
        failed_repeat = set()
        af_filtered = set()
        if report_raw:
            report_raw = report_raw[0]
            failed_polish_search = set()
            af_filtered_search = set()
            failed_repeat_search = set()
            with open(report_raw) as fr:
                for line in fr:
                    if 'Failed\tpolish by wbc bg\t' in line:
                        chr_, pos = line.strip().split('\t', 3)[0:2]
                        failed_polish_search.add(chr_+":"+pos)
                    elif 'Failed\trepeat filter\t' in line:
                        chr_, pos = line.strip().split('\t', 3)[0:2]
                        failed_repeat_search.add(chr_+":"+pos)
                    elif 'Failed\tAF filtered\t' in line:
                        chr_, pos = line.strip().split('\t', 3)[0:2]
                        af_filtered_search.add(chr_+":"+pos)
            for each_fail in final_failed:
                if each_fail.split(',')[0] in failed_polish_search:
                    failed_polish.add(each_fail)
                elif each_fail.split(',')[0] in af_filtered_search:
                    af_filtered.add(each_fail)
                elif each_fail.split(',')[0] in failed_repeat_search:
                    failed_repeat.add(each_fail)
        filtered_rescued_detail = failed_repeat|af_filtered|failed_polish
        filtered_rescued_num = len(filtered_rescued_detail)
        if filtered_rescued_num >= 1:
            print("共挽救回 {} 个由于polish|repeat|af被过滤的突变: {}".format(filtered_rescued_num, sample))

        # save result
        TP = len(final_success) + filtered_rescued_num
        FN = len(final_failed) - filtered_rescued_num
        FP = len(final_more)
        TN = panel_size - introduced_total_len - (len(final_failed) - filtered_rescued_num)
        total = TP + TN + FN + FP
        final_failed_detail = final_failed - filtered_rescued_detail
        result_dict[os.path.basename(sample)] = dict(
            TP=TP,
            FP=FP,
            TN=TN,
            FN=FN,
            total=total,
            sensitivity=1.*TP/(TP+FN),
            specificity=1.*TN/(TN+FP),
            ppv=1.*TP/(TP+FP),
            npv=1.*TN/(TN+FN),
            accuracy=1.*(TP+TN)/total,
            # more detail saved
            introduced_num=len(introduced_list),
            report_num=len(report) + filtered_rescued_num,
            rescued_postion_num = len(rescued_more),
            rescued_polish_num = len(failed_polish),
            rescued_repeat_num = len(failed_repeat),
            rescued_af_num = len(af_filtered),
            panel_size = panel_size,
            negative_num = panel_size - introduced_total_len,
            failed_detail = '|'.join(x.strip() for x in final_failed_detail) if final_failed_detail else 'None',
            false_detail = "|".join(x.strip() for x in final_more) if final_more else 'None',
        )
    return result_dict


def report(match_result="*_runPipeline", allowed_distance=50, panel_size=1789799, vardict=True):
    method = 'vardict' if vardict else 'varscan'
    if not os.path.exists('{}_stat_result.xls'.format(method)):
        pass
    else:
        print('Use existed {}_stat_result.xls for plotting'.format(method))
        plot_stat_result('{}_stat_result.xls'.format(method))
        return
    dirs = glob.glob(match_result)
    print(dirs)
    result = dict()
    original_existed_mutation = get_original_mutation('pipeline_result', vardict=vardict)
    for each in dirs:
        result_dict = summary(each, original_existed_mutation, allowed_distance=allowed_distance,
                              panel_size=panel_size, vardict=vardict)
        result.update(result_dict)
    data = pd.DataFrame(result).T

    # sort index
    tmp_list = []
    for sample in data.index:
        af = sample.split('seed')[0].split('af')[1]
        if 'snv' in sample:
            mut_id = sample.split('af')[0].split('snv')[1]
        elif 'Snv' in sample:
            mut_id = sample.split('af')[0].split('Snv')[1]
        elif 'indel' in sample:
            mut_id = sample.split('af')[0].split('indel')[1]
        elif 'Indel' in sample:
            mut_id = sample.split('af')[0].split('Indel')[1]
        tmp_list.append((sample, af, mut_id))
    data['AF'] = [float(x[1].replace('0','0.',1)) for x in tmp_list]
    tmp_list = sorted(tmp_list, key=lambda x:(x[1], x[2]), reverse=True)
    order_sample = [x[0] for x in tmp_list]
    # print(tmp_list)

    # sort column
    col_order = [
        "AF", "TP", "FP", "TN", "FN", "total", "sensitivity", "specificity", "ppv", "npv", "accuracy",
        # more detail saved
        "introduced_num", "report_num", "rescued_postion_num", "rescued_polish_num", "rescued_repeat_num",
        "rescued_af_num", "panel_size", "negative_num", "failed_detail", "false_detail",
    ]
    data = data.loc[order_sample, col_order]
    data.index.name = 'sample'
    data.to_csv('{}_stat_result.xls'.format(method), sep='\t', header=True, index=True)
    plot_stat_result('{}_stat_result.xls'.format(method))

    sum_data = data.loc[:, ["AF", "TP", "FP", "TN", "FN", "total"]].groupby('AF').sum()
    sum_data['sensitivity'] = 1.*sum_data['TP']/(sum_data['TP']+sum_data['FN'])
    sum_data['specificity'] = 1.*sum_data['TN']/(sum_data['TN']+sum_data['FP'])
    sum_data['ppv'] = 1.*sum_data['TP']/(sum_data['TP']+sum_data['FP'])
    sum_data['npv'] = 1.*sum_data['TN']/(sum_data['TN']+sum_data['FN'])
    sum_data['accuracy'] = 1.*(sum_data['TP']+sum_data['TN'])/sum_data['total']
    sum_data.to_csv('{}_summary.xls'.format(method), sep='\t', header=True, index=True)
    # plot_stat_result('{}_summary.xls'.format(method))


def plot_stat_result(result_table):
    data = pd.read_table(result_table, header=0, )
    curve_data = data.loc[:, ['sensitivity', 'specificity', 'ppv', 'npv', 'accuracy']]
    curve_data.plot(rot=90, figsize=(12, 9))
    af_list = list()
    ind_list = [0,]
    colors = [
        'lightgray',
        'lightpink',
        'lightsalmon',
        'lightseagreen',
        'lightslategray',
        'lightsteelblue',
        'lightcyan',
        'lightskyblue',
        'lightgoldenrodyellow',
        'lightgreen',
    ]
    for ind, af in enumerate(data['AF']):
        if ind > 0:
            if af not in af_list:
                plt.axvspan(ind_list[-1], ind, facecolor=colors.pop(), alpha=0.5)
                plt.text((ind_list[-1]+ind)/2.0, 0.8, 'AF='+str(af_list[-1]))
                ind_list.append(ind)
        af_list.append(af)
    else:
        plt.axvspan(ind_list[-1], ind, facecolor=colors.pop(), alpha=0.5)
        plt.text((ind_list[-1] + ind) / 2.0, 0.8, 'AF=' + str(af_list[-1]))
    plt.savefig(result_table[:-3]+'pdf')
    plt.close()


def introduce_command(func):
    import argparse
    import inspect
    import json
    import time
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
    introduce_command(report)
