import os
import glob
import re


def get_introduced_mutation(out_dir):
    changed_reads = glob.glob('{}/*reads_to_replace.txt'.format(out_dir))
    introduced = set()
    for name in changed_reads:
        chr_, pos, ref, alt = name.split('__')[1:5]
        introduced.add('{}:{},{}>{}\n'.format(chr_, pos, ref, alt))
    return introduced


def get_introduced_mutation2(out_dir):
    vcf = glob.glob('{}/Introduced_mutation.vcf'.format(out_dir))[0]
    introduced = set()
    with open(vcf) as fr:
        for line in fr:
            chr_, pos, _, ref, alt = line.split('\t')[0:5]
            introduced.add('{}:{},{}>{}\n'.format(chr_, pos, ref, alt))
    return introduced


def generate_vcf_list(out_dir="filtered_mutation"):
    vcfs = glob.glob('{}/vcf_*'.format(out_dir))
    out_dict = dict()
    for vcf in vcfs:
        vcf_id = os.path.basename(vcf).split('_')[1]
        out_name = os.path.join(os.path.dirname(vcf), vcf_id+'.list')
        with open(vcf) as fr, open(out_name, 'w') as fw:
            lines = (x.strip().split() for x in fr if x.strip())
            for line in lines:
                chr_, pos, _, ref, alt = line
                fw.write('{}:{},{}>{}\n'.format(chr_, pos, ref, alt))
        out_dict[vcf_id] = out_name
    return out_dict


def get_original_mutation(pipeline_out):

    report_files = glob.glob('{}/*/report_vardict/*_report.xls'.format(pipeline_out))
    for each in report_files:
        cmd = "less {} | grep -v '#' | cut -f8,9,13,14 | sed 1d | sed 's/\t/:/;s/\t/,/;s/\t/>/' > {}/report_detect.list"
        cmd = cmd.format(each, os.path.dirname(each))
        os.system(cmd)
    cmd = 'cat {}/*/report_vardict/report_detect.list | sort | uniq > {}/all_original_mutation.list'.format(pipeline_out, pipeline_out)
    print(cmd)
    os.system(cmd)
    return "{}/all_original_mutation.list".format(pipeline_out)


def generate_report_list(out_dir, original_existed_mutation, allowed_distance=100):
    sample_result_dirs = glob.glob("{}/*af*seed*".format(out_dir))
    # print(sample_result_dirs)
    result_dict = dict()
    for sample in sample_result_dirs:
        report_file = glob.glob('{}/report_vardict/*_report.xls'.format(sample))
        if not report_file:
            print("Find no *_report.xls in {}".format(sample))
            continue
        else:
            report_file = report_file[0]
        cmd = "less {} | grep -v '#' | cut -f8,9,13,14 | sed 1d | sed 's/\t/:/;s/\t/,/;s/\t/>/' > {}/report_detect.list"
        cmd = cmd.format(report_file, os.path.dirname(report_file))
        os.system(cmd)
        # vcf_id = re.match('.*(\d+)_af.*seed.*',  os.path.basename(sample)).groups()[0]
        # introduced_list = open(vcf_dir_dict[vcf_id]).readlines()
        report_list = open('{}/report_detect.list'.format(os.path.dirname(report_file))).readlines()
        origin_list = open(original_existed_mutation).readlines()
        report = set(report_list) - set(origin_list)
        simulate_out_dir = out_dir.split('_runPipeline')[0]
        # introduced_list = get_introduced_mutation(simulate_out_dir + '/' + os.path.basename(sample) + '_reads')
        introduced_list = get_introduced_mutation2(simulate_out_dir + '/' + os.path.basename(sample))
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
        rescued_fail = set()
        rescued_more = set()
        for each_fail in failed:
            chr_, pos = each_fail.split(',')[0].split(":")
            for every_fail in more:
                chr_2, pos2 = every_fail.split(',')[0].split(":")
                if chr_ == chr_2 and abs(int(pos) - int(pos2)) < allowed_distance:
                    rescued_fail.add(each_fail)
                    rescued_more.add(every_fail)
        final_more = more - rescued_more
        final_failed = failed - rescued_fail
        final_success = success | rescued_fail
        # ----
        result_dict[os.path.basename(sample)] = [
            len(introduced_list),
            len(report),
            len(final_success),
            len(final_failed),
            len(final_more),
            # len(failed),
            # len(more),
        ]
        # find those filtered by polish
        report_raw = glob.glob('{}/report_vardict/*.variations.RAW.xls'.format(sample))
        failed_polish = set()
        af_filtered = set()
        if report_raw:
            report_raw = report_raw[0]
            tmp_search = set()
            af_filtered_search = set()
            with open(report_raw) as fr:
                for line in fr:
                    if 'Failed\tpolish by wbc bg\t' in line or 'Failed\trepeat filter\t' in line:
                        chr_, pos = line.strip().split('\t', 3)[0:2]
                        tmp_search.add(chr_+":"+pos)
                    elif 'Failed\tAF filtered\t' in line:
                        chr_, pos = line.strip().split('\t', 3)[0:2]
                        af_filtered_search.add(chr_+":"+pos)
            for each_fail in final_failed:
                if each_fail.split(',')[0] in tmp_search:
                    failed_polish.add(each_fail)
                elif each_fail.split(',')[0] in af_filtered_search:
                    af_filtered.add(each_fail)
            if af_filtered:
                print("af filtered in sample {}:".format(os.path.basename(sample)), len(af_filtered))
        result_dict[os.path.basename(sample)].append(len(failed_polish))
        result_dict[os.path.basename(sample)].append(len(af_filtered))

        if report:
            result_dict[os.path.basename(sample)].append(len(final_more)/float(len(report)))
            result_dict[os.path.basename(sample)].append(len(final_success)/float(len(report)))
        else:
            result_dict[os.path.basename(sample)].append(0)
            result_dict[os.path.basename(sample)].append(0)
        if introduced_list:
            result_dict[os.path.basename(sample)].append(len(final_success)/float(len(introduced_list)))
        else:
            result_dict[os.path.basename(sample)].append(0)
        result_dict[os.path.basename(sample)].append(final_failed)
        result_dict[os.path.basename(sample)].append(final_more)
        result_dict[os.path.basename(sample)].append(failed_polish)
    return result_dict


def main():
    dirs = glob.glob('*_runPipeline')
    print(dirs)
    result = dict()
    original_existed_mutation = get_original_mutation('pipeline_result')
    for each in dirs:
        result_dict = generate_report_list(each, original_existed_mutation, allowed_distance=50)
        result.update(result_dict)
    fw = open('stat_result.xls', 'w')
    header = ['sample', 'introduced_num', 'report_num', 'success_num', 'failed_num', 'false_num',
              'failed_polish|repeat_num', 'af_filtered',
              'false_num/report_num', 'success_num/report_num', 'success_num/introduced_num',
              'failed_detail(failed_polish excluded)',
              'false_detail',
              # 'failed_polish_detail'
              ]
    fw.write('\t'.join(header) + '\n')
    tmp_list = []
    for sample in result.keys():
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
    # print(tmp_list)
    tmp_list = sorted(tmp_list, key=lambda x:(x[1], x[2]), reverse=True)
    order_sample = [x[0] for x in tmp_list]
    for sample in order_sample:
        tmp_list = result[sample][:-3]
        final_failed = result[sample][-3] if result[sample][-3] else {'None'}
        final_more = result[sample][-2] if result[sample][-2] else {'None'}
        failed_polish = result[sample][-1] if result[sample][-1] else {'None'}
        tmp_set = final_failed-failed_polish if (final_failed-failed_polish) else {'None'}
        tmp_list.append('; '.join(x.strip() for x in tmp_set))
        tmp_list.append('; '.join(x.strip() for x in final_more))
        # tmp_list.append('; '.join(x.strip() for x in failed_polish))
        fw.write('{}\t{}\n'.format(sample, '\t'.join(str(x) for x in tmp_list)))
    fw.close()


if __name__ == '__main__':
    if not os.path.exists('stat_result.xls'):
        main()
    else:
        print('Use existed stat_result.xls for plotting')
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import pandas as pd
    data = pd.read_table('stat_result.xls', header=0)
    data.loc[data['report_num']<30, 'report_num'] = 0  # small report num will greatly affect plotting
    curve_data = data.loc[:, ['false_num/report_num', 'success_num/report_num', 'success_num/introduced_num']]
    data['success_num + failed_polish|repeat_num + af_filtered'] = data['success_num']+data['failed_polish|repeat_num']+data['af_filtered']
    # data['failed_num - failed_polish|repeat_num - af_filtered'] = data['failed_num'] - data['failed_polish|repeat_num'] - data['af_filtered']
    curve_data['(success_num + failed_polish|repeat_num + af_filtered/report_num)'] = data['success_num + failed_polish|repeat_num + af_filtered']/data['report_num']
    curve_data['(success_num + failed_polish|repeat_num + af_filtered/introduced_num)'] = data['success_num + failed_polish|repeat_num + af_filtered']/data['introduced_num']
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
    for ind, sample in enumerate(data.iloc[:, 0]):
        af = sample.split('af')[1].split('seed')[0]
        if ind > 0:
            if af not in af_list:
                plt.axvspan(ind_list[-1], ind, facecolor=colors.pop(), alpha=0.5)
                plt.text((ind_list[-1]+ind)/2.0, 0.8, 'AF='+af_list[-1])
                ind_list.append(ind)
        af_list.append(af)
    else:
        plt.axvspan(ind_list[-1], ind, facecolor=colors.pop(), alpha=0.5)
        plt.text((ind_list[-1] + ind) / 2.0, 0.8, 'AF=' + af_list[-1])
    plt.savefig('plot_stat.pdf')
