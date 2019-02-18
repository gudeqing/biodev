import os
from glob import glob
import json
import plotly.graph_objs as go
from plotly.offline import plot as plt
import pandas as pd
import subprocess


def parse_star_final_out_log(logfile, outdir=None):
    if not outdir:
        outdir = os.path.dirname(logfile)
    sample = os.path.basename(logfile).split('.', 1)[0]
    out_table = os.path.join(outdir, '{}.alignment_summary.json'.format(sample))
    result = dict()
    with open(logfile) as f:
        for line in f:
            if line.strip().endswith(":") or not line.strip():
                continue
            desc, data = line.split('|')
            if "Number of input reads" in desc:
                result['total'] = int(data.strip())
            elif 'Average input read length' in desc:
                result['average_read_length'] = float(data.strip())
            elif 'Uniquely mapped reads number' in desc:
                result['unique_mapped'] = int(data.strip())
            elif 'Uniquely mapped reads %' in desc:
                result['unique_mapped_ratio'] = float(data.strip().strip('%'))
            elif 'Average mapped length' in desc:
                result['average_mapped_read_length'] = float(data.strip())
            elif 'Number of splices' in desc:
                result['spliced'] = int(data.strip())
            elif 'Number of splices: Annotated' in desc:
                result['annotated_spliced'] = int(data.strip())
            elif 'Number of reads mapped to multiple loci' in desc:
                result['multiple_mapped'] = int(data.strip())
            elif '% of reads mapped to multiple loci' in desc:
                result['multiple_mapped_ratio'] = float(data.strip().strip('%'))
            elif 'Number of reads mapped to too many loci' in desc:
                result['many_mapped(>10)'] = int(data.strip())
            elif '% of reads mapped to too many loci' in desc:
                result['many_mapped_ratio'] = float(data.strip().strip('%'))
            elif ' Number of chimeric reads' in desc:
                result['chimeric'] = int(data.strip())
            else:
                pass
    result['mapped_ratio'] = result['unique_mapped_ratio'] + result['multiple_mapped_ratio'] + result['many_mapped_ratio']
    result['mapped'] = result['unique_mapped'] + result['multiple_mapped'] + result['many_mapped(>10)']
    with open(out_table, 'w') as f:
        json.dump(result, f, indent=2)
    return result


def parse_samtools_flagstat_result(stat_file):
    result_list = list()
    with open(stat_file) as f:
        for line in f:
            result_list.append(line.split(' ', 3))
    return result_list


def parse_samtools_depth_result(stat_file, cov_limit=2000, step=10, outdir=None):
    depth_dict = dict()
    with open(stat_file) as f:
        for line in f:
            chr_, pos, depth = line.strip().split()
            depth_dict.setdefault(int(depth), 0)
            depth_dict[int(depth)] += 1

    merge_dict = {'0': depth_dict[0]}
    # print(sum(depth_dict.values()))

    max_depth = max(depth_dict.keys())
    if cov_limit >= max_depth:
        cov_limit = max_depth

    start = 0
    while start < cov_limit:
        start += step
        key = '{}-{}'.format(start-step+1, start)
        merge_dict.setdefault(key, 0)
        for ind in range(start-step+1, start+1):
            if ind in depth_dict:
                merge_dict[key] += depth_dict.pop(ind)

    key = '{}-{}'.format(start+1, max_depth)
    merge_dict.setdefault(key, 0)
    for ind in range(start+1, max_depth+1):
        if ind in depth_dict:
            merge_dict[key] += depth_dict.pop(ind)

    if not outdir:
        outdir = os.path.dirname(stat_file)
    sample = os.path.basename(stat_file).split('.', 1)[0]
    out_name = os.path.join(outdir, '{}.pos.depth.distribution.json'.format(sample))
    with open(out_name, 'w') as f:
        json.dump(merge_dict, f, indent=2)

    return merge_dict


def stat_target_bam(bed, bam, rRNA_bed=None, overlap=0.05,  rRNA_overlap=0.6, bedtools='bedtools', samtools='samtools',
                    outdir=None, threads=6, cov_limit=5000, step=10):
    if not outdir:
        outdir = os.path.dirname(bam)
        if not outdir:
            outdir = '.'
    sample = os.path.basename(bam).split('.', 1)[0]
    out = '{}/{}.target_region_bam.stat.txt'.format(outdir, sample)
    cmd = '{} intersect '.format(bedtools)
    cmd += '-a {} '.format(bam)
    cmd += '-b {} '.format(bed)
    cmd += '-wa '
    cmd += '-F {} '.format(overlap)
    cmd += '| '
    cmd += '{} flagstat '.format(samtools)
    cmd += '--threads {} '.format(threads/2)
    cmd += '- '
    cmd += '> {} '.format(out)
    print(cmd)
    p1 = subprocess.Popen(cmd, shell=True)

    # get rRNA bam and stat
    if rRNA_bed:
        out_rRNA = '{}/{}.rRNA_bam.stat.txt'.format(outdir, sample)
        cmd = '{} intersect '.format(bedtools)
        cmd += '-a {} '.format(bam)
        cmd += '-b {} '.format(rRNA_bed)
        cmd += '-wa -s '
        cmd += '-F {} '.format(rRNA_overlap)
        cmd += '| '
        cmd += '{} flagstat '.format(samtools)
        cmd += '--threads {} '.format(threads/2)
        cmd += '- '
        cmd += '> {} '.format(out_rRNA)
        print(cmd)
        p2 = subprocess.Popen(cmd, shell=True)

    # run samtools depth
    depth_stat = os.path.join(outdir, '{}.pos.depth'.format(sample))
    cmd = '{} depth '.format(samtools)
    cmd += '-a '
    cmd += '-d {} '.format(cov_limit)
    cmd += '-b {} '.format(bed)
    cmd += '{} '.format(bam)
    cmd += '> {} '.format(depth_stat)
    print(cmd)
    subprocess.check_call(cmd, shell=True)
    if p1.wait() != 0:
        raise Exception("Failed to split out target bam based on bed!")
    if rRNA_bed:
        if p2.wait() != 0:
            raise Exception("Failed to split out rRNA bam based on bed!")
    chr_stat = os.path.join(outdir, '{}.chromosome.alignment.stat.txt'.format(sample))
    cmd = '{} idxstats {} > {}'.format(samtools, bam, chr_stat)
    subprocess.check_call(cmd, shell=True)

    final_log = glob(os.path.join(os.path.dirname(bam), '*.Log.final.out'))[0]
    original_bam_summary = parse_star_final_out_log(final_log, outdir=outdir)
    target_bam_summary = parse_samtools_flagstat_result(out)
    target_ratio = int(target_bam_summary[4][0]) / original_bam_summary['mapped'] / 2
    depth_dict = parse_samtools_depth_result(depth_stat, cov_limit=cov_limit, step=step)
    if rRNA_bed:
        rRNA_bam_summary = parse_samtools_flagstat_result(out_rRNA)
        rRNA_ratio = int(rRNA_bam_summary[4][0]) / original_bam_summary['mapped'] / 2
        original_bam_summary['rRNA_ratio'] = round(rRNA_ratio * 100, 2)

    summary = os.path.join(outdir, '{}.target_region_alignment.summary.txt'.format(sample))
    with open(summary, 'w') as f:
        f.write('reads mapped to target region : {}/{}={:.2%}\n'.format(
            target_bam_summary[4][0],
            original_bam_summary['mapped']*2,
            target_ratio)
        )
        target_region_length = sum(depth_dict.values())
        none_zero_depth_region_length = target_region_length - depth_dict['0']
        cover_ratio = none_zero_depth_region_length/target_region_length
        f.write('non-zero depth target region: {}/{}={:.2%}\n'.format(
            none_zero_depth_region_length,
            target_region_length,
            cover_ratio
        ))
    original_bam_summary['target_capture_ratio'] = round(target_ratio*100, 2)
    original_bam_summary['target_cover_ratio'] = round(cover_ratio*100, 2)
    out_table = os.path.join(outdir, '{}.alignment_summary.json'.format(sample))
    with open(out_table, 'w') as f:
        json.dump(original_bam_summary, f, indent=2)

    return target_ratio, cover_ratio


if __name__ == '__main__':
    def introduce_command(func):
        import argparse
        import inspect
        import json
        import time
        if isinstance(func, type):
            description = func.__init__.__doc__
        else:
            description = func.__doc__
        parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
        func_args = inspect.getfullargspec(func)
        arg_names = func_args.args
        arg_defaults = func_args.defaults
        arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
        for arg, value in zip(arg_names, arg_defaults):
            if arg == 'self':
                continue
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
        if func_args.varkw is not None:
            print("warning: **keywords args is not supported, and will be neglected! ")
        args = parser.parse_args().__dict__
        with open("Argument_detail.json", 'w') as f:
            json.dump(args, f, indent=2, sort_keys=True)
        start = time.time()
        func(**args)
        print("total time: {}s".format(time.time() - start))


    introduce_command(stat_target_bam)
