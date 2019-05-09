import os
import re


def get_all_fastq_abs_path(path_lst: tuple, exp: str = '.*-(.*?)_combined_R[12].fastq.gz'):
    # ./180824_13_180905/Sample_R18054231-180824R-Pool-02-T180701R1L2/R18054231-180824R-Pool-02-T180701R1L2_combined_R2.fastq.gz
    result_dict = dict()
    for path in path_lst:
        for root, dirs, files in os.walk(path):
            for each in files:
                match = re.fullmatch(exp, each)
                if match:
                    sample = match.groups()[0]
                    result_dict.setdefault(sample, [[], []])
                    if each.endswith('R1.fastq.gz'):
                        result_dict[sample][0].append(os.path.join(root, each))
                    else:
                        result_dict[sample][1].append(os.path.join(root, each))

    with open('fastq.info', 'w') as f:
        os.mkdir('rawdata')
        os.chdir('rawdata')
        for sample, lst in result_dict.items():
            read1 = sorted(lst[0])
            read2 = sorted(lst[1])
            f.write('{}\t{}\t{}\n'.format(sample, ';'.join(read1), ';'.join(read2)))
            # make link
            os.mkdir(sample)
            for each in read1:
                os.symlink(each, os.path.join(sample, os.path.basename(each)))
            for each in read2:
                os.symlink(each, os.path.join(sample, os.path.basename(each)))


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())

