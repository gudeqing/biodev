import csv
import glob
import sys
import os


def coverage(files, out):
    files = glob.glob(files)
    files = sorted(files)
    print(files)
    gens = [csv.reader(open(x)) for x in files]
    names = [os.path.basename(x).split('.target_bed_coverage_metrics')[0] for x in files]
    merged = [['metrics']+names]
    for lines in zip(*gens):
        lst = lines[0][2:4]
        for line in lines[1:]:
            lst += [line[3]]
        # print(lst)
        merged.append(lst)
    with open(out, 'w') as f:
        for each in merged:
            f.write('\t'.join(each)+'\n')


def mapping(files, out):
    files = glob.glob(files)
    files = sorted(files)
    print(files)
    gens = [csv.reader(open(x)) for x in files]
    names = [os.path.basename(x).split('.mapping_metrics')[0] for x in files]
    merged = [['metrics'] + names]
    for i, lines in enumerate(zip(*gens)):
        if not lines[0][0].endswith('SUMMARY'):
            break
        if i == 0:
            # remove 100
            lines = [x[:-1] for x in lines]
        lst = [lines[0][2], lines[0][-1]]
        for line in lines[1:]:
            lst += [line[-1]]
        # print(lst)
        merged.append(lst)
    with open(out, 'w') as f:
        for each in merged:
            f.write('\t'.join(each) + '\n')


if sys.argv[1] == 'cov':
    coverage(sys.argv[2], sys.argv[3])
elif sys.argv[1] == 'map':
    mapping(sys.argv[2], sys.argv[3])
else:
    print('please use "cov" or "map" for the second argument')
