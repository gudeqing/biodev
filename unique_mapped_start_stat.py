import pysam
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-b', type=str, required=True, help='bam file')
parser.add_argument('-o', type=str, default=None, help="output file")
args = parser.parse_args()
bam_file = args.b
stat_out = args.o
if stat_out is None:
    stat_out = os.path.basename(bam_file) + '.unique_mapped_start.stat.txt'
pos_reads_dict = dict()
read_num = 0
with pysam.AlignmentFile(bam_file, "rb") as bam, open(stat_out, 'w') as fw:
    for read in bam.fetch():
        if read.is_proper_pair:
            read_num += 1
            key = read.reference_name + ":" +str(read.is_reverse) + ':' + str(read.pos)
            pos_reads_dict.setdefault(key, 0)
            pos_reads_dict[key] += 1
    unique_ratio = len(pos_reads_dict)/float(read_num)
    print("total unique mapped start positions: {}({:.2%})".format(
        len(pos_reads_dict), unique_ratio)
    )
    fw.write('# {}/{} = {:.2%}\n'.format(len(pos_reads_dict), read_num, unique_ratio))
    for pos, num in pos_reads_dict.items():
        fw.write("{}\t{}\n".format(pos, num))
