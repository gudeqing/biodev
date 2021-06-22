import sys
import subprocess
import argparse


def open_vcf(vcf, mode='r'):
    if vcf.endswith('.gz'):
        import gzip
        return open(gzip.open(vcf, mode=mode))
    else:
        return open(vcf, mode=mode)


# def get_sample_id(vcf):
#     with open_vcf(vcf) as f:
#         for line in f:
#             if line.startswith('#CHROM'):
#                 lst = line.strip().split('\t')
#                 return lst[-1]
#

def get_target_vcf(jc_vcf, out_vcf, sample_num=1):
    with open_vcf(jc_vcf) as fr, open(out_vcf, 'w') as fw:
        for line in fr:
            if line.startswith('##'):
                fw.write(line)
            else:
                lst = line.strip().split('\t')
                newline = '\t'.join(lst[:9+sample_num]) + '\n'
                fw.write(newline)


def joint_call(my_vcf_list, reference, out_vcf, other_vcf_list=None, threads=0,
               dbsnp=None, call_conf=30, emit_conf=30, emit_mode='variant',
               genotype_model='coalescent', max_alt_alleles=100):
    if genotype_model not in ['coalescent', 'multinomial']:
        raise Exception("genotype model only supports ['coalescent', 'multinomial'] ")
    if emit_mode not in ["variant", "confident", "all"]:
        raise Exception('emit_mode only supports ["variant", "confident", "all"] ')
    if not other_vcf_list:
        other_vcf_list = []
    # my_samples = [x for x in get_sample_id(my_vcfs)]
    cmd = 'sentieon driver '
    if threads > 0:
        cmd += '-t {threads} '.format(threads=threads)
    cmd += '-r {reference} '.format(reference=reference)
    cmd += '--algo GVCFtyper '
    if dbsnp:
        cmd += '--dbsnp {dbsnp} '.format(dbsnp=dbsnp)
    cmd += '--call_conf {call_conf} '.format(call_conf=call_conf)
    cmd += '--emit_conf {emit_conf} '.format(emit_conf=emit_conf)
    cmd += '--genotype_model {genotype_model} '.format(genotype_model=genotype_model)
    cmd += '--max_alt_alleles {max_alt_alleles} '.format(max_alt_alleles=max_alt_alleles)
    cmd += '-v ' + ' -v '.join(my_vcf_list + other_vcf_list)
    cmd += ' all.joint.vcf'
    print(cmd)
    subprocess.check_call(cmd, shell=True)
    get_target_vcf("all.joint.vcf", out_vcf, sample_num=len(my_vcf_list))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-my_vcfs', type=str, required=True, nargs='+', help='input gVcf files')
    parser.add_argument('-other_vcfs', type=str, required=False, nargs='+', help='input gVcf files')
    parser.add_argument('-r', type=str, required=True, help='reference fasta file')
    parser.add_argument('-o', type=str, required=True, help='output vcf file')
    parser.add_argument('-t', type=int, required=False, default=0, help='threads number')
    parser.add_argument('-d', type=str, required=False, help='dbSNP file')
    parser.add_argument('--call_conf', type=int, default=30, help='call confidence level(default=30)')
    parser.add_argument('--emit_conf', type=int, default=30, help='emit confidence level(default=30)')
    parser.add_argument('--emit_mode', type=str, default='variant', help='emit mode, one of [variant, confident, all]')
    parser.add_argument('--genotype_model', type=str, default='coalescent', help='genotype model, one of [coalescent, multinomial]')
    parser.add_argument('--max_alt_alleles', type=int, default=100, help='maximum number of alternate alleles(default=100)')
    args = parser.parse_args()
    try:
        joint_call(args.my_vcfs, args.r, args.o, other_vcf_list=args.other_vcfs, threads=args.t,
                   dbsnp=args.d, call_conf=args.call_conf, emit_conf=args.emit_conf,
                   genotype_model=args.genotype_model, max_alt_alleles=args.max_alt_alleles)
    except Exception as e:
        subprocess.call('rm -fr *vcf', shell=True)
        print(e)
        exit(1)
    finally:
        subprocess.call('rm all.*vcf*', shell=True)
