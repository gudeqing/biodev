# coding=utf-8
import argparse
import subprocess
import shlex
import pandas as pd
import re



def run_blast(target, query, top=1, out="diamond.out.txt", format=6,
              sensitive="sensitive", p=12, evalue=1e-3, blast_type='blastp',
              diamond="/mnt/ilustre/users/sanger-dev/app/bioinfo/align/diamond-0.8.35/diamond"):
    """
    比对的目的是为了能够分配一个已知的TF_id，后续靶基因预测可以根据这个已知的TF_id获得转录因子靶向Motif.
    :param target: file for build database file
    :param query: input query file
    :param top: report alignments within this percentage range of top alignment score (overrides --max-target-seqs)
    :param out: output file
    :param format: output format
    :param sensitive:  sensitive, (default: fast), more-sensitive
    :param p: number of CPU threads
    :param evalue: maximum e-value to report alignments
    :param blast_type: blastp or blastx
    :param diamond: where is diamond
    :return:
    """
    # make db cmd
    build_db_cmd = "{diamond} ".format(diamond=diamond)
    build_db_cmd += "makedb --in {} ".format(target)
    build_db_cmd += "-d {} ".format('seqdb')
    # blast cmd
    blast_cmd = "{diamond} ".format(diamond=diamond)
    blast_cmd += "{blast_type} ".format(blast_type=blast_type)
    blast_cmd += "-d {db} ".format(db="seqdb")
    blast_cmd += "-q {query} ".format(query=query)
    blast_cmd += "-o {out} ".format(out=out)
    blast_cmd += "-k {hit_num} ".format(hit_num=top)
    blast_cmd += "-p {threads} ".format(threads=p)
    blast_cmd += "-f {format} ".format(format=format)
    blast_cmd += "-e {evalue} ".format(evalue=evalue)
    blast_cmd += "--{} ".format(sensitive)
    # run cmd
    subprocess.check_call(shlex.split(build_db_cmd))
    subprocess.check_call(shlex.split(blast_cmd))


def parse_args():
    parser = argparse.ArgumentParser(description="""
    Check http://github.com/bbuchfink/diamond for updates.
    """, )
    parser.add_argument('-target', help="target file for db construction")
    parser.add_argument('-query', help="query file")
    parser.add_argument('-p', default=12, type=int, help="""
    number of CPU threads
    """)
    parser.add_argument('-o', help="output file", default="diamond.out.txt")
    parser.add_argument('-k', default=1, type=int, help="""
    report alignments within this percentage range of top alignment score (overrides --max-target-seqs)
    """)
    parser.add_argument('-f', default=6, type=int, help=""" 
    output format:
    5   = BLAST XML
	6   = BLAST tabular
	100 = DIAMOND alignment archive (DAA)
	101 = SAM
    """)
    parser.add_argument('-e', default=0.001, type=float, help="""
    maximum e-value to report alignments
    """)
    parser.add_argument('-sensitive', default="sensitive", help="""
    sensitive, (default: fast), more-sensitive
    """)
    parser.add_argument('-diamond', help="where is diamond", default="~/app/program/Python/bin/diamond")
    parser.add_argument('-blast_type', default="blastp", help="blastp or blastx")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    run_blast(
        target=args.target,
        query=args.query,
        top=args.k,
        out=args.o,
        format=args.f,
        sensitive=args.sensitive,
        p=args.p,
        evalue=args.e,
        blast_type=args.blast_type,
        diamond=args.diamond,
    )
    