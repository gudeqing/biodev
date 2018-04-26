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
