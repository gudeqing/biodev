import re


def fasta_stat(fasta_file, targets=None, strip_version=True):
    """
    提取gene和transcript序列信息，都会用这个函数
    :param fasta_file:
    :return: dict, ｛seq_id: sequence｝
    """
    seq = dict()
    if targets:
        target_set = set(x.strip() for x in open(targets))
    else:
        target_set = set()
    match_name = re.compile(r'>([^\s]+)').match
    with open(fasta_file, 'r+') as fasta:
        j = 0
        seq_id, sequence = '', ''
        for line in fasta:
            if line.startswith('#') or not line.strip():
                continue
            if line.startswith('>'):
                j += 1
                if j > 1:
                    seq[seq_id] = sequence
                    sequence = ''
                # get seq name
                seq_id = match_name(line).group(1)
                if strip_version:
                    seq_id = seq_id.rsplit('.', 1)[0]
            else:
                sequence += line.strip()
        else:
            # save the last sequence
            seq[seq_id] = sequence
    if not seq:
        print('提取序列信息为空')
    print("从{}共统计出{}条序列信息".format(fasta_file, len(seq)))
    with open('seq_stat.txt', 'w') as f:
        f.write('name\tlength\n')
        for seq_id, sequence in seq.items():
            if targets is None:
                f.write(seq_id + '\t' + str(len(sequence)) + '\n')
            else:
                if seq_id in target_set:
                    f.write(seq_id + '\t' + str(len(sequence)) + '\n')
    return seq


def extract_chromosome_fasta(genome_fasta):
    with open(genome_fasta) as fr:
        fw = None
        for line in fr:
            if line.startswith('>'):
                if fw:
                    fw.close()
                chr_name = line[1:].split()[0]
                fw = open(chr_name+'.fa', 'w')
                fw.write(f'>{chr_name}\n')
            else:
                fw.write(line)
        fw.close()


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['fasta_stat', 'extract_chromosome_fasta'])
