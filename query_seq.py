# coding=utf-8
import re
import sqlite3
__author__ = 'gdq'


def get_cds_seq(cds_path):
    """
    从已经下载好的cds序列文件中提取转录本对应的cds序列。
    :param cds_path: cds序列文件的绝对路径
    :return: dict，转录本id：{"name": cds_id, "sequence": cds_sequence,
                            "sequence_length": len(cds_sequence)}
    """
    cds_dict = dict()
    cds_pattern_match = re.compile(r'>([^\s]+)').match
    with open(cds_path, 'r+') as f:
        j = 0
        trans_id, cds_id, cds_sequence = '', '', ''
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                j += 1
                if j > 1:
                    seq_len = len(cds_sequence)
                    cds_dict[trans_id] = dict(name=cds_id, sequence=cds_sequence,
                                              sequence_length=seq_len)
                    cds_dict[cds_id] = dict(name=cds_id, sequence=cds_sequence,
                                            sequence_length=seq_len)
                cds_id = cds_pattern_match(line).group(1)
                if '.' in cds_id:
                    trans_id = cds_id[:cds_id.rfind('.')]
                else:
                    trans_id = cds_id
                # cds_id and trans_id will be both saved as gtf may use either one of them
                cds_sequence = ''
            else:
                cds_sequence += line.strip()
        else:
            seq_len = len(cds_sequence)
            cds_dict[trans_id] = dict(name=cds_id, sequence=cds_sequence,
                                      sequence_length=seq_len)
            cds_dict[cds_id] = dict(name=cds_id, sequence=cds_sequence,
                                    sequence_length=seq_len)
    if not cds_dict:
        print('提取cds序列信息为空')
    print("共统计出{}条转录本的cds信息".format(len(cds_dict)))
    return cds_dict


def get_pep_seq(pep_path):
    """
    get transcript's pep info, including protein sequence
    :param pep_path:
    :return: dict, trans_id={"name": pep_id, "sequence": pep_sequence,
                             "sequence_length": len(pep_sequence)}
    """
    pep_dict = dict()
    j, trans_id, trans_id_else, pep_sequence, pep_id = 0, '', '', '', ''
    pep_pattern = re.compile(r'>([^\s]+)')
    trans_pattern = re.compile(r'transcript:([^\s]+)')

    with open(pep_path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                j += 1
                if j > 1:
                    seq_len = len(pep_sequence)
                    pep_dict[trans_id] = dict(name=pep_id, sequence=pep_sequence,
                                              sequence_length=seq_len)
                    pep_dict[trans_id_else] = dict(name=pep_id, sequence=pep_sequence,
                                                   sequence_length=seq_len)
                pep_id = pep_pattern.match(line).group(1)
                trans_id = trans_pattern.search(line).group(1)
                if '.' in trans_id:
                    trans_id_else = trans_id[:trans_id.rfind('.')]
                else:
                    trans_id_else = trans_id
                pep_sequence = ''
            else:
                pep_sequence += line.strip()
        else:
            seq_len = len(pep_sequence)
            pep_dict[trans_id] = dict(name=pep_id, sequence=pep_sequence,
                                      sequence_length=seq_len)
            pep_dict[trans_id_else] = dict(name=pep_id, sequence=pep_sequence,
                                           sequence_length=seq_len)
    if not pep_dict:
        print('提取蛋白序列信息为空')
    print("共统计出{}条转录本的蛋白序列信息".format(len(pep_dict)))
    return pep_dict


def fasta2dict(fasta_file):
    """
    提取gene和transcript序列信息，都会用这个函数
    :param fasta_file:
    :return: dict, ｛seq_id: sequence｝
    """
    seq = dict()
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
                if '(' in seq_id:
                    seq_id = seq_id.split('(')[0]
                elif '|' in seq_id:
                    seq_id = seq_id.split('|')[0]
                else:
                    pass
            else:
                sequence += line.strip()
        else:
            # save the last sequence
            seq[seq_id] = sequence
    if not seq:
        print('提取序列信息为空')
    print("从{}共统计出{}条序列信息".format(fasta_file, len(seq)))
    return seq


# table_names = ['GeneSeqs', 'TranscriptSeqs', 'CDSSeqs', 'PEPSeqs']
# seq_dicts = [gene_dict, transcript_dict, cds_dict, pep_dict]


def build_seq_database(seq_dicts):
    """
    :param seq_dicts: {table_name, seq_dict}
    :return:
    """
    conn = sqlite3.connect('refrna_seqs.db')
    cursor = conn.cursor()
    # Create table

    for table_name in seq_dicts:
        cursor.execute('DROP TABLE IF EXISTS {}'.format(table_name))
        cursor.execute('CREATE TABLE {} (seq_id text, sequence text)'.format(table_name))
        seq_dict = seq_dicts[table_name]
        for seq_id in seq_dict:
            seq = seq_dict[seq_id]  # seq_id is transcript_id or gene_id
            if type(seq) == dict:
                seq = seq_dict[seq_id]['sequence']
            cursor.execute("INSERT INTO {} VALUES ('{}', '{}')".format(table_name, seq_id, seq))
    conn.commit()
    conn.close()


if __name__ == "__main__":
    base_path = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/"
    pep_path = base_path + "cds/Mus_musculus.GRCm38.pep.all.fa"
    cds_path = base_path + "cds/Mus_musculus.GRCm38.cds.all.fa"
    gene_path = "/mnt/ilustre/users/sanger-dev/workspace/20170724/Single_gene_fa_2/GeneFa/output/gene.fa"
    transcript_path = "/mnt/ilustre/users/sanger-dev/workspace/20170706/Single_rsem_stringtie_mouse_total_2/Express1/TranscriptAbstract/output/exons.fa"
    pep = get_pep_seq(pep_path)
    cds = get_cds_seq(cds_path)
    gene = fasta2dict(gene_path)
    transcript = fasta2dict(transcript_path)
    seq_dicts = dict(pep=pep, cds=cds, gene=gene, transcript=transcript)
    build_seq_database(seq_dicts)

