# coding=utf-8
from Bio import SeqIO
import unittest
author = 'gdq'


def split_fasta(fasta, chunk_size=100, out_prefix='my_fasta'):
    with open(fasta) as f:
        seq_num = 0
        chunk_num = 0
        out_name = out_prefix + '_' + str(chunk_num)
        file_objects = {chunk_num: open(out_name, 'w')}
        for line in f:
            if line.startswith('>'):
                seq_num += 1
                chunk_num = seq_num//(chunk_size+1)
                if chunk_num not in file_objects:
                    file_objects[chunk_num-1].close()
                    out_name = out_prefix + '_' + str(chunk_num)
                    file_objects[chunk_num] = open(out_name, 'w')
            file_objects[chunk_num].write(line)
        else:
            file_objects[chunk_num].close()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_my(self):
        split_fasta('ref_index.transcripts.fa')


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_my'))
    unittest.TextTestRunner(verbosity=2).run(suite)
