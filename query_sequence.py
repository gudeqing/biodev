# -*- coding: utf-8 -*-
import sqlite3
from biocluster.workflow import Workflow
__author__ = 'gdq'


class QuerySeqWorkflow(Workflow):
    """
    query sequence of transcript or gene or cds or pep
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(QuerySeqWorkflow, self).__init__(wsheet_object)
        options = [dict(name="seq_id", type="string"),
                   dict(name="seq_type", type="string"),
                   dict(name="seq_db", type="string"), ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def query_seq(self):
        conn = sqlite3.connect(self.option('seq_db'))
        cursor = conn.cursor()
        table_name = self.option('seq_type')
        seq_id = self.option('seq_id')
        cursor.execute("SELECT * FROM {} WHERE seq_id='{}'".format(table_name, seq_id))
        seq_id, sequence = cursor.fetchone()
        with open(self.output_dir + '/tmp_seq.fa') as f:
            f.write('>'+seq_id+'\n')
            f.write(sequence+'\n')

    def run(self):
        self.query_seq()
        super(QuerySeqWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([[".", "", "临时的查询序列结果文件"], ])
        super(QuerySeqWorkflow, self).end()
