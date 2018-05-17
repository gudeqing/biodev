import re
import argparse
from collections import defaultdict
__author__ = 'gdq'


class ParseNCBITaxonomy(object):
    """
    parse taxonomy file from ncbi
    input files are from ftp://ftp.ncbi.nih.gov/pub/taxonomy/
    """
    def __init__(self, node_file, name_file, merged_file):
        self.version = 0.1
        self.node_file = node_file
        self.name_file = name_file
        self.merged_file = merged_file

    def parse_nodes(self, node_file=None):
        """
        < nodes.dmp > This file represents taxonomy nodes.
        The description for each node includes the following fields:
        tax_id                                  -- node id in GenBank taxonomy database
        parent tax_id                           -- parent node id in GenBank taxonomy database
        rank                                    -- rank of this node (superkingdom, kingdom, ...)
        embl code                               -- locus-name prefix; not unique
        division id                             -- see division.dmp file
        inherited div flag  (1 or 0)            -- 1 if node inherits division from parent
        genetic code id                         -- see gencode.dmp file
        inherited GC  flag  (1 or 0)            -- 1 if node inherits genetic code from parent
        mitochondrial genetic code id           -- see gencode.dmp file
        inherited MGC flag  (1 or 0)            -- 1 if node inherits mitochondrial gencode from parent
        GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
        hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
        comments                                -- free-text comments and citations
        :rtype:dict
        :param node_file: nodes.dmp
        :return: dict, tax_id: {'parent_tax_id':123, 'rank':123,...}
        """
        fileds = ['tax_id', 'parent_tax_id', 'rank', 'embl_code', 'division_id',
                  'inherited_div_flag', 'genetic_code_id', 'inherited_GC_flag',
                  'mitochondrial_genetic_code_id', 'inherited_MGC_flag',
                  'GenBank_hidden_flag', 'hidden_subtree_root', 'comments']
        node_dict = dict()
        if not node_file:
            node_file = self.node_file
        with open(node_file) as f:
            for line in f:
                if not line.strip():
                    continue
                tmp_list = line.strip().replace('\t', '').split('|')[:-1]
                # tmp_dict = dict(zip(fileds, tmp_list))
                tmp_dict = dict(zip(fileds[0:3], tmp_list))
                tax_id = tmp_dict.pop('tax_id')
                node_dict[tax_id] = tmp_dict
        return node_dict

    def parse_names(self, name_file=None, target='scientific name'):
        """
        < names.dmp > Taxonomy names file has these fields:
        tax_id					-- the id of node associated with this name
        name_txt				-- name itself
        unique name				-- the unique variant of this name if name not unique
        name class				-- (synonym, common name, ...)
        :param name_file: names.dmp
        :param target: str, target name to be extracted
        :return: dict, tax_id: {'name_txt': xxx, 'name_class': xxx}
        """
        fileds = ['tax_id', 'name_txt', 'unique_name', 'name_class']
        name_dict = dict()
        if not name_file:
            name_file = self.name_file
        with open(name_file) as f:
            for line in f:
                if not line.strip():
                    continue
                tmp_list = line.strip().replace('\t', '').split('|')[:-1]
                if tmp_list[-1] == target:
                    tmp_dict = dict(zip(fileds, tmp_list))
                    tax_id = tmp_dict.pop('tax_id')
                    name_dict[tax_id] = tmp_dict
        return name_dict

    def parse_merged(self, merged_file=None):
        """
        parse merged.dmp file into dict
        :param merged_file:
        Merged nodes file (merged.dmp):
        old_tax_id         -- id of nodes which has been merged
        new_tax_id         -- id of nodes which is result of merging
        ----------------------------------------------
        12      |       74109   |
        30      |       29      |
        36      |       184914  |
        37      |       42      |
        ----------------------------------------------
        :return: {old_tax_id: new_tax_id}
        """
        if not merged_file:
            merged_file = self.merged_file
        old2new = dict()
        with open(merged_file) as f:
            for line in f:
                if not line.strip():
                    continue
                tmp_list = line.strip('\n').replace('\t', '').split('|')[:-1]
                old2new[tmp_list[0]] = tmp_list[1]
        return old2new

    def get_children_nodes(self, parent_nodes, node_file=None, node_dict=None):
        """
        Obtain all child node for a specific parent node
        :param node_file: nodes.dmp
        :param parent_nodes: str or list, example ['32199',]
        :param node_dict: result from self.parse_nodes.
        :return: set, {child_tax_id1, child_tax_id2,...}
        """
        if not node_dict:
            if not node_file:
                node_file = self.node_file
            node_dict = self.parse_nodes(node_file)

        if isinstance(parent_nodes, str):
            parent_nodes = [parent_nodes]
        children_set = set(parent_nodes)
        children_set_update = children_set.update
        while True:
            children = set()
            for n in node_dict:
                if node_dict[n]['parent_tax_id'] in parent_nodes:
                    children.add(n)
            if children:
                children_set_update(children)
                parent_nodes = set(children)
            else:
                break
        return children_set

    def children_set_calculator(self, reg, node_dict=None):
        """
        Only + and - are allowed for set calculation.
        :param reg: such as '2759-4751-33090' or simply '2759'
        :param node_dict: node dict
        :return: set
        """
        operations = re.findall(r'\d+([-+])', reg)
        node_list = re.split('[-+]', reg)
        if not operations:
            return self.get_children_nodes(node_list, node_dict=node_dict)
        else:
            children_set = self.get_children_nodes(node_list[0], node_dict=node_dict)
            for oper, node in zip(operations, node_list[1:]):
                if oper == '+' or oper == "&":
                    children_set = children_set & self.get_children_nodes(node, node_dict=node_dict)
                elif oper == '-':
                    children_set = children_set - self.get_children_nodes(node, node_dict=node_dict)
            return children_set

    def get_tax_id_of_taxon(self, node_file=None, taxon='all'):
        """
        query tax_id of a taxon
        :param taxon: could be anyone of:  class cohort family forma genus infraclass infraorder kingdom order
         parvorder phylum species 'species group' 'species subgroup' subclass subfamily subgenus subkingdom
         suborder subphylum subspecies subtribe superclass superfamily superkingdom superorder
         superphylum tribe varietas
        :param node_file: nodes.dmp
        :return: dict
        """
        if not node_file:
            node_file = self.node_file
        node_dict = self.parse_nodes(node_file)
        taxon_tax_id_dict = defaultdict(set)
        for n in node_dict:
            one_taxon = node_dict[n]['rank']
            taxon_tax_id_dict[one_taxon].add(n)
        if taxon == 'all':
            return taxon_tax_id_dict
        else:
            if taxon in taxon_tax_id_dict:
                return taxon_tax_id_dict[taxon]
            else:
                print(taxon + ' not found')

    @staticmethod
    def prot_accession2tax_id(prot2tax_file):
        """
        parse prot.accession2taxid.gz into dict.
        ---------------------------------------------------
        prot2tax_file:
        accession       accession.version       taxid   gi
        APZ74649        APZ74649.1      36984   1137646701
        AQT41667        AQT41667.1      1686310 1150388099
        WP_080502060    WP_080502060.1  95486   1169627919
        ......
        ---------------------------------------------------
        :param prot2tax_file:
        :return: dict, {prot_id: tax_id, ...}, Protein with tax_id=0 will be filtered.
        """
        with open(prot2tax_file) as f:
            _ = f.readline()
            prot2tax_dict = dict()
            for line in f:
                if not line.strip():
                    continue
                tmp_list = line.strip('\n').split('\t')
                if not tmp_list[2] == '0':
                    # protein are not classified.
                    prot2tax_dict[tmp_list[1]] = tmp_list[2]
        return prot2tax_dict

    @staticmethod
    def fasta_parser(fasta):
        """
        yield seq_name, sequence
        :param fasta:
        :return: seqname, sequence('\n' were kept)
        """
        with open(fasta) as f:
            match_pattern = re.compile(r'>([^\s]+)').match
            seq_name, sequence, seq_num, desc = '', '',  0, ''
            for line in f:
                if (not line.strip()) or line.startswith('#'):
                    continue
                if line.startswith('>'):
                    seq_num += 1
                    if seq_num >= 2:
                        yield seq_name, desc, sequence
                    sequence = ''
                    seq_name = match_pattern(line).group(1)
                    desc = line.lstrip('>'+seq_name).strip()
                else:
                    sequence += line  # "\n" still contained
            else:
                yield seq_name, desc, sequence

    @staticmethod
    def get_seq_names_from_fasta(fasta):
        match_pattern = re.compile(r'>([^\s]+)').match
        with open(fasta) as f:
            seq_names = set()
            for line in f:
                if line.startswith('>'):
                    m = match_pattern(line)
                    if m:
                        name = m.group(1)
                        seq_names.add(name)
        return seq_names

    @staticmethod
    def get_kingdom_taxid():
        """
        List tax_id of kingdom and super-kingdom
        :return: dict
        """
        # three kingdom
        kingdom_dict = dict(
            fungi_tax_id='4751',
            viridiplantae_tax_id='33090',
            metazoa_tax_id='33208',
            # five super-kingdom
            bacteria_tax_id='2',
            archaea_tax_id='2157',
            eukaryota_tax_id='2759',
            viruses_tax_id='10239',
            viroids_tax_id='12884',
            protist='2759-4751-33090-33208')
        return kingdom_dict

    def nr_split(self, name2taxid_dict, prot2tax_file, fasta_file, name_file=None,
                 node_file=None, merged_file=None):
        """
        Split fasta file of NR.fa according to input name2taxid_dict. For each tax_id in the dict,
        its children tax nodes will be find first, and proteins associated with these children tax
        nodes will be found subsequently. finally, fasta for these proteins will be extracted.
        :param name2taxid_dict: dict, such as {'scientific_name':tax_id, ...}
        :param prot2tax_file: unzipped file of prot.accession2taxid.gz
        :param fasta_file: NR fasta file from ncbi
        :param name_file: names.dmp file from ncbi
        :param node_file: nodes.dmp file from ncbi
        :param merged_file: merged.dmp file from ncbi
        :return: None, but fasta files will be generated.
        """
        f = open('not_classified_protein.list', 'w')
        f.write('not_classified_id\twhy_not_classified\n')

        prot2tax_dict = self.prot_accession2tax_id(prot2tax_file)
        name_dict = self.parse_names(name_file=name_file, target='scientific name')
        old2new_dict = self.parse_merged(merged_file=merged_file)

        # to reduce the size of dict, and replace the old id with new id
        seq2tax = dict()
        seq_set = self.get_seq_names_from_fasta(fasta_file)
        for key in seq_set:
            tmp_tax_id = prot2tax_dict.get(key)
            if not tmp_tax_id:
                f.write(key + '\tIt was not found in prot.accession2taxid\n')
                continue
            else:
                # if tmp_tax_id is old, then new tax id will replace it
                new_tax_id = old2new_dict.get(tmp_tax_id)
                if new_tax_id:
                    tmp_tax_id = new_tax_id
                tmp_tax_name = name_dict.get(tmp_tax_id)
                if not tmp_tax_name:
                    sci_name = 'None'
                    print('tax_id {} of {} was not in names.dmp'.format(tmp_tax_id, key))
                else:
                    sci_name = tmp_tax_name['name_txt']
                seq2tax[key] = (tmp_tax_id, sci_name)

        # get tax ids for each target and generate file objects
        target_children_taxid_dict = dict()
        file_obj_dict = dict()
        node_dict = self.parse_nodes(node_file=node_file)
        for tax_name in name2taxid_dict:
            tmp_tax_name = tax_name.split('_tax_id')[0]
            file_obj_dict[tax_name] = [open(tmp_tax_name+'.fa', 'w'), open(tmp_tax_name+'.prot_taxId_scientificName.txt', 'w')]
            tmp_tax_id = name2taxid_dict[tax_name]
            target_children_taxid_dict[tax_name] = self.children_set_calculator(tmp_tax_id, node_dict=node_dict)

        # begin splitting
        for seq_name, desc, sequence in self.fasta_parser(fasta_file):
            if seq_name in seq2tax:
                prot_tax_id, sci_name = seq2tax.pop(seq_name)
            else:
                continue
            # write info
            classified_num = 0
            for tax_name in target_children_taxid_dict:
                if prot_tax_id in target_children_taxid_dict[tax_name]:
                    # if scientific name and tax id are wanted, using the following line
                    # file_obj_dict[tax_name].write('>'+seq_name+' ['+prot_tax_id+':'+sci_name+'] '+desc+'\n '+sequence)
                    file_obj_dict[tax_name][0].write('>'+seq_name+' '+desc+'\n'+sequence)
                    file_obj_dict[tax_name][1].write('{}\t{}\t{}\n'.format(seq_name, prot_tax_id, sci_name))
                    classified_num += 1
            if classified_num == 0:
                f.write(seq_name + '\tits tax_id {} was not classified as any of target '
                                   'taxonomies\n'.format(prot_tax_id))
        f.close()
        # close file objects
        for obj_name in file_obj_dict:
            file_obj_dict[obj_name][0].close()
            file_obj_dict[obj_name][1].close()


if __name__ == "__main__":
    import time
    start_time = time.time()
    parser = argparse.ArgumentParser(description="parse taxonomy from ftp://ftp.ncbi.nih.gov/pub/taxonomy/")
    parser.add_argument('-node', type=str, required=True, help="nodes.dmp file from ncbi")
    parser.add_argument('-name', type=str, required=True, help='names.dmp file from ncbi')
    parser.add_argument('-merged', type=str, required=True, help='merged.dmp file from ncbi')
    parser.add_argument('-pro2tax', type=str, required=True, help="unzipped file of prot.accession2taxid.gz")
    parser.add_argument('-nr', metavar="NR.fa", type=str, required=True, help="NR fasta file from ncbi")
    args = parser.parse_args()
    worker = ParseNCBITaxonomy(args.node, args.name, args.merged)
    target2taxid_dict = worker.get_kingdom_taxid()
    worker.nr_split(target2taxid_dict, args.pro2tax, args.nr, name_file=args.name,
                    node_file=args.node, merged_file=args.merged)
    end_time = time.time()
    print("It takes {} seconds to do the work".format(end_time-start_time))

