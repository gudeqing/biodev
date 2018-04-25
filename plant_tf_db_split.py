import argparse
import sys
import json
import re

plant_db_pep_file = sys.argv[1]
our_annot_json = sys.argv[2]

# get our species info
with open(our_annot_json) as f:
    our_annot_dict = json.load(f)
    our_species = {x.lower().replace('_', ' ') for x in our_annot_dict.keys()}

# get plant db species info and split fasta
match_pattern = re.compile(r'>.*?\s(.*?)\|.*')
plant_db_species = set()
with open(plant_db_pep_file) as fr:
    first_line = fr.readline()
    species = match_pattern.match(first_line).groups()[0].lower()
    plant_db_species.add(species)
    fw = open(species.replace(" ", "_")+'.pep.fa', 'w')
    fw.write(first_line)
    seq = ''
    for line in fr:
        find = match_pattern.match(line)
        if find:
            tmp = find.groups()[0].lower()
            fw.write(seq)
            if tmp != species:
                plant_db_species.add(tmp)
                fw.close()
                fw = open(tmp.replace(" ", "_") + '.pep.fa', 'w')
            species = tmp
            fw.write(line)
            seq = ''
        else:
            seq += line
not_in_our_db_species = plant_db_species - our_species
# print("plant_db_species ({})".format(len(plant_db_species)), plant_db_species)
# print("our_db_species ({})".format(len(our_species)), our_species)
print('not_in_our_db ({}):'.format(len(not_in_our_db_species)), not_in_our_db_species)
