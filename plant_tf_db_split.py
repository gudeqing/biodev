import argparse
import sys
import json
import re

# plant_db_pep_file = sys.argv[1]
# our_annot_json = sys.argv[2]
plant_db_pep_file = "PlantTFDB-all_TF_pep.fas"
our_annot_json = "annot_species.json"
species_short_name = "species_short_name.list"

# get our species info
with open(our_annot_json) as f:
    our_annot_dict = json.load(f)
    our_species = {x.lower().replace('_', ' ') for x in our_annot_dict.keys()}

# get species short name
with open(species_short_name) as f:
    short_name_dict = dict()
    for line in f:
        if not line.strip():
            continue
        short_name, sp = line.strip().split('\t')
        short_name_dict[sp.lower()] = short_name

# get plant db species info and split fasta
match_pattern = re.compile(r'>.*?\s(.*?)\|.*')
plant_db_species = set()
with open(plant_db_pep_file) as fr:
    first_line = fr.readline()
    species = match_pattern.match(first_line).groups()[0].lower()
    plant_db_species.add(species)
    if species not in short_name_dict:
        print('{} not found in pep.list'.format(species))
        short_name = species.replace(' ', '_')
    else:
        short_name = short_name_dict[species]
    fw = open(short_name + '.pep.fa', 'w')
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
                if tmp.lower() not in short_name_dict:
                    print('{} not found in pep.list'.format(species))
                    short_name = tmp.replace(' ', '_')
                else:
                    short_name = short_name_dict[tmp]
                fw = open(short_name + '.pep.fa', 'w')
            species = tmp
            fw.write(line)
            seq = ''
        else:
            seq += line
    else:
        fw.write(seq)
not_in_our_db_species = plant_db_species - our_species
print("plant_db_species ({})".format(len(plant_db_species)), plant_db_species)
print("our_db_species ({})".format(len(our_species)), our_species)
print('not_in_our_db ({}):'.format(len(not_in_our_db_species)), not_in_our_db_species)

