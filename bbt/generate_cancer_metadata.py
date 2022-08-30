import os
import pandas as pd
import shortuuid
import random


meta_dict = {
    'File name': {'type': str, 'range': [], 'desc': '文件路径，如果/sample1/r1.fq.gz'},
    'Primary_site': {'type': str, 'range': ['Colorectal', 'Pancreas', 'Stomach', 'Liver', 'Breast', 'Lung'], 'desc': '原发灶位置'},
    'Disease_type': {'type': str, 'range': ['Colon Adenocarcinoma', 'Rectum Adenocarcinoma', 'Pancreatic Adenocarcinoma', 'Stomach Adenocarcinoma',
                                            'Liver Hepatocellular Carcinom', 'Breast Invasive Carcinoma',
                                            'Lung Adenocarcinoma', 'Lung Squamous Cell Carcinoma'], 'desc': ''},
    'Case_id': {'type': str, 'range': [], 'desc': ''},
    'Age': {'type': int, 'range':list(range(1,101)), 'desc': ''},
    'Vital_status': {'type': str, 'range': ['Alive', 'Dead', 'Unknown'], 'desc': ''},
    'Race': {'type': str, 'range': ['white', 'american indian or alaska native', 'black or african american', 'asian', 'native hawaiian or other pacific islander', 'Unknown', 'other'], 'desc': ''},
    'Gender': {'type': str, 'range': ['Female', 'Male', 'Unknown'], 'desc': ''},
    'Ethnicity': {'type': str, 'range': ['Unknown'], 'desc': ''},
    'Primary_diagnosis': {'type': str, 'range': ['Unknown'], 'desc': ''},
    'Age_at_diagnosis': {'type': str, 'range':list(range(1,101)), 'desc': ''},
    'Tissue_type': {'type': str, 'range': ['Tumor', 'Normal', 'Abnormal', 'Unknown'], 'desc': ''},
    'Sample_type': {'type': str, 'range': ['Solid Tissue', 'Plasma', '3D Organoid', 'Serum', 'Peripheral Whole Blood', 'cell'], 'desc': ''},
    'Sample_id': {'type': str, 'range': [], 'desc': ''},
    'Data_category': {'type': str, 'range': ['Copy Number Variation', 'Simple Nucleotide Variation', 'Transcriptome Profiling', 'DNA Methylation', ], 'desc': ''},
    'Reference_genome': {'type': str, 'range': ['GRCh37/hg19', 'GRCh37/hg19', ''], 'desc': ''},
    'Quality_scale': {'type': str, 'range': ['illumina 1.8'], 'desc': ''},
    'Data_format': {'type': str, 'range': ['FASTQ.GZ', 'VCF'], 'desc': ''},
    'Platform': {'type': str, 'range': ['Illumina'], 'desc': ''},
    'Method of detection': {'type': str, 'range': ['WXS', 'RNA-Seq', 'Bisulfite-Seq', 'Targeted Sequencing'], 'desc': ''}
}

disease_to_primary_site = {
    'Colon Adenocarcinoma': 'Colorectal',
    'Rectum Adenocarcinoma': 'Colorectal',
    'Pancreatic Adenocarcinoma': 'Pancreas',
    'Stomach Adenocarcinoma': 'Stomach',
    'Liver Hepatocellular Carcinom': 'Liver',
    'Breast Invasive Carcinoma': 'Breast',
    'Lung Adenocarcinoma': 'Lung',
    'Lung Squamous Cell Carcinoma': 'Lung'
}
outdir = 'dataset'
os.makedirs(outdir)
su = shortuuid.ShortUUID(alphabet="0123456789")
rows = []
for disease in meta_dict['Disease_type']['range']:
    for i in range(random.randint(30, 80)):
        row = []
        ref = meta_dict['Reference_genome']['range'][random.randint(0, 2)]
        category = meta_dict['Data_category']['range'][random.randint(0, 3)]
        sample_id = 'S' + str(su.random(length=6))
        if ref and category=="Simple Nucleotide Variation":
            row.append(f'/{outdir}/{sample_id}.vcf')
        else:
            ref = ''
            row.append(f'/{outdir}/{sample_id}.r1.fq.gz')
        with open(f'{row[0][1:]}', 'w') as f:
            f.write('this is empty file for category function test\n')
        row.append(disease_to_primary_site[disease])
        row.append(disease)
        # case id
        row.append(su.random(length=7))
        age = random.randint(45, 77)
        row.append(age)
        row.append(meta_dict['Vital_status']['range'][random.randint(0,2)])
        row.append(['white', 'asian'][random.randint(0, 1)])
        row.append(meta_dict['Gender']['range'][random.randint(0,2)])
        # Ethnicity
        row.append('Unknown')
        # Primary_diagnosis
        row.append('Unknown')
        # Age_at_diagnosis
        row.append(age)
        row.append(meta_dict['Tissue_type']['range'][random.randint(0, 2)])
        row.append(meta_dict['Sample_type']['range'][random.randint(0, 2)])
        row.append(sample_id)
        row.append(category)
        row.append(ref)
        row.append('illumina 1.8')
        if ref and category == 'Simple Nucleotide Variation':
            file_format = 'VCF'
        else:
            file_format = 'FASTQ.GZ'
        row.append(file_format)
        row.append(meta_dict['Platform']['range'][0])
        if file_format == 'VCF':
            row.append('WXS')
        else:
            row.append(meta_dict['Method of detection']['range'][random.randint(0, 3)])
        # add row of read1
        rows.append(row.copy())
        if not ref:
            # add row of read2
            row[0] = f'/{outdir}/{sample_id}.r2.fq.gz'
            rows.append(row)
            with open(f'{row[0][1:]}', 'w') as f:
                f.write('this is empty file for category function test\n')
df = pd.DataFrame(rows, columns=list(meta_dict.keys()))
df.to_excel(f'{outdir}/Demo.metadata.xlsx', index=False)
