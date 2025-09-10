#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import numpy as np
import re

# 1. Get parameters from cirro pipeline call
ds = PreprocessDataset.from_running()
ds.logger.info("List of starting params")
ds.logger.info(ds.params)

ds.logger.info('checking ds.files')
files = ds.files
ds.logger.info(files.head())
ds.logger.info(files.columns)

# 2. Add samplesheet parameter and set equal to ds.samplesheet
ds.logger.info("Checking samplesheet parameter")
ds.logger.info(ds.samplesheet)


# we want to create a samplesheet that looks like this for the fastq entry
# group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
# PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,fastq,library_id:S1;lane:001,/path/to/PATIENT1-T_S1_L001_R1_001.fastq.gz;/path/to/PATIENT1-T_S1_L001_R2_001.fastq.gz
# PATIENT1,PATIENT1,PATIENT1-N,normal,dna,fastq,library_id:S1;lane:002,/path/to/PATIENT1-T_S1_L002_R1_001.fastq.gz;/path/to/PATIENT1-T_S1_L002_R2_001.fastq.gz
# PATIENT1,PATIENT1,PATIENT1-T-RNA,tumor,rna,fastq,library_id:S1;lane:002,/path/to/PATIENT1-T_S1_L002_R1_001.fastq.gz;/path/to/PATIENT1-T_S1_L002_R2_001.fastq.gz

# sample id comes from ds.samplesheet['sample'] but with _ replaced with .
# group id and subject id are the same and they should be whatever comes after GBM1
# sample type is normal for PBMC and tumor for everything else
# filetype is fastq
# sequence type dna if fastq filename contains "merged" -- THIS IS HARDCODED BUT NOT ALWAYS TRUE FOR OTHER RUNS
# filepath is the concatenation of fastq_1 and fastq_2 with a ; in between
# there is not library ID or lane information in ds.samplesheet so we will just use row index instead so every sample is unique

pd.set_option('display.max_columns', None)
samplesheet = pd.DataFrame(columns=["group_id","subject_id","sample_id","sample_type","sequence_type","filetype","info","filepath"])
samplesheet['filepath'] = ds.samplesheet['fastq_1'].astype(str) + ';' + ds.samplesheet['fastq_2'].astype(str)
print(samplesheet)
samplesheet['sample_id'] = [re.sub('_','.', sample_id) for sample_id in ds.samplesheet['sample'].astype(str)]
samplesheet['group_id'] = [re.split(r'[._]', sample_id)[1] for sample_id in samplesheet['sample_id']]
samplesheet['subject_id'] = samplesheet['group_id']
samplesheet['info'] = ['library_id:S' + str(i+1) + ';lane:' + str(i+1) for i in range(len(samplesheet))]
print(samplesheet)
samplesheet['sample_type'] = ['normal' if 'PBMC' in sample_id else 'tumor' for sample_id in samplesheet['sample_id']]
samplesheet['filetype'] = ['fastq' for _ in samplesheet['filepath']]
print(samplesheet)

# this run specifically -- NOT ALWAYS TRUE
samplesheet['sequence_type'] = ['dna' if 'merged' in path else 'rna' for path in samplesheet['filepath']]
print(samplesheet)

# add -RNA to sample id for rows that are rna
samplesheet.loc[samplesheet['sequence_type'] == 'rna', 'sample_id'] = samplesheet['sample_id'] + '-RNA'
print(samplesheet)

samplesheet.to_csv('samplesheet.csv', index=False)
ds.add_param("input", "samplesheet.csv")

for key in param_list:  # list() avoids modifying during iteration
    ds.remove_param(key, force=True)

ds.logger.info(ds.params)