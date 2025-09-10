#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import numpy as np
import re
pd.set_option('display.max_columns', None)

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

# dcast ds.files based on sample column, and read column indicates fastq_1 and fastq_2, use file as the content column
# this assumes that each sample has exactly one fastq_1 and one fastq_2
pivoted_file_list = ds.files.pivot(index=['sample','process'], columns='read', values='file').reset_index()
# rename column names to be more descriptive
pivoted_file_list.columns = ['sample', 'process', 'fastq_1', 'fastq_2']
pivoted_file_list['fastq_1_basename'] = pivoted_file_list['fastq_1'].apply(lambda x: x.split('/')[-1])
pivoted_file_list['fastq_2_basename'] = pivoted_file_list['fastq_2'].apply(lambda x: x.split('/')[-1])
print(pivoted_file_list)

samplesheet = pd.DataFrame(columns=["group_id","subject_id","sample_id","sample_type","sequence_type","filetype","info","filepath"])
samplesheet['filepath'] = pivoted_file_list['fastq_1'].astype(str) + ';' + pivoted_file_list['fastq_2'].astype(str)
samplesheet['sample_id'] = [re.sub('_','.', sample_id) for sample_id in pivoted_file_list['sample'].astype(str)]
samplesheet['sequence_type'] = ['dna' if 'dna' in process else 'rna' for process in pivoted_file_list['process']]


samplesheet['group_id'] = [re.split(r'[._]', sample_id)[1] for sample_id in samplesheet['sample_id']]
samplesheet['subject_id'] = samplesheet['group_id']
samplesheet['info'] = ['library_id:S' + str(i+1) + ';lane:' + str(i+1) for i in range(len(samplesheet))]
samplesheet['sample_type'] = ['normal' if 'PBMC' in sample_id else 'tumor' for sample_id in samplesheet['sample_id']]
samplesheet['filetype'] = ['fastq' for _ in samplesheet['filepath']]
samplesheet.loc[samplesheet['sequence_type'] == 'rna', 'sample_id'] = samplesheet['sample_id'] + '-RNA'
print(samplesheet)

samplesheet.to_csv('samplesheet.csv', index=False)
ds.add_param("input", "samplesheet.csv")

ds.logger.info(ds.params)