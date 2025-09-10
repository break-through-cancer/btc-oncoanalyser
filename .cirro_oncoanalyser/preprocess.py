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

param_list = ["dna_cram_path","dna_cram_path_normal","rna_bam_path"]
# DNA cram files are like this -- data/preprocessing/recalibrated/GBM1.DFCI4.S1.C4/GBM1.DFCI4.S1.C4.recal.cram
# DNA normal cram are like this -- data/preprocessing/recalibrated/GBM1.DFCI4.PBMC/GBM1.DFCI4.PBMC.recal.cram
# rna bam files are like this -- data/star_salmon/GBM1_DFCI4_S1_C4.markdup.sorted.bam

# we want to create a samplesheet that looks like this:
# group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
# PATIENT1,PATIENT1,PATIENT1-N,normal,dna,bam,/path/to/PATIENT1-N.dna.bam
# PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
# PATIENT1,PATIENT1,PATIENT1-T-RNA,tumor,rna,bam,/path/to/PATIENT1-T.rna.bam

# group id and subject id are the same and they should be whatever comes after GBM1
# sample id is what comes after subject id
# sample type is normal for PBMC and tumor for everything else
# sequence type is dna for cram files and rna for bam files
# filetype is cram for cram files and bam for bam files
# filepath is the full path to the file which we canb get by doing ds.params.get('dna_cram_path') or ds.params.get('rna_bam_path')

samplesheet = pd.DataFrame(columns=["group_id","subject_id","sample_id","sample_type","sequence_type","filetype","filepath"])
samplesheet['filepath'] = np.array([ds.params.get(param) for param in param_list]).flatten()
print(np.array([ds.params.get(param) for param in param_list]).flatten())
print(samplesheet)

samplesheet['sample_type'] = ['normal' if 'PBMC' in path else 'tumor' for path in samplesheet['filepath']]
samplesheet['group_id'] = [re.split(r'[._]', path.split('/')[-1])[1] for path in samplesheet['filepath']]
samplesheet['subject_id'] = samplesheet['group_id']
samplesheet['sample_id'] = [re.split(r'[._]', path.split('/')[-1])[0] for path in samplesheet['filepath']]
print(samplesheet)

# this run specifically has cram for dna and bam for rna
samplesheet['filetype'] = ['cram' if 'cram' in path else 'bam' for path in samplesheet['filepath']]
samplesheet['sequence_type'] = ['dna' if 'cram' in path else 'rna' for path in samplesheet['filepath']]
print(samplesheet)

samplesheet.to_csv('samplesheet.csv', index=False)
ds.add_param("input", "samplesheet.csv")

for key in param_list:  # list() avoids modifying during iteration
    ds.remove_param(key, force=True)

ds.logger.info(ds.params)