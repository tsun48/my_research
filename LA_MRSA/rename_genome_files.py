#!/usr/bin/env python

# %%
import os
import sys
import re
import pandas as pd
import glob
import numpy as np
import argparse

# data = pd.read_csv(f"{os.getenv('HOME')}/CCIR24/NCBI_CC398_metadata_v3.csv")
metadata_csv = f"{os.getenv('HOME')}/Summer2024/ccirwork/NCBI_CC398_metadata_v3_Oct16.csv"
src_dir = f"{os.getenv('HOME')}/Downloads/Oct16_ncbi_dataset_onedir"
dst_dir = f"{os.getenv('HOME')}/Downloads/Oct16_ncbi_dataset_onedir_newname"

#============================================================
parser = argparse.ArgumentParser()
parser.add_argument('--metadata_csv', type=str, default=metadata_csv)
parser.add_argument('--src_dir', type=str, default=src_dir)
parser.add_argument('--dst_dir', type=str, default=dst_dir)
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()
#============================================================
metadata_csv = args.metadata_csv
src_dir = args.src_dir
dst_dir = args.dst_dir

print(f"metadata_csv= {metadata_csv}")
# data = pd.read_csv(metadata_csv, sep=",")
data = np.loadtxt(f"{metadata_csv}", delimiter=',', skiprows=1, dtype=str)

print(f"data.shape= {data.shape}")

# sys.exit(0)
#
# nameMap = {}
# for i in range(len(data["ST"])):
#     nameMap[ data["Accession_number"][i] ] = data["Linkage_name"][i]
#     print(f'nameMap[{data["Accession_number"][i]}]= {nameMap[ data["Accession_number"][i] ]}')

nameMap = {}
for i in range(data.shape[0]):
    nameMap[ data[i][0] ] = data[i][12]
    print(f'nameMap[{data[i][0]}]= {nameMap[ data[i][0] ]}')

newFileMap = {}

#%%
# src_dir = f"{os.getenv('HOME')}/CCIR24/allgenomes_fna"
os.makedirs(dst_dir, exist_ok=True)

orig_file_list = glob.glob(f"{src_dir}/*.fna")
print(f"len(orig_file_list)= {len(orig_file_list)}")

if len(orig_file_list) == 0:
    assert False, f"Error [rename_genome_files.py] no .fna files in src_dir= {src_dir}"

for file in orig_file_list:
    print(f"file= {file}")
    baseName = file.split("/")[-1]
    gca = baseName.split("_")[0]
    id = baseName.split("_")[1]
    accNum = f"{gca}_{id}"
    if accNum not in nameMap:
        assert False, f"error, accNum = {accNum} not found in nameMap"
    else:
        newFile = f"{dst_dir}/{nameMap[accNum]}.fna"
        print(f"old: {baseName}, new: {newFile}")
        if newFile not in newFileMap:
            # os.rename(file, newFile)
            os.system(f"rsync -avz {file} {newFile}")
            # faster way
            # os.system(f"/bin/mv -fv {file} {newFile}")
            newFileMap[newFile] = accNum
        else:
            assert False, f"error, accNum = {accNum}, newFile = {newFile}, newFileMap = {newFileMap[newFile]}, already exists "

#%%

'''

rename_genome_files.py --metadata_csv ~/Summer2024/ccirwork/NCBI_CC398_metadata_v3_Oct16.csv \
--src_dir ~/Downloads/Oct16_ncbi_dataset_onedir \
--dst_dir ~/Downloads/Oct16_ncbi_dataset_onedir_newname

'''

