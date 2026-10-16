#!/usr/bin/env python

# %%
import os
import sys
import re
import pandas as pd
import glob
import numpy as np
import argparse

metadata_csv = f" "
amr_gene_csv = f" "

#============================================================
parser = argparse.ArgumentParser()
parser.add_argument('--metadata_csv', type=str, default=metadata_csv)
parser.add_argument('--amr_gene_csv', type=str, default=amr_gene_csv)
parser.add_argument('--output_prefix', type=str, default='_out')
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()
#============================================================
metadata_csv = args.metadata_csv
amr_gene_csv = args.amr_gene_csv
output_prefix = args.output_prefix

host_map = {
    "africanwilddonkey": "other",
    "commondolphin": "other",
"cowwithbovinemastitis": "other",
"fallowdeer": "other",
"hedgehog": "pig",
"iberianibex": "other",
"milk": "milk",
"monkey": "other",
"plant": "other",
"turkeyvulture": "poultry",
"bird": "poultry",
"duck": "poultry",
"environment": "other",
"goat": "other",
"bovine": "bovine",
"sheep": "other",
"cow": "other",
"wildturkey": "poultry",
"rodent": "other",
"cat": "cat",
"dog": "dog",
"animal": "other",
"chicken": "poultry",
"cattle": "cattle",
"horse": "horse",
"bostaurus": "other",
"notapplicable": "other",
"domesticpig": "pig",
"pig": "pig",
"wildpig": "pig",
"human": "human",
}

for (k,v) in host_map.items():
    print(f"host_map[{k}]= {v}")

print("")
print(">>> for pd.read_csv(metadata_csv) to work, make sure run: conda deactivate; conda activate py3_12_4")
print("")

print(f"metadata_csv= {metadata_csv}")
print(f"amr_gene_csv= {amr_gene_csv}")

# metadata_df = np.loadtxt(f"{metadata_csv}", delimiter=',', skiprows=1, dtype=str)
metadata_df = pd.read_csv(metadata_csv, sep=",")
print(f"metadata_df.shape= {metadata_df.shape}")

amr_gene_df = pd.read_csv(amr_gene_csv, sep=",")
print(f"amr_gene_df.shape= {amr_gene_df.shape}")

def map_Host_common_name_to_HostCategoty(Host_common_name):
    # print(f"[map_Host_common_name_to_HostCategoty] Host_common_name= {Host_common_name}")
    res = ""
    if Host_common_name in host_map:
        res = host_map[Host_common_name]
    elif re.search("human", Host_common_name, re.IGNORECASE):
        res = "human"
    elif re.search("wild", Host_common_name, re.IGNORECASE):
        res = "wild"
    elif re.search("environment", Host_common_name, re.IGNORECASE):
        res = "environment"
    else:
        res = "misc"
    return res

Biosample_ID_map = {}
for i in range(metadata_df.shape[0]):
    Biosample_ID = metadata_df["Biosample_ID"][i]
    Biosample_ID_map[Biosample_ID] = {}
    for key in ["Country_Code", "Collection_year", "Host_common_name"]:
        Biosample_ID_map[Biosample_ID][key] = metadata_df[key][i]
        # print(f"ZZZ metadata_df[{key}][{i}]= {metadata_df[key][i]}")
        if key == "Host_common_name":
            Biosample_ID_map[Biosample_ID]["HostCategory"] = \
                map_Host_common_name_to_HostCategoty(Biosample_ID_map[Biosample_ID]["Host_common_name"])
# print(f"Biosample_ID_map= {Biosample_ID_map}")

if False:
    for Biosample_ID in Biosample_ID_map.keys():
        print(f"{Biosample_ID}", end="")
        for (k, v) in Biosample_ID_map[Biosample_ID].items():
            print(f" {k} {v}", end="")
        print("")

#------------------------------------------------------------------------------------
# freq tab by amr
#------------------------------------------------------------------------------------
amr_names = list(amr_gene_df.columns)[1:]
print(f"amr_names= {amr_names}")
cnt_by_amr = {}
for amr_name in amr_names:
    cnt_by_amr[amr_name] = {'HumanY': 0, 'HumanN': 0, 'LivestockY': 0, 'LivestockN': 0}
for i in range(amr_gene_df.shape[0]):
    NAME = amr_gene_df["NAME"][i]
    NAME_short = NAME.split('_')[0]
    assert NAME_short in Biosample_ID_map, f"Error {NAME_short} not in metadata"
    for amr_name in amr_names:
        # SAMN02391427 Country_Code NLD Collection_year 2008 Host_common_name rodent
        Host_common_name = Biosample_ID_map[NAME_short]["Host_common_name"]
        Collection_year = Biosample_ID_map[NAME_short]["Collection_year"]
        Country_Code = Biosample_ID_map[NAME_short]["Country_Code"]
        amr_val = str(amr_gene_df[amr_name][i])
        if re.search(r"human", Host_common_name, re.IGNORECASE):
            if amr_val == "1":
                cnt_by_amr[amr_name]['HumanY'] += 1
            else:
                cnt_by_amr[amr_name]['HumanN'] += 1
        else:
            if amr_val == "1":
                cnt_by_amr[amr_name]['LivestockY'] += 1
            else:
                cnt_by_amr[amr_name]['LivestockN'] += 1
# print(f"cnt_by_amr=\n{cnt_by_amr}")
freq_df = pd.DataFrame.from_records(cnt_by_amr).transpose()
# print(f"freq_df=\n{freq_df}")

# new_order = [3, 14, 15, 16, 17, 5, 6, 7, 8, 9, 10, 11, 12, 13, 2, 1, 4]  # New order of rows
new_order = ['blaZ', 'mecA', 'tetK', 'tetL', 'tetM', 'dfrA', 'dfrG', 'dfrK', 'ermA', 'ermA_SDS', 'ermB', 'ermC', 'ermC_LP', 'ermT', 'aadD', 'aacA-aphD', 'cfr', 'sccmec', 'tn916', 'chs', 'scn', 'sak', 'sea_gene', 'sep_gene']
freq_df_reordered = freq_df.reindex(new_order) 

freq_df_reordered.to_csv(f"{output_prefix}_freq_by_amr.csv", index=True)
print(f"saving to {output_prefix}_freq_by_amr.csv")
#------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------
# WIP by year/HostCategory/country
#------------------------------------------------------------------------------------
byType = "Collection_year"
for byType in ["Collection_year", "HostCategory", "Country_Code"]:
    cnt_by_key = {}
    # keys = list(sorted(metadata_df[byType].unique()))
    # print(f"keys= {keys}")
    # for key in keys:
    #     cnt_by_key[key] = {'Human': 0, 'Livestock': 0}
    for i in range(amr_gene_df.shape[0]):
        NAME = amr_gene_df["NAME"][i]
        NAME_short = NAME.split('_')[0]
        assert NAME_short in Biosample_ID_map, f"Error {NAME_short} not in metadata"
        # SAMN02391427 Country_Code NLD Collection_year 2008 Host_common_name rodent
        Host_common_name = Biosample_ID_map[NAME_short]["Host_common_name"]
        # Collection_year = Biosample_ID_map[NAME_short]["Collection_year"]
        # Country_Code = Biosample_ID_map[NAME_short]["Country_Code"]
        # print(f"i= {i} NAME_short= {NAME_short} Collection_year= {Collection_year}")
        key = Biosample_ID_map[NAME_short][byType]
        if key not in cnt_by_key:
            cnt_by_key[key] = {'Human': 0, 'Livestock': 0}
        if re.search(r"human", Host_common_name, re.IGNORECASE):
            cnt_by_key[key]['Human'] += 1
        else:
            cnt_by_key[key]['Livestock'] += 1
    freq_df = pd.DataFrame.from_records(cnt_by_key).transpose()
    freq_df.to_csv(f"{output_prefix}_freq_by_{byType}.csv", index=True)
    print(f"saving to {output_prefix}_freq_by_{byType}.csv")

#------------------------------------------------------------------------------------

print("End: amr_gene_to_freq_blast.py")

'''

conda deactivate; conda activate py3_12_4

cd ~/Output_Nov5ccir_2rpt/pt04_pathgonwatch
python ~/Summer2024/ccirwork/amr_gene_to_freq_reordered.py \
  --output_prefix nov5_431sample_reordered_alldata \
  --metadata_csv ~/Summer2024/ccirwork/NCBI_ST398_metadata_no_NA_4844row_13col.csv \
  --amr_gene_csv ~/Summer2024/ccirwork/amr_gene_speci_pred_merged_out.csv

  
cat pathogenwatch-saureus-rg0d8um0bcri-nov5_431sample-amr-genes.csv | lf.sh | grep tetM | head
cat pathogenwatch-saureus-rg0d8um0bcri-nov5_431sample-amr-genes.csv | awk -F',' '(NR==1 || $15==1){print}' > pathogenwatch-saureus-rg0d8um0bcri-nov5_431sample-amr-genes_tetM_1.csv

cat amr_gene_speci_pred_merged_out.csv | awk -F',' '(NR==1 || $15==1){print}' > amr_gene_speci_pred_merged_out_tetM_1.csv

cd ~/Output_Nov5ccir_2rpt/pt04_pathgonwatch
python ~/Summer2024/ccirwork/amr_gene_to_freq_reordered.py \
  --output_prefix nov5_431sample_reordered_alldata_LA \
  --metadata_csv ~/Summer2024/ccirwork/NCBI_ST398_metadata_no_NA_4844row_13col.csv \
  --amr_gene_csv ~/Summer2024/ccirwork/amr_gene_speci_pred_merged_out_tetM_1.csv
'''