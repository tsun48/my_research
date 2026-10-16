#!/usr/bin/env python

# %%
import os
import sys
import re
import pandas as pd
import glob
import numpy as np
import argparse
sys.path.append(f"{os.getenv('HOME')}/Summer2024/ccirwork")
from mylib import map_Host_common_name_to_HostCategoty

# amr_gene_csv="~/Output_Nov5ccir_2rpt/pt04_pathgonwatch/pathogenwatch-saureus-rg0d8um0bcri-nov5_431sample-amr-genes.csv"
# speci_pred_csv="~/Output_Nov5ccir_2rpt/pt04_pathgonwatch/pathogenwatch-saureus-rg0d8um0bcri-nov5_431sample-species-prediction.csv"
# output_prefix="merged_out"

#============================================================
parser = argparse.ArgumentParser()
parser.add_argument('--amr_gene_csv', type=str, default="")
parser.add_argument('--speci_pred_csv', type=str, default="")
parser.add_argument('--sample_sccmec_type_csv', type=str, default="")
parser.add_argument('--tn916_match_csv', type=str, default="")
parser.add_argument('--output_prefix', type=str, default='merged_out')
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()
#============================================================
amr_gene_csv = args.amr_gene_csv
speci_pred_csv = args.speci_pred_csv
sample_sccmec_type_csv = args.sample_sccmec_type_csv
tn916_match_csv = args.tn916_match_csv
output_prefix = args.output_prefix

# for (k,v) in host_map.items():
#     print(f"host_map[{k}]= {v}")

print("")
print(">>> for pd.read_csv() to work, make sure run: conda deactivate; conda activate py3_12_4")
print("")

print(f"amr_gene_csv= {amr_gene_csv}")
print(f"speci_pred_csv= {speci_pred_csv}")

amr_gene_df = pd.read_csv(amr_gene_csv, sep=",")
print(f"amr_gene_df.shape= {amr_gene_df.shape}")

speci_pred_df = pd.read_csv(speci_pred_csv, sep=",")
print(f"speci_pred_df.shape= {speci_pred_df.shape}")

tn916_match_df = pd.read_csv(tn916_match_csv, sep=",")
print(f"tn916_match_df.shape= {tn916_match_df.shape}")

merged_df = amr_gene_df
merged_df['biosample_id_short'] = [ string[0:14] for string in merged_df['NAME'] ]
merged_df['Host_common_name'] = [ string.split('_')[1] for string in merged_df['NAME'] ]
merged_df['Country_Code'] = [ string.split('_')[2] for string in merged_df['NAME'] ]
merged_df['Collection_year'] = [ string.split('_')[3] for string in merged_df['NAME'] ]
merged_df['host_class'] = [ map_Host_common_name_to_HostCategoty(string.split('_')[1]) for string in merged_df['NAME'] ]
speci_pred_4col = speci_pred_df[ ["Genome Name", "Matching Hashes","p-Value","Mash Distance"] ]

merged_df = pd.merge(merged_df, speci_pred_4col, left_on='NAME', right_on='Genome Name', how='left', suffixes=('', ''))
merged_df = merged_df.drop(columns=['Genome Name'])
merged_df['pct_mhash'] = [ float(string.split('/')[0])/float(string.split('/')[1]) for string in speci_pred_df['Matching Hashes'] ]

# sccmec
sample_sccmec_type_df = pd.read_csv(sample_sccmec_type_csv, sep=",", dtype=str,
  na_filter=False, converters={'regions': lambda x: '-' if x == '' else x})
print(f"sample_sccmec_type_df.shape= {sample_sccmec_type_df.shape}")
merged_df = pd.merge(merged_df, sample_sccmec_type_df, left_on='NAME', right_on='sample', how='left', suffixes=('', ''))
merged_df = merged_df.drop(columns=['sample'])

merged_df['sccmecType'] = [ string if string != '-' else "NA" for string in merged_df['type'] ]
merged_df['sccmec'] = [ 0 if string == '-'  else 1 for string in merged_df['type'] ]
merged_df['sccmecSubType'] = merged_df['subtype']
merged_df['sccmecRegions'] = merged_df['regions']
# merged_df['sccmecRegions'] = [ f"<{string}>" if (string != '' and string != '-') else 'NA' for string in merged_df['regions'] ]
merged_df = merged_df.drop(columns=['type', 'subtype', 'regions'])

merged_df = pd.merge(merged_df, tn916_match_df, left_on='NAME', right_on='sampleId', how='left', suffixes=('', ''))
merged_df['tn916'] = [ 1 if string != '-' else 0 for string in merged_df['sseqid'] ]
merged_df = merged_df.drop(columns=['sampleId', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'pident', 'evalue'])

# WIP
# from IPython import embed; ns = locals().copy(); ns.update(globals()); embed(user_ns=ns);

merged_df.to_csv(f"{args.output_prefix}.csv", index=False)
print(f"[amr_gene_csv_merge.py] saving to {args.output_prefix}.csv")

print("End: amr_gene_csv_merge.py")

'''

conda deactivate; conda activate py3_12_4

cd ~/Output_Nov5ccir_2rpt/sccmec
cat *fna.tsv | tr ',' ';'| tr '\t' ',' | dos2unix \
  | dedup.pl -c 0| cut -d',' -f1-3,6 | sed 's/_out_//' | sed 's|\.fna||' \
   > ~/Summer2024/ccirwork/sample_sccmec_type.csv

cat ~/Summer2024/ccirwork/sample_sccmec_type.csv | cut -d',' -f2,3,4| awk '{printf "<%s>\n", $1}' |dedup.pl
1 <-,Vc,Vc>
1 <IV,-,>
1 <IV,IVa,IVa>
1 <type,subtype,regions>
4 <IV,IVd,IVd>
4 <V,Vb,Vb>
5 <multiple,-,>
29 <V,-,>
73 <IV,multiple,IVa;IVn>
161 <V,Vc,Vc>
179 <-,-,>

conda deactivate; conda activate py3_12_4
cd ~/Summer2024/ccirwork
python ~/Summer2024/ccirwork/amr_gene_csv_merge.py \
  --output_prefix amr_gene_speci_pred_merged_out \
  --amr_gene_csv   ~/Output_Nov5ccir_2rpt/pt04_pathgonwatch/pathogenwatch-saureus-rg0d8um0bcri-nov5_431sample-amr-genes.csv \
  --speci_pred_csv ~/Output_Nov5ccir_2rpt/pt04_pathgonwatch/pathogenwatch-saureus-rg0d8um0bcri-nov5_431sample-species-prediction.csv \
  --sample_sccmec_type_csv ~/Summer2024/ccirwork/sample_sccmec_type.csv \
  --tn916_match_csv ~/Summer2024/ccirwork/431sample_tn916_match.csv

cat amr_gene_speci_pred_merged_out.csv | csv_shape.pl
# 432x32 including header

cat amr_gene_speci_pred_merged_out.csv | lf.sh | less

'''

