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

# amr_gene_csv="~/Downloads/pathogenwatch-saureus-rg0d8um0bcri-nov5_431sample-amr-genes.csv"
# color_code_csv="~/Summer2024/ccirwork/amr_color_selection.csv"
# output_prefix="merged_out"

#============================================================
parser = argparse.ArgumentParser()
parser.add_argument('--amr_gene_csv', type=str, default="")
parser.add_argument('--color_code_csv', type=str, default="")
parser.add_argument('--output_prefix', type=str, default='merged_out')
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()
#============================================================
amr_gene_csv = args.amr_gene_csv
color_code_csv = args.color_code_csv
output_prefix = args.output_prefix

print("")
print(">>> for pd.read_csv() to work, make sure run: conda deactivate; conda activate py3_12_4")
print("")

print(f"amr_gene_csv= {amr_gene_csv}")
print(f"color_code_csv= {color_code_csv}")

amr_gene_df = pd.read_csv(amr_gene_csv, sep=",")
print(f"amr_gene_df.shape= {amr_gene_df.shape}")

color_code_df = pd.read_csv(color_code_csv, sep=",")
print(f"color_code_df.shape= {color_code_df.shape}")

# from IPython import embed; ns = locals().copy(); ns.update(globals()); embed(user_ns=ns);

new_df = amr_gene_df.astype( {k: 'str' for k in list(amr_gene_df.columns) } )
# insert a column in the front
new_df.insert(0, 'biosample_id_short', [ string[0:14] for string in amr_gene_df['NAME'] ])

for i in range(amr_gene_df.shape[0]):
    for fld in amr_gene_df.columns[1:]:
        new_df.loc[i, fld] = color_code_df[fld][color_code_df['AMR'] == amr_gene_df[fld][i]].item()

new_df.to_csv(f"{args.output_prefix}.csv", index=False)
print(f"[amr_gene_csv_add_color.py] saving to {args.output_prefix}.csv")

print("End: amr_gene_csv_add_color.py")

'''

conda deactivate; conda activate py3_12_4
cd ~/Summer2024/ccirwork
python ~/Summer2024/ccirwork/amr_gene_csv_add_color.py \
  --output_prefix amr_gene_color \
  --amr_gene_csv   ~/Downloads/pathogenwatch-saureus-rg0d8um0bcri-nov5_431sample-amr-genes.csv \
  --color_code_csv ~/Summer2024/ccirwork/amr_color_selection.csv

cat amr_gene_speci_pred_merged_out.csv | csv_shape.pl

'''

