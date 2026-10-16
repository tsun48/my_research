#!/usr/bin/env python

# %%
import os
import sys
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
pd.set_option('display.max_rows', None)

#============================================================
parser = argparse.ArgumentParser()
parser.add_argument('--genSizeMean', type=str, default='')
parser.add_argument('--genSizeSTD', type=str, default='')
parser.add_argument('--output_prefix', type=str, default='qc_by_quast')
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()
#============================================================

quast_report_tsv = args.rest[0]; # sys.argv[1]
# out_passed_names_csv = args.rest[1]; # sys.argv[2]

# quast_report_tsv = "~/CCIR24/Oct16_quast_newname/transposed_report.tsv"
print(f"quast_report_tsv= {quast_report_tsv}")

# quastRep = pd.read_csv(quast_report_csv, sep=',')
# pandas._libs.lib.maybe_convert_objects
# TypeError: Cannot convert numpy.ndarray to numpy.ndarray
# %%

# somehow this only work in py3_12_4
quastRep = pd.read_csv(quast_report_tsv, sep='\t')
N50 =  quastRep['N50']
smallContigPr = 100 * (quastRep.iloc[0:]['Total length (>= 0 bp)'] - quastRep.iloc[0:]['Total length (>= 1000 bp)'])/(quastRep.iloc[0:]['Total length (>= 0 bp)'])
quastRep['Prop. small contigs'] = smallContigPr
totalLength = quastRep['Total length']

if args.genSizeMean != '' and args.genSizeSTD != '':
    genSizeSTD = float(args.genSizeSTD)
    genSizeMean = float(args.genSizeMean)
    print(f"overriding genSizeSTD= {genSizeSTD}")
    print(f"overriding genSizeMean= {genSizeMean}")
else:
    genSizeSTD = np.std(totalLength)
    genSizeMean = np.mean(totalLength)
    print(f"creating genSizeSTD= {genSizeSTD}")
    print(f"creating genSizeMean= {genSizeMean}")

merged_df2 = quastRep

merged_df2['Pass QC'] = 1
merged_df2['FailReason'] = 'Passed'

#%%
for i in range(len(merged_df2['N50'])):
    # print(f"i = {i}, N50 = {merged_df2['N50'][i]}")
    if int(merged_df2['N50'][i]) < 10000:
        merged_df2['Pass QC'][i] = 0
        merged_df2['FailReason'][i] = 'N50_lt_10000'
#%%

for i in range(len(merged_df2['Prop. small contigs'])):
    # print(f"i = {i}, N50 = {merged_df2['N50'][i]}")
    if int(merged_df2['Prop. small contigs'][i]) > 2.5:
        merged_df2['Pass QC'][i] = 0
        merged_df2['FailReason'][i] = 'prop_small_contigs_gt_2.5'
#%%

for i in range(len(merged_df2['Total length'])):
    # print(f"i = {i}, N50 = {merged_df2['N50'][i]}")
    if int(merged_df2['Total length'][i]) > genSizeMean+2*genSizeSTD or int(merged_df2['Total length'][i]) < genSizeMean-2*genSizeSTD:
        merged_df2['Pass QC'][i] = 0
        merged_df2['FailReason'][i] = 'total_length_exceed_2std'

#%%
print(merged_df2.columns)
print(merged_df2[merged_df2['Pass QC']==0])
# merged_df2.to_csv('_merged_df2.csv', index=False)
merged_df2.to_csv(f"{args.output_prefix}_pass_reason.csv", index=False)

QC_final_passed_names_df = pd.DataFrame.from_dict({'Assembly': merged_df2['Assembly'][merged_df2['Pass QC']==1]})

# QC_final_passed_names_df.to_csv("QC_final_passed_names.csv", index=False, header=False)
# QC_final_passed_names_df.to_csv(out_passed_names_csv, index=False, header=False)
QC_final_passed_names_df.to_csv(f"{args.output_prefix}_passed_names.csv", index=False, header=False)

print("End: QCByQuast.py")

'''

conda_deactivate; conda activate py3_12_4
cd ~/Summer2024/ccirwork
QCByQuast.py ~/CCIR24/Oct16_quast_newname/transposed_report.tsv QC_final_passed_names.csv
wc QC_final_passed_names.csv
  435   435 12480 QC_final_passed_names.csv

cd ~/Nov5CCIR/nov5_as_start_458sample_quast
~/Summer2024/ccirwork/QCByQuast.py ~/Nov5CCIR/nov5_as_start_458sample_quast/transposed_report.tsv \
  nov5_as_start_458sample_QC_passed_names.csv
wc nov5_as_start_458sample_QC_passed_names.csv
  431   431 12403 nov5_as_start_458sample_QC_passed_names.csv

# with fixed mean/std of total length
conda deactivate; conda activate py3_12_4
cd ~/Nov5CCIR/nov5_as_start_458sample_quast
~/Summer2024/ccirwork/QCByQuast.py \
  --genSizeMean 2820883 \
  --genSizeSTD 88396 \
  --output_prefix nov5_as_start_458sample_with_mean_std_qc \
  ~/Nov5CCIR/nov5_as_start_458sample_quast/transposed_report.tsv \
  2>&1 | cat | tee _nov5_as_start_458sample_with_mean_std_qc.log
cat _nov5_as_start_458sample_with_mean_std_qc.log | grep gen
cat nov5_as_start_458sample_with_mean_std_qc_pass_reason.csv| cut -d',' -f26 | dedup.pl

'''

