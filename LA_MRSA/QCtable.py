#!/usr/bin/env python

#%%
import os
import sys
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
pd.set_option('display.max_rows', None)
#%%
filepath1 = "~/CCIR24/Oct16_quast_newname/transposed_report.tsv"
filepath2 = "~/CCIR24/pathogenwatch_Oct16/pathogenwatch-saureus-1nag68k8s5o7-oct16ncbinewname-species-prediction.csv"
output_prefix = "_QCtable_out"

if len(sys.argv) > 2:
    filepath1 = sys.argv[1]
    filepath2 = sys.argv[2]
if len(sys.argv) > 3:
    output_prefix = sys.argv[2]

print(f"filepath1={filepath1}")
print(f"filepath2={filepath2}")

quastRep = pd.read_csv(filepath1, sep='\t')
spID = pd.read_csv(filepath2)

#%%
# print(quastRep)
# print(spID)
#%%

# totalContig = quastRep[0]
N50 =  quastRep['N50']
# smallContig = (quastRep.iloc[0][1:] - quastRep.iloc[1][1:])
# smallContigPr = 100 * smallContig/(quastRep.iloc[0][1:])
smallContigPr = 100 * (quastRep.iloc[0:]['Total length (>= 0 bp)'] - quastRep.iloc[0:]['Total length (>= 1000 bp)'])/(quastRep.iloc[0:]['Total length (>= 0 bp)'])
mashD = spID['Mash Distance']
matchingHashesSplit = spID['Matching Hashes'].str.split('/')
matchingHashesLabel = [s[0] for s in matchingHashesSplit]
matchingHashes = spID['Matching Hashes']

totalLength = quastRep['Total length']

# contamination = pd.Series(spID[11].values, index=spID[1])

# print(contamination)
#%%
quastRep['Prop. small contigs'] = smallContigPr

df1 = quastRep
df2 = spID[['Genome Name', 'Matching Hashes', 'Mash Distance']]

merged_df2 = pd.merge(df1, df2, left_on='Assembly', right_on='Genome Name', how='left')

print(merged_df2.shape)
#%%
merged_df2.head()
#%%

merged_df2 = merged_df2.drop(columns=['Genome Name'])
#%%
merged_df2['Pass QC'] = 1
#%%
merged_df2.columns
#%%
#merged_df2.to_csv("QC_final.csv", index=False)
#%%
# df = pd.concat({'N50': N50,
#               '# contigs': totalContig,
#               'prop. small contigs': smallContigPr},axis=1)

# df.to_csv("_qualitycontrol_Oct16.csv")
#%% from slack

# N50_labels_list = N50.index.to_list() # corresponding to spID_Genome_Name
# N50_values_list = N50.to_list()
# totalContig_list = totalContig.to_list()
# smallContigPr_list = smallContigPr.to_list()
# spID_Genome_Name = spID['Genome Name'].to_list()
# spID_Matching_Hashes = spID['Matching Hashes'].to_list()

# N50_Matching_Hashes_list = []
# for i in range(len(N50_values_list)):
#     N50_Matching_Hashes_list.append( spID_Matching_Hashes [ spID_Genome_Name.index(N50_labels_list[i]) ] )

# N50_new_df = pd.DataFrame(
#     {'Genome Name': N50_labels_list,
#      'N50': N50_values_list,
#      '# of contigs': totalContig_list,
#      'Prop. of contigs smaller than 1kb': smallContigPr_list,
#      'Matching Hashes': N50_Matching_Hashes_list
#     })

# print(N50_new_df)

#                N50_Genome_Name  N50_Value N50_Matching_Hashes
# 0  SAMN17368098_Horse_CHE_2019   101763.0            985/1000
# 1  SAMN18356604_Human_MEX_2018    58609.0            995/1000
# 2  SAMEA4553050_Human_IRL_2011   252101.0            988/1000
# 3  SAMN08584184_Horse_DEU_1905   149998.0           1000/1000
# 4  SAMN17368178_Horse_CHE_2020   148226.0            976/1000

#%% # of contigs plot - not using anymore

# values, bins, bars = plt.hist(totalContig, bins=30, edgecolor="black", color="gray")
# plt.title("C", loc='left', fontsize=15)
# # plt.text(50, 170, '10KB', fontsize=6, color='red')
# plt.xlabel('# of contigs')
# plt.ylabel('Frequency')
# plt.bar_label(bars, fontsize=6)
# plt.ticklabel_format(style='plain')
# # ax.xaxis.set_major_locator(ticker.MultipleLocator(200000)
# plt.show()
#%% prop of small contigs plot

values, bins, bars = plt.hist(smallContigPr, bins=30, edgecolor="black", color="gray")
plt.title("D", loc='left', fontsize=20)
plt.text(2.6, 170, '2.5%', fontsize=8, color='red')
plt.xlabel('Proportion of genome in contigs smaller than 1kb', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.bar_label(bars, fontsize=6)
plt.ticklabel_format(style='plain')
plt.axvline(x=2.5, color='red', linestyle='--', linewidth=1)
plt.show()
#%%

values, bins, bars = plt.hist(N50, bins=30, edgecolor="black", color="gray")
plt.title("A", loc='left', fontsize=20)
# plt.text(2, 170, '10KB', fontsize=8, color='red')
plt.xlabel('N50', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.bar_label(bars, fontsize=6)
plt.axvline(x=10000, color='red', linestyle='--', linewidth=1)
plt.ticklabel_format(style='plain')
plt.show()
#%%

values, bins, bars = plt.hist(mashD, bins=30, edgecolor="black", color="gray")
plt.title("C", loc='left', fontsize=20)
#plt.text(2, 170, '10KB', fontsize=8, color='red')
plt.xlabel('Mash Distance', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.bar_label(bars, fontsize=6)
#plt.axvline(x=10000, color='red', linestyle='--', linewidth=1)
plt.ticklabel_format(style='plain')
plt.show()
#%%

values, bins, bars = plt.hist(matchingHashes, bins=15, edgecolor="black", color="gray")
#plt.title("C", loc='left', fontsize=15)
#plt.text(2, 170, '10KB', fontsize=8, color='red')
plt.xlabel('Matching Hashes')
plt.ylabel('Frequency')
plt.bar_label(bars, fontsize=6)
# plt.axvline(x=10000, color='red', linestyle='--', linewidth=1)
# plt.ticklabel_format(style='plain')
plt.show()
#%%

genSizeSTD = np.std(totalLength)
genSizeMean = np.mean(totalLength)
print(f"genSizeSTD={genSizeSTD}")
print(f"genSizeMean={genSizeMean}")

#%%
values, bins, bars = plt.hist(totalLength, bins=15, edgecolor="black", color="gray")
plt.title("B", loc='left', fontsize=20)
plt.xlabel('Genome size', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.bar_label(bars, fontsize=6)
plt.axvline(x=genSizeMean-2*genSizeSTD, color='red', linestyle='--', linewidth=1)
plt.axvline(x=genSizeMean+2*genSizeSTD, color='red', linestyle='--', linewidth=1)
plt.ticklabel_format(style='plain')
plt.tick_params(axis='both', labelsize=8)
plt.show()

#%%
for i in range(len(merged_df2['N50'])):
    # print(f"i = {i}, N50 = {merged_df2['N50'][i]}")
    if int(merged_df2['N50'][i]) < 10000:
        merged_df2['Pass QC'][i] = 0
#%%

for i in range(len(merged_df2['Prop. small contigs'])):
    # print(f"i = {i}, N50 = {merged_df2['N50'][i]}")
    if int(merged_df2['Prop. small contigs'][i]) > 2.5:
        merged_df2['Pass QC'][i] = 0
#%%

for i in range(len(merged_df2['Total length'])):
    # print(f"i = {i}, N50 = {merged_df2['N50'][i]}")
    if int(merged_df2['Total length'][i]) > genSizeMean+2*genSizeSTD or int(merged_df2['Total length'][i]) < genSizeMean-2*genSizeSTD:
        merged_df2['Pass QC'][i] = 0

#%%
merged_df2['Matching Hashes'].value_counts()
#%%
merged_df2.columns
#%%

# merged_df2.to_csv("QC_final_3.csv", index=False)
merged_df2.to_csv(f"{output_prefix}.csv", index=False)

# %%

"""

conda deactivate; conda activate py3_12_4

filepath1="~/CCIR24/Oct16_quast_newname/transposed_report.tsv"
filepath2="~/CCIR24/pathogenwatch_Oct16/pathogenwatch-saureus-1nag68k8s5o7-oct16ncbinewname-species-prediction.csv"
output_prefix="_QCtable_out"

cd ~/Summer2024/ccirwork
python ~/Summer2024/ccirwork/QCtable.py ${filepath1} ${filepath2} ${output_prefix} 2>&1 | cat | tee _QCtable.log

cat _QCtable.log | grep genSize

"""
