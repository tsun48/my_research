#!/usr/bin/env python

# %%
import os
import sys
import re
import pandas as pd
import glob
import numpy as np
import argparse
from Bio import AlignIO
from collections import Counter

#SKIP sys.path.append(f"{os.getenv('HOME')}/Summer2024/ccirwork")
#SKIP from mylib import map_Host_common_name_to_HostCategoty
#SKIP
#SKIP in_fasta_file = "core_gene_alignment.fasta"
#SKIP out_phy_file = "_tmpout_alignment.phy"
#SKIP
#SKIP alignment = AlignIO.read(in_fasta_file, "fasta")
#SKIP # <<class 'Bio.Align.MultipleSeqAlignment'> instance (432 records of length 2021720) at 147d4e35c4c0>
#SKIP len(alignment)
#SKIP rec = alignment[0]
#SKIP type(rec)
#SKIP # Bio.SeqRecord.SeqRecord
#SKIP rec.seq
#SKIP len(str(rec.seq))
#SKIP rec.id
#SKIP rec.name

# def concatenate_core_genome_alignments(core_names, output_dir, hc_threshold):
from joblib import Parallel, delayed
import networkx as nx
from collections import defaultdict
import numpy as np
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools as iter
from tqdm import tqdm
import random

def update_col_counts(col_counts, s):
    s = np.array(bytearray(s.lower().encode()), dtype=np.int8)
    s[(s!=97) & (s!=99) & (s!=103) & (s!=116)] = 110
    col_counts[0,s==97] += 1
    col_counts[1,s==99] += 1
    col_counts[2,s==103] += 1
    col_counts[3,s==116] += 1
    col_counts[4,s==110] += 1
    return (col_counts)

def calc_hc(col_counts):
    with np.errstate(divide='ignore', invalid='ignore'):
        col_counts = col_counts/np.sum(col_counts,0)
        hc = -np.nansum(col_counts[0:4,:]*np.log(col_counts[0:4,:]), 0)
    return(np.sum((1-col_counts[4,:]) * hc)/np.sum(1-col_counts[4,:]))

# ~/anaconda3/envs/panaroo_py39/lib/python3.9/site-packages/panaroo/generate_output.py
# def concatenate_core_genome_alignments(core_names, output_dir, hc_threshold):

# panaroo_dir = f"{os.getenv('HOME')}/Nov5CCIR/Nov13_431sample_plus_1outgroup_panaroo"
# included_genes_txt_file = f"{os.getenv('HOME')}/Summer2024/ccirwork/included_genes.txt"
# output_prefix = "included"

if len(sys.argv) > 3:
    panaroo_dir = sys.argv[1]
    included_genes_txt_file = sys.argv[2]
    output_prefix = sys.argv[3]
else:
    print(f"Usage: included_core_gene_to_aln_and_phy.py panaroo_dir included_genes_txt_file output_prefix")
    sys.exit(1)

# output_dir = f"{os.getenv('HOME')}/Nov5CCIR/Nov13_431sample_plus_1outgroup_panaroo"
output_dir = panaroo_dir
alignments_dir = output_dir + "/aligned_gene_sequences/"
alignment_filenames = os.listdir(alignments_dir)
core_names = [ s.rstrip() for s in open(included_genes_txt_file, "r").readlines() ]
print(f"len(core_names)= {len(core_names)}")
core_filenames = [ x for x in alignment_filenames if x.split('.')[0] in core_names ]
# Read in all these alignments
gene_alignments = []
isolates = set()
for idx, filename in enumerate(core_filenames):
    print(f"reading idx= {idx} filename= <{filename}>")
    gene_name = os.path.splitext(os.path.basename(filename))[0]
    alignment = AlignIO.read(alignments_dir + filename, "fasta")
    gene_dict = {}
    for record in alignment:
        if len(gene_dict)<1:
            gene_length = len(record.seq)
            col_counts = np.zeros((5,gene_length), dtype=float)
        col_counts = update_col_counts(col_counts, str(record.seq))
        if record.id[:3] == "_R_":
            record.id = record.id[3:]
        genome_id = record.id.split(";")[0]
        if genome_id in gene_dict:
            if str(record.seq).count("-") < str(gene_dict[genome_id][1]).count("-"):
                gene_dict[genome_id] = (record.id, record.seq)
        else:
            gene_dict[genome_id] = (record.id, record.seq)
        isolates.add(genome_id)
    gene_alignments.append((gene_name, gene_dict, gene_length, calc_hc(col_counts)))
print(f"len(gene_alignments)= {len(gene_alignments)}")
# Combine them
isolate_aln = []
for idx, iso in enumerate(isolates):
    print(f"idx= {idx} iso= {iso}")
    seq = ""
    for gene in gene_alignments:
        if iso in gene[1]:
            seq += gene[1][iso][1]
        else:
            seq += "-" * gene[2]
    isolate_aln.append(SeqRecord(seq, id=iso, description=""))
print(f"len(isolate_aln)= {len(isolate_aln)}")

# Write out the two output files
n_gene = len(gene_alignments)
n_iso = len(isolate_aln)
SeqIO.write(isolate_aln, output_dir + f"/{output_prefix}_iso{n_iso}_gene{n_gene}_alignment.aln", "fasta")
print(f"save to " + output_dir + f"/{output_prefix}_iso{n_iso}_gene{n_gene}_alignment.aln")

# Convert to .phy file
alignment = AlignIO.read(output_dir + f"/{output_prefix}_iso{n_iso}_gene{n_gene}_alignment.aln", "fasta")
AlignIO.write(alignment, output_dir + f"/{output_prefix}_iso{n_iso}_gene{n_gene}_alignment.phy", "phylip")
print(f"save to " + output_dir + f"/{output_prefix}_iso{n_iso}_gene{n_gene}_alignment.phy")
# write_alignment_header(gene_alignments, output_dir, "core_alignment_header.embl")

print("End: included_core_gene_to_aln_and_phy.py")

"""

# included_core_gene_to_aln_and_phy.py panaroo_dir included_genes_txt_file output_prefix

cd ~/Nov13_431sample_plus_1outgroup_panaroo
python ~/Summer2024/ccirwork/included_core_gene_to_aln_and_phy.py \
  ~/Nov5CCIR/Nov13_431sample_plus_1outgroup_panaroo \
  ~/Summer2024/ccirwork/included_genes.txt \
  "contour_plot_included_genes"

mkdir -p ~/Nov5CCIR/contour_plot_included_genes_iso432_gene721_alignment_raxml_ml1
cd       ~/Nov5CCIR/contour_plot_included_genes_iso432_gene721_alignment_raxml_ml1
ln -s ../Nov13_431sample_plus_1outgroup_panaroo/ contour_plot_included_genes_iso432_gene721_alignment.phy
ncpu=` cpu_count.py | awk '{print int(.9*$1)}' `; echo "ncpu=${ncpu}"
time raxmlHPC-PTHREADS -T ${ncpu} -s contour_plot_included_genes_iso432_gene721_alignment.phy \
     -# 1 -m GTRGAMMA -p 12345 \
     -n contour_plot_included_genes_iso432_gene721_alignment_raxml_ml1 \
     2>&1 | cat | tee _contour_plot_included_genes_iso432_gene721_alignment_raxml_ml1.log

cd ~/Nov5CCIR
zip_dir.sh contour_plot_included_genes_iso432_gene721_alignment_raxml_ml1

"""

