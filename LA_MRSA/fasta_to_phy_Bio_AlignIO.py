#!/usr/bin/env python
import sys
import os
from Bio import AlignIO
from collections import Counter

in_fasta_file = sys.argv[1]
out_phy_file = sys.argv[2]
# in_fasta_file = "core_gene_alignment.fasta"
# out_phy_file = "output_alignment.phy"

print(f"in_fasta_file= {in_fasta_file}")
print(f"out_phy_file= {out_phy_file}")

# Read the input FASTA file
# alignment = AlignIO.read("input_alignment.fasta", "fasta")
# alignment = AlignIO.read("combined_DNA_CDS.fasta", "fasta")
# alignment = AlignIO.read("core_gene_alignment.aln", "fasta")
# alignment = AlignIO.read("core_gene_alignment.aln", "clustal")
# alignment = AlignIO.read("core_gene_alignment.aln", "clustal")
alignment = AlignIO.read(in_fasta_file, "fasta")
print(f"alignment= {alignment}")

# Get sequence names
seq_names = [record.id for record in alignment]
for index, record in enumerate(alignment):
    if seq_names.count(record.id) > 1:
        record.id = f"{record.id}_{index}"  # or any other logic to ensure uniqueness

# Count occurrences
name_counts = Counter(seq_names)
# Identify duplicates
duplicates = [name for name, count in name_counts.items() if count > 1]
print(f"duplicates= {duplicates}")

if duplicates:
    print("Duplicate sequence names found:", duplicates)
else:
    # Proceed to write the alignment if no duplicates
    AlignIO.write(alignment, out_phy_file, "phylip")

# Write the output in PHYLIP format
# AlignIO.write(alignment, "output_alignment.phy", "phylip")
# AlignIO.write(alignment, "output_alignment.phy", "phylip")
print("End: fasta_to_phy_Bio_AlignIO.py")

'''

cd ~/CCIR24/Oct16_panaroo_core
cat *core_gene_alignment.fasta | grep '>'|sed 's/>//' | awk '{print length($1)}'|dedup.pl
# 12 13
# 26 14
# 397 12

ln -s core_gene_alignment.aln core_gene_alignment.fasta

replace.pl '_PHYLIP_ID_WIDTH = 10' '_PHYLIP_ID_WIDTH = 14' ` find ~/anaconda3/envs/panaroo_py39/lib/python3.9/site-packages/Bio -name '*.py' `

cd ~/CCIR24/Oct16_panaroo_core
python ~/Summer2024/ccirwork/fasta_to_phy_Bio_AlignIO.py

clustalw2 -INFILE=core_gene_alignment_filtered.fasta -OUTFILE=output_alignment_clustalw2.phy -OUTPUT=PHYLIP

muscle -in core_gene_alignment_filtered.fasta -out output_alignment_muscle.phy -f phylip

'''
