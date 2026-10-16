#!/usr/bin/env python

import os
import sys
import re
import argparse
from Bio import SeqIO

#============================================================
parser = argparse.ArgumentParser()
parser.add_argument('--input_file', type=str, default="")
parser.add_argument('--fna_record_idx', type=int, default=1) # 1-based index
parser.add_argument('--start_pos', type=int, default=1) # 1-based index
parser.add_argument('--end_pos', type=int, default=-1) # not-including end_pos
parser.add_argument('--output_file', type=str, default="")
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()
#============================================================
input_file = args.input_file
fna_record_idx = args.fna_record_idx
start_pos = args.start_pos
end_pos = args.end_pos
output_file = args.output_file
if len(args.rest) >= 4:
    input_file = args.rest[0]
    start_pos = int(args.rest[1])
    end_pos = int(args.rest[2])
    output_file = args.rest[3]
print(f"input_file= {input_file}")
print(f"output_file= {output_file}")
print(f"start_pos= {start_pos}")
print(f"end_pos= {end_pos}")

# Define the input .fna file, start position, and end position
# input_file = "genome.fna"
# start_pos = 100000
# end_pos = 150000

# Parse the records
record_idx = 0
record = None
for i_record in SeqIO.parse(input_file, "fasta"):
    record_idx += 1
    print(i_record.id)  # Access the sequence ID
    # print(i_record.seq)  # Access the sequence itself
    # print(i_record.description)  # Access the full description line
    if record_idx == fna_record_idx:
        record = i_record
assert record is not None, f"Error [] fna_record_idx= {fna_record_idx} not found in {record_idx} records in fna"

# from IPython import embed; ns = locals().copy(); ns.update(globals()); embed(user_ns=ns);

# Read the sequence from the .fna file
# record = SeqIO.read(input_file, "fasta")

# Extract the desired subsequence
# 1-based index, 100..200 means [100,199]
if end_pos == -1:
    subsequence = record.seq[(start_pos - 1):]
else:
    subsequence = record.seq[(start_pos - 1):end_pos]

# Create a new SeqRecord object with the extracted sequence
new_record = SeqIO.SeqRecord(subsequence, id=record.id, description=record.description)

# Save the extracted sequence to a new file
# output_file = "extracted_sequence.fna"
SeqIO.write(new_record, output_file, "fasta")

print(f"[fna_cut.py] saved to {output_file}")
print("End: fna_cut.py")

"""

fna_cut.py input_file start_pos end_pos output_file

cd ~/Nov5CCIR/blast
fna_cut.py ~/Nov5CCIR/newsample_raw_ncbi_dataset_fna/GCA_002089035.2_ASM208903v2_genomic.fna 1941431 1941781 \
  _cut_1941431_1941781_GCA_002089035.2_ASM208903v2_genomic.fna

"""

