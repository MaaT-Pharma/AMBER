#!/usr/bin/env python

# Copyright 2020 MaaT Pharma
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#############
# This script extracts the sequences with a length greater than or equal to a length threshold from a FASTA file.
# python filter_FASTA_by_seq_length.py in.fasta out.fasta 1000
#############

from Bio import SeqIO
import sys, os 

if len(sys.argv) == 4 : 
	fasta_file = sys.argv[1]
	output_file = sys.argv[2]
	length = int(sys.argv[3])
	
	output_file = open(output_file, "w") 

	if os.path.isfile(fasta_file) : 
		with open(fasta_file, 'r') as ff:
			for seq_record in SeqIO.parse(ff, "fasta"):
				seq_length = len(seq_record.seq) - seq_record.seq.count("N")
				if (seq_length >= length) : 
					SeqIO.write(seq_record, output_file, "fasta")


	output_file.close()