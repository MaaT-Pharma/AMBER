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
# This script copies the content of each file named "GC*_cds_from_genomic.fna.gz" in a folder named "cds", extracts the 16S rDNA sequences from the corresponding "GC*_rna_from_genomic.fna.gz" in a folder named "rna" located at the same path as the "cds" directory, and copies the 16S rDNA sequences in each newly copied cds file.
# python merge_CDS_16S.py /path/to/cds/directory/cds /path/to/output/directory/cds_16S 
#############

from Bio import SeqIO
import os, sys, gzip

if len(sys.argv) == 3 : 
	cds_dir = sys.argv[1]
	output_dir = sys.argv[2]


	if os.path.isdir(cds_dir) : 

		cds_dir_files = os.listdir(cds_dir)

		for cds_file in cds_dir_files : 
			cds=cds_dir+"/"+cds_file
			if os.path.isfile(cds) and "GC" in os.path.basename(cds) :

				rna_file = cds.replace("cds", "rna")

				output_path = output_dir + "/" + cds_file.replace("cds", "cds_16S")
				print(cds)
				output_file = gzip.open(output_path, "a") 
				
				with gzip.open(cds, 'r') as gfp:
					for line in gfp : 
						output_file.write(line)
	
				if os.path.isfile(rna_file):
					with gzip.open(rna_file, 'r') as rfp:
						for seq_record in SeqIO.parse(rfp, "fasta"):
							if "[product=16S ribosomal RNA]" in seq_record.description : 
								SeqIO.write(seq_record, output_file, "fasta")

				output_file.close()
    