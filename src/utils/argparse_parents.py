#!/usr/bin/env python

# Modifications Copyright 2020 MaaT Pharma
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


import argparse

HELP_GOLD_STANDARD_FILE = "Gold standard - ground truth - file"
HELP_FILE = "File containing precision and recall for each genome"
HELP_QUERY_FILE = "Binning file"
HELP_QUERY_FILES = "Binning files"
HELP_FASTA_FILE = "FASTA or FASTQ file with sequences of gold standard (required if gold standard file misses column _LENGTH)"
HELP_FILTER = "Filter out [FILTER]%% smallest bins (default: 0)"
HELP_GENOMES_FILE = "File with list of genomes to be removed"
HELP_LABEL = "Binning name"
HELP_LABELS = "Comma-separated binning names"
HELP_KEYWORD = "Keyword in the second column of file with list of genomes to be removed (no keyword=remove all genomes in list)"
HELP_MAP_BY_RECALL = "Map genomes to bins by maximizing completeness"
HELP_USE_SEQ_COUNTS= "Compute metrics using sequence counts instead of base pairs and map expected bins (genomes) to predicted bins by maximizing the proportion of shared sequences as compared to the expected sequence counts (if activated, -m is overrided)"
HELP_THRESHOLDS_COMPLETENESS = "Comma-separated list of min. completeness thresholds (default %%: 50,70,90)"
HELP_THRESHOLDS_CONTAMINATION = "Comma-separated list of max. contamination thresholds (default %%: 10,5)"

PARSER_GS = argparse.ArgumentParser(add_help=False)
PARSER_GS.add_argument("bin_file", help=HELP_QUERY_FILE)
PARSER_GS.add_argument("-g", "--gold_standard_file", help=HELP_GOLD_STANDARD_FILE, required=True)
PARSER_GS.add_argument("-f", "--fasta_file", help=HELP_FASTA_FILE)

PARSER_MULTI = argparse.ArgumentParser(add_help=False)
PARSER_MULTI.add_argument('file', nargs='?', type=argparse.FileType('r'), help=HELP_FILE)
PARSER_MULTI.add_argument('-p', '--filter', help=HELP_FILTER)
PARSER_MULTI.add_argument('-l', '--label', help=HELP_LABEL, required=False)

PARSER_MULTI2 = argparse.ArgumentParser(add_help=False)
PARSER_MULTI2.add_argument("bin_files", nargs='+', help=HELP_QUERY_FILES)
PARSER_MULTI2.add_argument("-g", "--gold_standard_file", help=HELP_GOLD_STANDARD_FILE, required=True)
PARSER_MULTI2.add_argument("-f", "--fasta_file", help=HELP_FASTA_FILE)
PARSER_MULTI2.add_argument('-l', '--labels', help=HELP_LABELS, required=False)
PARSER_MULTI2.add_argument('-p', '--filter', help=HELP_FILTER)
PARSER_MULTI2.add_argument('-r', '--remove_genomes', help=HELP_GENOMES_FILE, required=False)
PARSER_MULTI2.add_argument('-k', '--keyword', help=HELP_KEYWORD, required=False)
