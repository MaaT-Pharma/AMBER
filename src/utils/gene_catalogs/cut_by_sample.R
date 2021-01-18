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
# This script cuts a large count table by column to create sample-specific files (with low RAM usage). 
# Before running it, the lines of the count table must be split into several files. It must then be launched sequentially on all the files.
# For instance (bash):
# split -l 20583 --numeric-suffixes counts.tsv counts_split_
# nb_samples=1145
# for i in `ls counts_split_*`; do Rscript cut_by_sample.R $i $nb_samples; done
#
# The "counts.tsv" file must contain the following tab-separated columns: geneID <SAMPLE_ID_1> <SAMPLE_ID_2> ... <SAMPLE_ID_nb_samples>.
###############

library(data.table)

args <- commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Two argument must be supplied (input file and number of samples)", call.=FALSE)
}

output_prefix <- "counts_"
nb_samples <- as.numeric(args[2])

col_samples <- nb_samples + 1

##############
# Copies the content of the current file into as many files as the number of samples 
#

for (currSample in 2:col_samples) {

    if (currSample > 2) {
      rm(abundance_profile, output_file, sample_file)
    } 

    abundance_profile <- fread(args[1], select=c(1,currSample), nThread=16)   

    if (currSample-1 < 10) {
      sample_file <- paste0("0",currSample-1)
    } else {
      sample_file <- currSample-1
    }
    
    output_file <- paste0(output_prefix, sample_file,".tsv") 
                        
    fwrite(abundance_profile, file=output_file, sep="\t",row.names=F,quote=F,append=T,nThread=16,col.names=F)
}

# 
###############
