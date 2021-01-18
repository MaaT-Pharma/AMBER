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
# This script retrieves NCBI ftp links from RefSeq or GenBank accession numbers.
#############

library(reutils)
library(stringr)

refSeqGenBankAccessionNb <- c("GCF_XXX.1", "GCA_XXX.1")

for (acc in refSeqGenBankAccessionNb) {
  e <- content(esummary(esearch(acc, db="assembly"), db="assembly"), "parsed")
  if ("GCF" %in% str_split(acc, "_")[[1]]) { 
    e_cds <- str_replace(e$FtpPath_Assembly_rpt, "_assembly_report.txt", "_cds_from_genomic.fna.gz")
    e_rna <- str_replace(e$FtpPath_Assembly_rpt, "_assembly_report.txt", "_rna_from_genomic.fna.gz")
    e_md5 <- paste0(dirname(e$FtpPath_Assembly_rpt),"/md5checksums.txt")
    e_genome <- str_replace(e$FtpPath_Assembly_rpt, "_assembly_report.txt", "_genomic.fna.gz")
    e_data <- as.character(c(e[c("LastMajorReleaseAccession", "Organism")],e_genome, e_cds, e_rna, e_md5))
  }
  if ("GCA" %in% str_split(acc, "_")[[1]]) {
    e_cds <- paste0(e$FtpPath_GenBank,"/",basename(e$FtpPath_GenBank), "_cds_from_genomic.fna.gz")
    e_rna <- paste0(e$FtpPath_GenBank,"/",basename(e$FtpPath_GenBank), "_rna_from_genomic.fna.gz")
    e_md5 <- paste0(e$FtpPath_GenBank,"/md5checksums.txt")
    e_genome <- paste0(e$FtpPath_GenBank,"/",basename(e$FtpPath_GenBank), "_genomic.fna.gz")
    e_data <- as.character(c(e[c("LastMajorReleaseAccession", "Organism")],e_genome, e_cds, e_rna, e_md5))
  }
  print(e_data)
}



