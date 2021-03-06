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
# This script creates a Gold Standard with multi-assignments and one with single assignments from gene clustering results obtained with CD-HIT. For multi-assignments, the representative sequence of each cluster is also assigned to the species of all the genes in its cluster. 
#############

library(seqinr)
library(dplyr)
library(stringr)

##############
# Input files: 
# - "SGC.fasta" and "SGC.clstr" : clustering results generated by CD-HIT
# - "taxonomyTable.tsv" : Accession Species
# - "correspGenomeGene.tsv" : correspondence between each accession number and each gene of the catalog : Accession GeneName
# For instance, this table can be generated by the following bash code parsing the compressed FASTA files of the current directory (e.g. "GC*_cds_16S_from_genomic.fna.gz" which can be obtained using merge_CDS_16S.py):
# refSeqGenBankAccessionNb="GCF_XXX.1 GCA_XXX.1" ; for i in $refSeqGenBankAccessionNb; do zcat ${i}* | grep ">" | awk -v acc=$i '{print acc "\t" substr($1,2,length($1))}' ; done
#

catalog <- read.fasta(file="SGC.fasta") 
catalog_clusters <- read.fasta(file="SGC.clstr",as.string=TRUE) 
acc2sp <- read.delim("taxonomyTable.tsv", quote="") 
acc2gene <- read.delim("correspGenomeGene.tsv", quote="",header=FALSE) 

#
############## 

catalog_names <- data.frame(Gene=getName(catalog),Length=getLength(catalog)) 
colnames(acc2gene) = c("Accession","Gene")
expectedAll <- inner_join(acc2sp, acc2gene, by="Accession") 

##############
# Generation of the "rows2add" table associating each representative gene to the species of all the genes it represents
#

seq_clstr <- grep("1\t", unlist(getSequence(catalog_clusters,as.string=T)))
annot_seq_clstr <- getAnnot(catalog_clusters[seq_clstr])

rows2add <- c()

for (seqC in 1:length(seq_clstr)) {
	
	all_seq_clstr <- unlist(strsplit(unlist(getSequence(catalog_clusters[seq_clstr[seqC]],as.string=T)),".\t"))
	rep_seq_clstr <- all_seq_clstr[grep("[*]", all_seq_clstr)]
	other_seq_clstr <- all_seq_clstr[grep("[*]", all_seq_clstr,invert=T)]
	gene_name_idx <- regexpr(">.*[.][.][.]",rep_seq_clstr)
	gene_name_rep <- str_sub(rep_seq_clstr,gene_name_idx[1]+5,gene_name_idx[1]-4+attr(gene_name_idx,"match.length"))
	gene_name_other <- lapply(other_seq_clstr[2:length(other_seq_clstr)], function(x)  str_sub(x,regexpr(">.*[.][.][.]",x)[1]+5,regexpr(">.*[.][.][.]",x)[1]-4+attr(regexpr(">.*[.][.][.]",x),"match.length")))
	acc_rep <- expectedAll[grep(paste0(gene_name_rep,"$"), expectedAll$Gene,ignore.case = T),]
	acc_other <- lapply(gene_name_other, function(x) expectedAll[grep(paste0(x,"$"), expectedAll$Gene,ignore.case = T),])

  for(oth in acc_other) {
    oth$Gene <- acc_rep$Gene
    rows2add <- rbind(rows2add, oth)
  }
}

# 
##############

##############
# Creation of the Gold Standard with multi-assignements
#

fullExpectedAll <- bind_rows(expectedAll[,c("Accession","Species","Gene")], rows2add[,c("Accession","Species","Gene")]) 
catalog_annot_species_id_full <- inner_join(catalog_names,fullExpectedAll[,c("Gene","Species")], by="Gene") 
colnames(catalog_annot_species_id_full) <- c("@@SEQUENCEID", "_LENGTH","BINID")
catalog_annot_species_id_full <- unique(catalog_annot_species_id_full) 

#
##############

##############
# Creation of the Gold Standard with single assignements
#

catalog_annot_species_id <- inner_join(catalog_names,expectedAll[,c("Gene","Species")], by="Gene") 
colnames(catalog_annot_species_id) <- c("@@SEQUENCEID", "_LENGTH","BINID")

#
##############

write.table(catalog_annot_species_id_full[,c(1,3,2)], file="GS.tsv", sep="\t", quote=F, row.names=F)
write.table(catalog_annot_species_id[,c(1,3,2)], file="GS_SA.tsv", sep="\t", quote=F, row.names=F)