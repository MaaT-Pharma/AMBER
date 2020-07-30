# Scope

This a fork from the [AMBER repository](https://github.com/CAMI-challenge/AMBER) (Assessment of Metagenome BinnERs) published by [Meyer *et al.* 2018](https://doi.org/10.1093/gigascience/giy069) and developed in the context of the [Critical Assessment of Metagenomic Interpretation](http://www.cami-challenge.org/) :

* Fernando Meyer, Peter Hofmann, Peter Belmann, Ruben Garrido-Oter, Adrian Fritz, Alexander Sczyrba, and Alice C. McHardy. (2018). **AMBER: Assessment of Metagenome BinnERs.** *GigaScience*, giy069. doi:[10.1093/gigascience/giy069](https://doi.org/10.1093/gigascience/giy069)

The metrics implemented in AMBER were used and described in the CAMI manuscript, thus you may also cite:
* Sczyrba, Hofmann, Belmann, et al. (2017). **Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software.** *Nature Methods*, 14, 11:1063–1071. doi:[10.1038/nmeth.4458](https://doi.org/10.1038/nmeth.4458)

This version handles overlapping binning results (covers). It computes adapted metrics (*e.g.* (A)RI and accuracy replaced by the 
Generalized Normalized Mutual Information - GNMI) and allows to map expected bins to predicted bins by their proportion of shared sequences as compared to the number of expected sequences. Please note that this is a lighter version than the original one, therefore only the main AMBER script is supported. As described in the last section, according to the selected options, additional files are generated (such as an interactive HTML version of individual heatmaps or tables  to reproduce figures). 

# Requirements

See [default.txt](requirements/default.txt) for all python dependencies.

To compute the GNMI metric : 

* Install the standalone program of GenConvMI implementation from the [GenConvMI](https://github.com/eXascaleInfolab/GenConvMI) repository.
* Place the gecmi executable in /ports.
* Convert each [Bioboxes](https://github.com/bioboxes/rfc/tree/master/data-format) input files into the [CNL](https://github.com/eXascaleInfolab/GenConvMI) format. Make sure that all unassigned sequences are represented as singleton bins in the CNL file (one line per unassigned sequence ID). 
* Ensure that each pair of Biobox and CNL files have the same basename and the ".binning" and ".cnl" extensions (*e.g.* binning_results_1.binning and binning_results_1.cnl).
* Place each CNL file in the same directory as its corresponding Biobox file.

An example of conversion R script is provided in /src/utils/convert_to_biobox_cnl_formats.R. 

# User Guide

## amber.py

~~~BASH
usage: amber.py [-h] -g GOLD_STANDARD_FILE [-f FASTA_FILE] [-l LABELS]
                [-p FILTER] [-r REMOVE_GENOMES] [-k KEYWORD] -o OUTPUT_DIR
                [-m] [-x MIN_COMPLETENESS] [-y MAX_CONTAMINATION] [-c]
                [-i] [-s] [-v]
                bin_files [bin_files ...]

Compute all metrics and figures for one or more binning files; output summary
to screen and results per binning file to chosen directory

positional arguments:
  bin_files             Binning files

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file with sequences of gold standard
                        (required if gold standard file misses column _LENGTH)
  -l LABELS, --labels LABELS
                        Comma-separated binning names
  -p FILTER, --filter FILTER
                        Filter out [FILTER]% smallest bins (default: 0)
  -r REMOVE_GENOMES, --remove_genomes REMOVE_GENOMES
                        File with list of genomes to be removed
  -k KEYWORD, --keyword KEYWORD
                        Keyword in the second column of file with list of
                        genomes to be removed (no keyword=remove all genomes
                        in list)
  -n MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum length of sequences
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to write the results to
  -m, --map_by_completeness
                        Map genomes to bins by maximizing completeness
  -x MIN_COMPLETENESS, --min_completeness MIN_COMPLETENESS
                        Comma-separated list of min. completeness thresholds
                        (default %: 50,70,90)
  -y MAX_CONTAMINATION, --max_contamination MAX_CONTAMINATION
                        Comma-separated list of max. contamination thresholds
                        (default %: 10,5)
  -c, --plot_heatmaps   Plot heatmaps of confusion matrices (can take some
                        minutes)
  -i, --covers          Compute GNMI (compatible with standard NMI values) and percentage of sequences or base pairs assigned at least to one predicted bin (unique), instead of ARI/accuracy
  -s, --use_seq_counts
                        Compute metrics using sequence counts instead of base pairs and map expected bins (genomes) to predicted bins by maximizing the proportion of shared sequences as compared to the expected sequence counts (if activated, -m is overrided)
  -v, --version         show program's version number and exit
~~~

## Example: 

~~~BASH
python amber.py -g gold_standard.binning \
-o output_dir \
binning_results_1.binning binning_results_2.binning \
-s -i -c \
-l "binning_results_1, binning_results_2" 
~~~

## List of output files : 

The main output directory will contain : 
* **summary.tsv**: contains a summary of the computed metrics described in the section below
* **avg_purity_completeness_per_bin.png + .pdf + .eps**: figure of average purity per bin vs. average completeness per bin 
* **avg_purity_completeness.png  + .pdf**: figure of purity vs. completeness per base pair (default) or sequence (-s option)
* **ari_vs_assigned.png + .pdf**: figure of adjusted Rand index weighed by number of base pairs or sequences counts (-s option) vs. percentage of assigned base pairs or sequences (not computed if -i is activated)
* **purity_completeness_per_bin**: figure of purity per bin vs completeness per expected bin
* **rankings.txt**: tools sorted by average purity per bin, average completeness per bin, and sum of average purity per bin and average completeness per bin 
* **boxplot_purity.png + .pdf + .eps**: boxplot of purity per bin
* **boxplot_completeness.png + .pdf + .eps**:  boxplot of completeness per bin, including unassigned expected bins (not mapped to any predicted bin) associated to a completeness of 0
* **boxplot_completeness_mapped.png + .pdf + .eps**: boxplot of completeness per bin, excluding unassigned expected bins (not mapped to any predicted bin)
* **df_boxplot_completeness_unasssignedFalse.csv**: table used to generate boxplot_completeness_mapped.png
* **df_boxplot_completeess_unassignedTrue.csv**: table used to generate boxplot_completeness.png
* **df_boxplot_purity_unassignedTrue.csv**:  table used to generate boxplot_purity.png
* **GNMI_max.tsv**: GNMI value for each pair of binning results

In the same directory, a subdirectory for each input binning results will be created with the following files:
* **rand_index.tsv**: contains value of (adjusted) Rand index and percentage of assigned bases or sequences. Rand index is both weighed by base pairs or sequences counts (not computed if -i is activated)
* **percentage_assigned.tsv** : contains percentage of (unique) assigned sequences (-s) or base pairs (default) (computed only if -i is activated)
* **purity_completeness.tsv**: contains purity per bin and completeness values per bin
* **purity_completeness_avg.tsv**: contains purity and completeness values averaged over expected bins. Includes standard deviation and standard error of the mean
* **purity_completeness_by_bpcount.tsv or purity_completeness_by_seqcount.tsv**: contains purity and completeness weighed by base pairs (default) or sequences counts (-s option)
* **heatmap.png + .pdf + .eps**: heatmap representing either (-s option) the proportion of shared sequences between the set of predicted bins and the set of expected bins as compared to the number of expected sequences; or (default) base pair assignments to predicted bins vs. their true origins from the underlying expected bin
* **heatmap_2.html + .png + .eps + .pkl**: interactive HTML version of heatmap representing the proportion of shared sequences between the set of predicted bins and the set of expected bins as compared to the number of expected sequences (only generated with the -s option)
* **df_confusion.tsv**: table used to create heatmaps


## List of metrics and abbreviations

* **avg_purity_per_bin**: purity averaged over expected bins
* **std_dev_purity_per_bin**: standard deviation of purity averaged over expected bins
* **sem_purity_per_bin**: standard error of the mean of purity averaged over expected bins
* **avg_completeness_per_bin**: completeness averaged over expected bins, including unassigned expected bins associated to a completeness of 0
* **std_dev_completeness_per_bin**: standard deviation of completeness averaged over expected bins
* **sem_completeness_per_bin**: standard error of the mean of completeness averaged over expected bins
* **avg_purity**: average purity per base pair (default) or per sequence (-s option)
* **avg_completeness**: average completeness per base pair (default) or per sequence (-s option)
* **rand_index_by_bp**: Rand index weighed by base pairs (not computed if -i is activated)
* **rand_index_by_seq**: Rand index weighed by sequence counts (not computed if -i is activated)
* **a_rand_index_by_bp**: adjusted Rand index weighed by base pairs (not computed if -i is activated)
* **a_rand_index_by_seq**: adjusted Rand index weighed by sequence counts (not computed if -i is activated)
* **percent_assigned**: percentage of sequences (-s option) or base pairs (default) that were assigned to bins
* **percent_assigned_unique**: percentage of unique sequences (-s option) or base pairs (default) that were assigned to bins (computed only if -i is activated)
* **accuracy**: accuracy (not computed if -i is activated) 
* **\>0.5compl<0.1cont**: number of bins with more than 50% completeness and less than 10% contamination
* **\>0.7compl<0.1cont**: number of bins with more than 70% completeness and less than 10% contamination
* **\>0.9compl<0.1cont**: number of bins with more than 90% completeness and less than 10% contamination
* **\>0.5compl<0.05cont**: number of bins with more than 50% completeness and less than 5% contamination
* **\>0.7compl<0.05cont**: number of bins with more than 70% completeness and less than 5% contamination
* **\>0.9compl<0.05cont**: number of bins with more than 90% completeness and less than 5% contamination