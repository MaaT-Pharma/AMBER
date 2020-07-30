#!/usr/bin/env python

import collections

TOOL = "tool"
PRECISION = "avg_purity"
RECALL = "avg_completeness"
AVG_PRECISION = "avg_purity_per_bin"
STD_DEV_PRECISION = "std_dev_purity_per_bin"
SEM_PRECISION = "sem_purity_per_bin"
AVG_RECALL = "avg_completeness_per_bin"
STD_DEV_RECALL = "std_dev_completeness_per_bin"
SEM_RECALL = "sem_completeness_per_bin"
RI_BY_BP = "rand_index_by_bp"
RI_BY_SEQ = "rand_index_by_seq"
ARI_BY_BP = "a_rand_index_by_bp"
ARI_BY_SEQ = "a_rand_index_by_seq"
PERCENTAGE_ASSIGNED = "percent_assigned"
PERCENTAGE_ASSIGNED_UNIQUE = "percent_unique_assigned"
ACCURACY = "accuracy"

abbreviations = collections.OrderedDict([(AVG_PRECISION, 'purity averaged over genome bins'),
                 (STD_DEV_PRECISION, 'standard deviation of purity averaged over genome bins'),
                 (SEM_PRECISION, 'standard error of the mean of purity averaged over genome bins'),
                 (AVG_RECALL, 'completeness averaged over genome bins'),
                 (STD_DEV_RECALL, 'standard deviation of completeness averaged over genome bins'),
                 (SEM_RECALL, 'standard error of the mean of completeness averaged over genome bins'),
                 (PRECISION, 'average purity per base pair or sequence counts'),
                 (RECALL, 'average completeness per base pair or sequence counts'),
                 ('rand_index_by_bp', 'Rand index weighed by base pairs'),
                 ('rand_index_by_seq', 'Rand index weighed by sequence counts'),
                 ('a_rand_index_by_bp', 'adjusted Rand index weighed by base pairs'),
                 ('a_rand_index_by_seq', 'adjusted Rand index weighed by sequence counts'),
                 ('percent_assigned', 'percentage of  sequences or base pairs that were assigned to bins'),
                 ('percent_unique_assigned', 'percentage of unique sequences or base pairs that were assigned to bins'),
                 ('accuracy', 'accuracy'),
                 ('>0.5compl<0.1cont', 'number of bins with more than 50% completeness and less than 10% contamination'),
                 ('>0.7compl<0.1cont', 'number of bins with more than 70% completeness and less than 10% contamination'),
                 ('>0.9compl<0.1cont', 'number of bins with more than 90% completeness and less than 10% contamination'),
                 ('>0.5compl<0.05cont', 'number of bins with more than 50% completeness and less than 5% contamination'),
                 ('>0.7compl<0.05cont', 'number of bins with more than 70% completeness and less than 5% contamination'),
                 ('>0.9compl<0.05cont', 'number of bins with more than 90% completeness and less than 5% contamination')])
