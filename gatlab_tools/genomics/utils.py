#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Utilities to work with genomics based python modules.
"""


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


P_INDEX = (9, 10, 8)
T_INDEX = (2, )
G_INDEX = (0, )
GT_INDEX = (0, 2)
F_INDEX = 11
S_INDEX = 3


class GenomicsUtilityException(Exception):
    """
    A simple class for all exceptions from genomics/utils.py
    """
    def __init__(self, m):
        message = "[genomics/utils] {}".format(m)
        super(GenomicsUtilityException, self).__init__(message)


def process_cds(cds_file, use_flag=False, t_index=None):
    """
    Parses a specially formatted file containing CDS information. The file
    has to of 'prepared_cds' format
    """
    tr_db = {}
    if t_index is None:
        t_index = T_INDEX
    try:
        with open(cds_file) as cds_handle:
            for line in cds_handle:
                parsed = line.strip().split('\t')
                tr_id = '|'.join([parsed[i] for i in t_index])
                if (parsed[S_INDEX] == 'composite' or
                    (use_flag and parsed[F_INDEX] != '*')):
                        continue
                tr_db[tr_id] = tuple([int(parsed[i]) for i in P_INDEX])
    except IOError:
        raise GenomicsUtilityException(
            "Could not read the CDS info: {}\n".format(cds_file))
    except IndexError:
        raise GenomicsUtilityException(
            "Could not parse the CDS info: {}\n".format(cds_file))
    return tr_db


def process_gene_cds(cds_file, use_flag=False):
    """
    Parses a specially formatted file containing CDS information. The
    difference from process_cds() is that it returns a dict of genes
    rather than transcripts
    """
    gene_tr_db = {}
    try:
        with open(cds_file) as cds_handle:
            for line in cds_handle:
                parsed = line.strip().split('\t')
                gene_id = '|'.join([parsed[i] for i in G_INDEX])
                tr_id = '|'.join([parsed[i] for i in T_INDEX])
                if (parsed[S_INDEX] == 'composite' or
                    (use_flag and parsed[F_INDEX] != '*')):
                    continue
                if gene_id not in gene_tr_db:
                    gene_tr_db[gene_id] = {}
                gene_tr_db[gene_id][tr_id] = tuple(
                    [int(parsed[i]) for i in P_INDEX])
    except IOError:
        raise GenomicsUtilityException(
            "Could not read the CDS info: {}\n".format(cds_file))
    except IndexError:
        raise GenomicsUtilityException(
            "Could not parse the CDS info: {}\n".format(cds_file))
    return gene_tr_db


def get_trs(tr_file):
    """
    Parses (any white space) and extracts TR-ids from a flat file
    """
    if not tr_file:
        return None
    tr_list = []
    try:
        with open(tr_file) as tr_f:
            for line in tr_f:
                tr_list.extend(line.strip().split())
    except IOError:
        raise GenomicsUtilityException(
            "Could not read the TR-IDs file: {}\n".format(tr_file))
    return tr_list
