#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Utilities to work with genomics based python modules.
"""


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017-2020, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.1"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


P_INDEX = (9, 10, 8)
T_INDEX = (2, )
G_INDEX = (0, )
GT_INDEX = (0, 2)
F_INDEX = 11
S_INDEX = 3
I5_INDEX = 6
I3_INDEX = 7

class GenomicsUtilityException(Exception):
    """
    A simple class for all exceptions from genomics/utils.py
    """
    def __init__(self, m):
        message = "[genomics/utils] {}".format(m)
        super(GenomicsUtilityException, self).__init__(message)


def process_cds(cds_file, use_flag=False, t_index=None,
                only_single=False, include_incomplete=False):
    """
    Parses a specially formatted file containing CDS information. The file
    has to of 'prepared_cds' format
    Args:
         use_flag: :obj:`bool` to specify if the flag information from the
                   CDS file should be used (True) or not (False). If it is
                   used, any transcript models which are not flagged as 'passed'
                   by a "*" will not be used
         t_index: :obj:`tuple` of integers to identify the columns to be
                  concatanated to produce the unique ID. Use T_, G_ or GT_INDEX
         only_single: :obj:`bool` to specify if only transcript models that
                      are labeled as 'SINGLE_PROTEIN' - sole transcript isoform
                      from a gene that can code for a protein
         include_incomplete: :obj:`bool` to specify if transcript models that
                             have either incomplete start or end (but not both)
                             should be included. In this case, the incomplete
                             terminus will be clipped off so that CDS length is
                             multiple of 3
    Returns:
         tr_db: :obj:`dict` with key referring to transcripts, and value as a
         tuple of (CDS_start, CDS_end, TR_len).

    """
    tr_db = {}
    if t_index is None:
        t_index = T_INDEX
    try:
        with open(cds_file) as cds_handle:
            for line in cds_handle:
                parsed = line.strip().split('\t')
                tr_id = '|'.join([parsed[i] for i in t_index])
                is_incomplete = parsed[I5_INDEX] != parsed[I3_INDEX]
                is_unusable = parsed[I5_INDEX] == parsed[I3_INDEX] == "False"
                if (parsed[S_INDEX] == 'composite' or is_unusable or
                    (use_flag and parsed[F_INDEX] != '*') or
                    (only_single and parsed[1][:6] != 'SINGLE') or
                    (not include_incomplete and is_incomplete)):
                        continue
                c_start, c_end, c_len = tuple([int(parsed[i]) for i in P_INDEX])
                if is_incomplete:
                    extra_bases = (c_end - c_start) % 3
                    if parsed[I5_INDEX] == "False":
                        c_start = c_start + extra_bases
                    else:
                        c_end = c_end - extra_bases
                tr_db[tr_id] = (c_start, c_end, c_len)
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
