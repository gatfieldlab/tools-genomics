#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Definition of protein structure related sequence alphabets
"""


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


SECONDARY_STRUCTURE = {
    'Str': 'Strand',  # beta-strand
    'Hlx': 'Helix',  # alpha-helix
    'Trn': 'Turn',  # turn
    'Uns': 'Unstructured',  # unstructured
    'Unk': 'Unknown'}  # unknown

TOPOLOGY = {
     'Trm': 'Transmembrane',
     'Itm': 'Intramembrane',
     'Cyt': 'Cytoplasmic',
     'Ext': 'Extracellular',
     'Lum': 'Lumenal',
     'Mti': 'Mito-intermembrane',
     'Mtm': 'Mito-matrix',
     'Unk': 'Unknown'}
