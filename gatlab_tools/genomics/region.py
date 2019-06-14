#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Region module to define and extract easily different regions from
transcripts' features as such CDS, UTRs etc
"""


import sys
import re
import operator
import argparse
from genomics import utils


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.2.2"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


REGIONS = ('5utr', 'cds', '3utr', 'transcript', 'custom')
CDS_START, CDS_END, TR_LEN = range(3)
BOX_START, BOX_END = range(2)
REGION_PATTERN = (r'^((?P<region>{region})|'
                  r'(?P<number>\d+)|'
                  r'(?P<star>(?P<pre>\*)?(?=({region}))'
                  r'({region})(?(pre)|\*?))(?P<op>[\+\-])'
                  r'(?P<offset>\d+%?))$'.format(region='|'.join(REGIONS)))
REGION_MATCH = re.compile(REGION_PATTERN)
OPS = {'+': operator.add, '-': operator.sub}


class RegionException(Exception):
    """
    A simple class for all exceptions from region.py
    """

    def __init__(self, m):
        message = "[calc_codon] {}".format(m)
        super(RegionException, self).__init__(message)


class Region(object):
    """
    Main class for defining a region. Examples:
    cds = [cds_start, cds_end)
    cds:3utr = [cds_start, 3utr_end)
    *cds = [cds_start, cds_start)
    *cds:cds* = [cds_start, cds_end)
    *cds-50:cds = [cds_start-50, cds_end)
    *cds-50:cds*+50 = [cds_start-50, cds_end+50)
    *cds+40:3utr*-20 = [cds_start+40, 3utr_end-20]
    transcript = [0, transcript_end)
    *transcript+10:transcript*-10 = [10, transcript_end-10)
    10:3utr*+10 = [10, 3utr_end-10) = [10, transcript_end-10)
    """

    def __init__(self, region):
        self.parse_words = []
        parsed_regions = region.split(":")
        for parsed_region in parsed_regions:
            match = REGION_MATCH.search(parsed_region)
            if match:
                if match.group('region'):
                    self.parse_words.append([match.group('region')])
                elif match.group('number'):
                    self.parse_words.append([match.group('number')])
                else:
                    self.parse_words.append([match.group('star'),
                                             match.group('op'),
                                             match.group('offset')])
            else:
                raise RegionException(
                    "'{}' could not be match to region regex".format(
                        parsed_region))

    def get_limits(self, tr_info, enforce=False):
        """
        Given a tr_info tuple (CDS_START, CDS_END, TR_LEN), or similarly
        (BOX_START, BOX_END, TR_LEN), calculates the
        limits of for the Region instance. The limits are 0-based , half closed
        intervals [start, end).
        """
        limits = []
        for sentence in self.parse_words:
            limits.append([])
            cur_word_limits = None
            for word in sentence:
                if word.isdigit():
                    limits[-1].append(int(word))
                elif word[-1] == '%' and word[:-1].isdigit():
                    if not cur_word_limits:
                        raise RegionException(
                            'Fatal error during % calc: {}'.format(sentence))
                    limits[-1].append((cur_word_limits[1] -
                                       cur_word_limits[0]) // 100 * int(word[:-1]))
                elif word in ['+', '-']:
                    limits[-1].append(word)
                else:
                    left, right = word[0] == '*', word[-1] == '*'
                    if left:
                        word = word[1:]
                    if right:
                        word = word[:-1]
                    if word == '5utr':
                        word_limits = (0, tr_info[CDS_START])
                    elif word == 'cds':
                        word_limits = (tr_info[CDS_START], tr_info[CDS_END])
                    elif word == 'custom':
                        word_limits = (tr_info[BOX_START], tr_info[BOX_END])
                    elif word == '3utr':
                        word_limits = (tr_info[CDS_END], tr_info[TR_LEN])
                    elif word == 'transcript':
                        word_limits = (0, tr_info[TR_LEN])
                    cur_word_limits = word_limits
                    if left:
                        limits[-1].append(word_limits[0])
                    if right:
                        limits[-1].append(word_limits[1])
                    if not (left or right):
                        limits[-1].extend(word_limits)
        result = []
        for limit in limits:
            if len(limit) == 1 and not isinstance(limit[0], int):
                raise RegionException('Fatal error: {}'.format(limits))
            while '-' in limit or '+' in limit:
                o_i, optr = next((i, o) for i, o in enumerate(limit) if
                                 o == '-' or o == '+')
                try:
                    res = OPS[optr](limit[o_i-1], limit[o_i+1])
                except KeyError:
                    raise RegionException('Fatal error: {}'.format(limits))
                limit = limit[:max(0, o_i-1)] + [res] + limit[min(len(limit),
                                                                  o_i+2):]
            result.append(limit)
        lmin = min([min(l) for l in result])
        lmax = max([max(l) for l in result])
        if enforce:
            lmin = max(0, lmin)
            lmin = min(tr_info[TR_LEN], lmin)
            lmax = min(tr_info[TR_LEN], lmax)
            lmax = max(0, lmax)
        return (lmin, lmax)

    def __str__(self):
        return "Region: {}".format(self.parse_words)


def main():
    """
    Runs the CLI for calculating the limits of a given region for the
    supplied transcripts from adatabase
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('region', help='Region of the transcript', type=str)
    parser.add_argument('transcript', help='Transcript ID', nargs='?')
    parser.add_argument('-t', '--tr-ids', type=str, help='A file with TR-IDs')
    parser.add_argument('-c', '--cds', type=str, default='prepared_cds.txt',
                        required=True)
    parser.add_argument('-f', '--flag', action='store_true', help='Does the '
                        'ĈDS database file have a flag column and if so, '
                        'should it be used?')
    parser.add_argument('-n', '--no-enforce', action='store_false', help='Do '
                        'not enforce the returned limits to be within 0 - '
                        'transcript length')
    args = parser.parse_args()
    try:
        tr_ids = utils.get_trs(args.tr_ids)
    except utils.GenomicsUtilityException as utils_exception:
        sys.stderr.write("{}".format(utils_exception))
        return 1
    if args.transcript:
        tr_ids.append(args.transcript)
    if not tr_ids:
        sys.stderr.write('No TR-IDs were given to extract the region from\n')
        return 1
    try:
        tr_db = utils.process_cds(args.cds, use_flag=args.flag)
    except utils.GenomicsUtilityException as utils_exception:
        sys.stderr.write("{}".format(utils_exception))
        return 1
    try:
        region = Region(args.region)
    except RegionException as region_exception:
        sys.stderr.write("{}\n".format(region_exception))
        return 1
    for tr_id in tr_ids:
        try:
            limits = region.get_limits(tr_db[tr_id], enforce=args.no_enforce)
        except KeyError:
            sys.stderr.write("Could not find '{}' in DB\n".format(tr_id))
        except RegionException as region_exception:
            sys.stderr.write("Could not extract limits for {}: {}\n".format(
                tr_id, region_exception))
        else:
            limits = "[{}, {})".format(*limits)
            sys.stdout.write("{}: {}\n".format(tr_id, limits))
    return 0
