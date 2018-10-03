#!/usr/local/bin/python
# Last updated: 1/29/2018
# Author: Matt Douglas

from __future__ import print_function
import os
import pysam
import sys
from collections import defaultdict

def optype(path, op='r'):
    """If the file is BAM formatted, read/write as binary."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.bam':
        op += 'b'
    return op


def get_strand(line):
    if line.is_reverse:
        return '-'
    else:
        return '+'


def parse_CIGAR(cigar):
    """Parse a pysam formatted CIGAR string return positions of each skip."""
    ref_pos = 0
    str_pos = 0
    ref_pos_list = [0]
    length_list = []

    cigar_type = [i[0] for i in cigar]
    cigar_len  = [i[1] for i in cigar]

    length = 0
    for n, i in enumerate(cigar_type):
        j = cigar_len[n]
        if i == 0: # match
            ref_pos += j
            length += j
        elif i == 2: # deletion
            ref_pos += j
            length += j
        elif i == 3: # skip
            ref_pos += j
            ref_pos_list.append(ref_pos)
            length_list.append(length)
            length = 0
        elif i == 4: # soft clipping
            length += j
    ref_pos_list.append(ref_pos)
    length_list.append(length)

    return ref_pos_list, length_list


def convert_alignment_to_tuple(samfile):
    for line in samfile.fetch():
        chrom = samfile.get_reference_name(line.reference_id)
        pos = line.pos + 1 # SAM coordinates are 1-based
        seq = line.query_sequence
        qual = line.query_qualities
        cigar = line.cigartuples
        strand = get_strand(line)
        ref_pos_list, length_list = parse_CIGAR(cigar)
        for i in range(len(ref_pos_list)-1):
            start = pos + ref_pos_list[i]
            end = start + length_list[i] - 1
            yield chrom, start, end, strand


def sort_by_pos(exons):
    """Sort tuples of exons by chromosome, then start position, then end
    position. NOTE: C. elegans uses roman numerals for chromosomes names.
    """
    numerals = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'X':10, 'MtDNA':11}
    try:
        return sorted(exons, key=lambda x: (numerals[x[0]], int(x[1]), int(x[2])))
    except KeyError:
        return sorted(exons, key=lambda x: (x[0], int(x[1]), int(x[2])))


def print_as_gff3(count_dict, outfile):
    sorted_feats = sort_by_pos(count_dict.keys())

    with open(outfile, 'w') as f:
        print('#gff3-version 3', file=f)
        for feature in sorted_feats:
            chrom, start, end, strand = feature
            count = count_dict[feature]
            line = chrom, '.', 'exon', start, end, count, strand, '.', '.'
            print(*line, sep='\t', file=f)


if __name__ == '__main__':
    count_dict = defaultdict(int)

    infile = sys.argv[1]
    outfile = sys.argv[2]
    samfile = pysam.AlignmentFile(infile, optype(infile, 'r'))

    for feature in convert_alignment_to_tuple(samfile):
        count_dict[feature] += 1

    print_as_gff3(count_dict, outfile)
