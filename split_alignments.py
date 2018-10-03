#!/home2/mattdoug/python3/bin/python3
# Last updated: 12/2/2018
# Author: Matt Douglas

from __future__ import print_function
import os, pysam, sys

def eprint(*args, **kwargs):
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def optype(path, op='r'):
    """If the file is BAM formatted, read/write as binary."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.bam':
        op += 'b'
    return op


def parse_CIGAR(cigar):
    """Parse a pysam formatted CIGAR string return positions of each skip."""
    ref_pos = 0
    str_pos = 0
    ref_pos_list = [0]
    str_pos_list = [0]
    cigar_list = []

    cigar_type = [i[0] for i in cigar]
    cigar_len  = [i[1] for i in cigar]

    prev = []
    for n, i in enumerate(cigar_type):
        j = cigar_len[n]
        if i == 0: # match
            ref_pos += j
            str_pos += j
            prev.append((i, j))
        elif i == 1: # insertion
            str_pos += j
            prev.append((i, j))
        elif i == 2: # deletion
            ref_pos += j
            prev.append((i, j))
        elif i == 3: # skip
            ref_pos += j
            ref_pos_list.append(ref_pos)
            str_pos_list.append(str_pos)
            cigar_list.append(prev)
            prev = []
        elif i == 4: # soft clipping
            str_pos += j
            prev.append((i, j))
    ref_pos_list.append(ref_pos)
    str_pos_list.append(str_pos)
    cigar_list.append(prev)

    return ref_pos_list, str_pos_list, cigar_list


def split_alignments(infile):
    for line in infile.fetch():
        pos = line.pos # SAM coordinates are 1-based
        seq = line.query_sequence
        qual = line.query_qualities
        cigar = line.cigar
        ref_pos_list, str_pos_list, cigar_list = parse_CIGAR(cigar)
        for i in range(len(str_pos_list) - 1):
            new_pos = pos + ref_pos_list[i]
            new_seq = seq[str_pos_list[i] : str_pos_list[i + 1]]
            if qual is not None:
                new_qual = qual[str_pos_list[i] : str_pos_list[i + 1]]
            else:
                new_qual = None
            # create the new line(s)
            new_line = line
            new_line.pos = new_pos
            new_line.query_sequence = new_seq
            new_line.query_qualities = new_qual
            new_line.cigar = cigar_list[i]
            new_line.template_length = 0
            yield new_line


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('USAGE: split_alignments.py input.bam output.bam')
        sys.exit(1)

    infile = pysam.AlignmentFile(sys.argv[1], optype(sys.argv[1], 'r'), check_sq=False)
    outfile = pysam.AlignmentFile(sys.argv[2], optype(sys.argv[2], 'w'), template=infile)

    for line in split_alignments(infile):
        outfile.write(line)
