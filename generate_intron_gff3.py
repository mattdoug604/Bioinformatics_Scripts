#!/home2/mattdoug/python3/bin/python3
# Last updated: 1/15/2018
# Author: Matt Douglas
# Purpose: Read through a SAM/BAM file and generate a GFF3 file of all the
#          supported introns. Alternatively, supply a GFF3 file of introns and
#          count read support for each one.

import pysam
import argparse, os, sys
from collections import defaultdict

def parse_commandline_arguments():
    parser = argparse.ArgumentParser(description='Count the number of supporting reads for each intron.')
    parser.add_argument('-a',
                        type=str,
                        nargs='+',
                        help='one or more SAM/BAM file of alignments')
    parser.add_argument('-i',
                        type=str,
                        nargs='?',
                        help='(optional) limit search to a list of introns in GFF3 format')
    parser.add_argument('-m',
                        type=int,
                        nargs='?',
                        default=0,
                        help='minimum required read support (default=0)')
    parser.add_argument('-s',
                        '--strand-only',
                        action='store_true',
                        help='discard any introns without a defined strand')
    args = parser.parse_args()

    if args.a is None:
        parser.print_help()
        sys.exit(1)

    return args.a, args.i, args.m, args.strand_only


def eprint(*args, **kwargs):
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def optype(path, op='r'):
    """If the file is BAM formatted, read/write as binary."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.bam':
        op += 'b'
    return op


def parse_CIGAR(chrom, pos, cigar):
    """Get the pos, end, and size of any introns from the CIGAR string.
    if(  cigar_type == 0): #match
    elif(cigar_type == 1): #insertions
    elif(cigar_type == 2): #deletion
    elif(cigar_type == 3): #skip
    elif(cigar_type == 4): #soft clipping
    elif(cigar_type == 5): #hard clipping
    elif(cigar_type == 6): #padding
    """
    introns = list()

    cigar_type = [i[0] for i in cigar]
    cigar_len = [i[1] for i in cigar]

    for i in [i for i, l in enumerate(cigar_type) if l == 3]:
        size = cigar_len[i]
        start = pos
        for j in range(len(cigar_type[:i])):
            if cigar_type[j] in [0, 2, 3]:
                start += cigar_len[j]
        end = start + size - 1

        introns.append([chrom, start, end])

    return introns


def find_introns(bamfile, count_dict):
    """Read though the BAM file and find all introns, as specified in the CIGAR
    string, with number of supporting reads.
    """
    for line in bamfile.fetch():
        chrom = bamfile.get_reference_name(line.rname)
        pos = line.pos + 1
        cigar = line.cigar
        if line.has_tag('XS'):
            strand = line.get_tag('XS')
        else:
            strand = '.'
        for intron in parse_CIGAR(chrom, pos, cigar):
            intron = tuple(intron + [strand])
            count_dict[intron] += 1

    return count_dict


def parse_gff3(path):
    intron_set = set()

    with open(path, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.split('\t')
                intron = line[0], int(line[3]), int(line[4]), line[6]
                intron_set.add(intron)

    return intron_set


def sort_features(introns):
    """Sort tuples of introns by chromosome, then start position, then end
    position. NOTE: C. elegans uses roman numerals for chromosomes names.
    """
    numerals = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'X':10, 'MtDNA':11}
    try:
        return sorted(introns, key=lambda x: (numerals[x[0]], int(x[1]), int(x[2])))
    except KeyError:
        return sorted(introns, key=lambda x: (x[0], int(x[1]), int(x[2])))


def output_as_gff3(count_dict):
    """Ouptut the results in GFF3 format."""
    print('##gff-version 3')
    for n, intron in enumerate(sort_features(count_dict)):
        chrom, start, end, strand = intron
        count = count_dict[intron]
        line = chrom, '.', 'intron', start, end, count, strand, '.', 'ID='+str(n+1)
        print(*line, sep='\t')


def run(bamfiles, gff_path):
    count_dict = defaultdict(int)

    for b in bamfiles:
        count_dict = find_introns(b, count_dict)
    eprint('Found {:,} introns in {} file(s)'.format(len(count_dict), len(bamfiles)))

    if gff_path is not None:
        intron_set = parse_gff3(gff_path)
        eprint('{:,} introns specified in file: {}'.format(len(intron_set), gff_path))
        # add in any introns that were in the GFF3 file, but not found
        for intron in intron_set:
            if intron not in count_dict:
                count_dict[intron] = 0
        # only report introns in the specified set
        for intron, _ in count_dict.items():
            if intron not in intron_set:
                del count_dict[intron]

    # discard any introns not meeting filtering criteria
    del_set = set()
    m, n = 0, 0
    for intron, count in count_dict.items():
        to_del = False
        if count < min_count:
            to_del = True
            m += 1
        if intron[3] not in ('+', '-') and strand_only:
            to_del
            n += 1
        if to_del:
            del_set.add(intron)
    for intron in del_set:
        del count_dict[intron]

    if m > 0:
        eprint('  Discarding {:,} introns with less than {} support'.format(m, min_count))
    if n > 0:
        eprint('  Discarding {:,} introns without a defined strand'.format(n))

    return count_dict


if __name__ == '__main__':
    bam_paths, gff_path, min_count, strand_only = parse_commandline_arguments()
    bamfiles = [pysam.AlignmentFile(b, optype(b, op='r')) for b in bam_paths]
    count_dict = run(bamfiles, gff_path)
    output_as_gff3(count_dict)
