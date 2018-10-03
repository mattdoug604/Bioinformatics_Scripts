#!/usr/local/bin/python
# Last updated: 21/12/2017
# Author: Matt Douglas

from __future__ import print_function, division
import argparse, re, os, sys, pysam
from collections import defaultdict
from itertools import chain
from time import time

class Progress(object):
    """A timer to monitor the number of lines read, and number of supporting
    alignments found. Reports these statistics ever X seconds, where X is some
    interval in seconds.

    Attributes:
        time_start: The time the object was created.
        time_prev: The last time the object reported.
        interval: Minimum period of time needed to pass before reporting again.
        num_lines: The number of lines in the BAM file (to track progess)
    """

    def __init__(self, num_lines):
        """Return a timer object."""
        self.time_prev = time()
        self.time_start = self.time_prev
        self.num_lines = num_lines

    def update(self, line_count, found_count):
        """If more than 'interval' seconds has passed. Print progress."""
        if time() - self.time_prev > 1:
            self.time_prev = time()
            self.p_complete = (line_count * 100) / self.num_lines
            eprint('\r {:,} lines read ({}% of total)... {:,} supporting alignments found...'
                   .format(line_count, '%.1f' % self.p_complete, found_count), end='')

    def elapsed(self):
        """Return the elapsed time in minutes."""
        return (time() - self.time_start) / 60


#####################
# Utility functions #
#####################
def eprint(*args, **kwargs):
    """Print to stderr."""
    if not quiet:
        print(*args, file=sys.stderr, **kwargs)


def optype(path, op='r'):
    """If the file is BAM formatted, read/write as binary."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.bam':
        op += 'b'
    return op


def format_intron(intron):
    """Print a intron tuple of (chromosome, start position, end position) into
    a readable format.
    """
    return intron[0] + ':' + str(intron[1]) + '-' + str(intron[2])


def count_lines(samfile):
    """Count the number of lines in a BAM file."""
    count = 0

    try:
        for line in pysam.idxstats(samfile).split('\n'):
            if len(line) > 0:
                count += eval( '+'.join( line.split('\t')[2:] ))
    except pysam.utils.SamtoolsError as e:
        eprint('[ERROR]', e)
        eprint('Is the BAM file indexed?')
        sys.exit(1)

    return count


################################
# Parse command line arguments #
################################
def parse_intron(intron_coord):
    """Parse a set of intron coordinates into a tuple of (chromosome, start
    position, end position).
    """
    intron_coord = intron_coord.strip().replace(',', '').replace('..', '-')
    chrom, start, end = re.split(':|-|_', intron_coord)
    return chrom, int(start), int(end)


def parse_GFF3(GFF3_file):
    """Parse a GFF3 file and return all intron coordinates."""
    for line in GFF3_file:
        if line[0] != '#':
            col = line.split('\t')
            chrom, start, end = col[0], int(col[3]), int(col[4])
            yield chrom, start, end


def parse_TSV(TSV_file):
    """Parse a tab-seperated file from HISAT2 and return all intron
    coordinates."""
    for line in TSV_file:
        col = line.split('\t')
        chrom, start, end = col[0], int(col[1])+2, int(col[2])
        yield chrom, start, end


def parse_commandline_arguments():
    """Parse command line arguments and return:
        1) A SAM formatted file
        2) A list of introns to search for
        3) A Bool specifying whether to search for alternative alignments or not
    """
    parsed_introns = []

    # parse command line options
    parser = argparse.ArgumentParser(description='Return alignments supporting one or more specified introns.')
    parser.add_argument('-a', type=str, nargs='?', help='a SAM or BAM file')
    parser.add_argument('-i', type=str, nargs='+', help="one or more introns on the command line in the format 'I:1234..1345'")
    parser.add_argument('-g', type=str, nargs='?', help='a GFF3 file of introns to search for')
    parser.add_argument('-t', type=str, nargs='?', help='a tab-seperated file of introns to search for')
    parser.add_argument('-r', action='store_true', help='report all alternative alignments, and paired alignments, for supporting reads')
    parser.add_argument('-o', type=str, nargs='?', help='output file (if not specified, each intron will have a seperate file)')
    parser.add_argument('-q', action='store_true', help='quiet mode (do not print progress)')
    args = parser.parse_args()

    # parse each intron specifed in the command line (if any)
    if args.i is not None:
        for i in args.i:
            parsed_introns.append(parse_intron(i))

    # parse each intron listed in a GFF3 file (if any)
    if args.g is not None:
        with open(args.g, 'r') as f:
            for i in parse_GFF3(f):
                parsed_introns.append(i)

    # parse each intron listed in a TSV file (if any)
    if args.t is not None:
        with open(args.t, 'r') as f:
            for i in parse_TSV(f):
                parsed_introns.append(i)

    # an alignment file and one source of introns must be specified
    if args.a is None or len(parsed_introns) < 1:
        parser.print_help()
        sys.exit(1)

    return parsed_introns, args.a, args.r, args.o, args.q


#######################################
# Functions for doing the actual work #
#######################################
def regions_to_search(parsed_introns, read_length=150):
    """Return a set of non-overlapping regions to search on each chromosome."""
    flatten = chain.from_iterable
    data = defaultdict(list)

    # process by chromosome
    for i in parsed_introns:
        data[i[0]].append([i[1], i[2]])

    # find overlapping regions on each chromosome
    for chrom, ranges in data.items():
        ranges = sorted(flatten(((start, 1), (end+read_length, -1)) for start, end in ranges))

        c, x = 0, 0
        for value, label in ranges:
            if c == 0:
                x = value
            c += label
            if c == 0:
                yield chrom, x-read_length, value


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

        introns.append((chrom, start, end))

    return introns


def find_supporting_alignments(samfile, parsed_introns):
    """Read though the BAM file and find all introns, as specified in the CIGAR
    string, with number of supporting reads.
    """
    alignments_matching = defaultdict(list)
    line_count = 0
    found_count = 0

    for chrom, start, end in regions_to_search(parsed_introns):
        for line in samfile.fetch(chrom, start, end):
            line_count += 1
            progress.update(line_count, found_count)
            # parse alignment
            chrom = samfile.get_reference_name(line.reference_id)
            pos = line.pos + 1  # SAM coordinates are 1-based
            cigar = line.cigar
            # find introns
            for intron in parse_CIGAR(chrom, pos, cigar):
                if intron in parsed_introns:
                    found_count += 1
                    alignments_matching[intron].append(line)

    eprint('\r {:,} lines read. {:,} supporting alignments found!{}'.format(line_count, found_count, ' '*20))

    return alignments_matching


def report_all_alignments(samfile, alignments_matching):
    """Report all alignments for reads, and paired-reads, supporting the
    specified introns
    """
    reads_supporting = defaultdict(set)
    alignments_matching_all = defaultdict(set)
    line_count = 0
    found_count = 0

    # get the reads that support the intron(s)
    for intron, lines in alignments_matching.items():
        for read in [i.qname for i in lines]:
            reads_supporting[read].add(intron)

    for line in samfile.fetch():
        line_count += 1
        progress.update(line_count, found_count)
        # parse alignment
        read = line.qname
        if read in reads_supporting:
            for intron in reads_supporting[read]:
                found_count += 1
                alignments_matching_all[intron].add(line)

    eprint('\r', ' '*79, end='')
    eprint('\r Reporting {:,} total alignments!'.format(found_count))

    return alignments_matching_all


######################
# Output the results #
######################
def print_to_individual_files(samfile, parsed_introns, alignments_matching, ext='sam'):
    for intron in parsed_introns:
        if intron in alignments_matching:
            output_path = '_'.join(map(str, intron)) + '.' + ext
            outfile = pysam.AlignmentFile(output_path, optype(output_path, 'w'), template=samfile)
            for line in alignments_matching[intron]:
                 outfile.write(line)
            outfile.close()
            eprint(' Wrote {:,} lines to: {}'.format(len(alignments_matching[intron]), output_path))
        else:
            eprint(' Could not find support for intron at {}'.format(format_intron(intron)))


def print_all_to_one_file(samfile, parsed_introns, alignments_matching, output_path):
    alignment_set = set()

    #  output file is either SAM or BAM, depending the output_path extension
    outfile = pysam.AlignmentFile(output_path, optype(output_path, 'w'), template=samfile)

    # alignments can support more than one intron, so remove duplicate entries
    for intron, lines in alignments_matching.items():
        for line in lines:
            alignment_set.add(line)

    # print the output
    for line in alignment_set:
        outfile.write(line)

    eprint(' Wrote {:,} lines to {}'.format(len(alignment_set), output_path))

    outfile.close()


#############
# Main loop #
#############
if __name__ == '__main__':
    # parse commandline arguments
    parsed_introns, input_path, report_all, output_path, quiet = parse_commandline_arguments()
    samfile = pysam.AlignmentFile(input_path, optype(input_path, 'r'))
    eprint('{:,} intron{} to search for.'.format(len(parsed_introns), ['s' if len(parsed_introns) != 1 else ''][0]))

    # get the number of lines, to keep track of progress
    num_lines = count_lines(input_path)
    progress = Progress(num_lines)

    # find supporting alignments
    eprint('Searching for supporting alignments:')
    alignments_matching = find_supporting_alignments(samfile, parsed_introns)
    if report_all:
        eprint('Searching for read mates and alternative alignments for supporting reads:')
        alignments_matching = report_all_alignments(samfile, alignments_matching)

    # output the results
    eprint('Writing output:')
    if output_path is not None:
        print_all_to_one_file(samfile, parsed_introns, alignments_matching, output_path)
    else:
        print_to_individual_files(samfile, parsed_introns, alignments_matching)

    samfile.close()
    eprint('Done! (runtime = {}min)'.format('%.2f' % progress.elapsed()))
