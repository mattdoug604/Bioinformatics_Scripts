#!/usr/local/bin/python

from __future__ import print_function
import sys
import pysam
import random

random.seed('mattdoug')   # use the same seed for reporducibility

# 737,824,345 alignments (wc -l)
# 737,824,345 (this progam)

def eprint(*args, **kwargs):
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def count_alignments(bamfile):
    """Return the number of alignments in the BAM file."""
    for line_num, _ in enumerate(bamfile.fetch()):
        pass

    return line_num + 1


def build_sample_set(sample_size, num_alignments):
    """Return a list of sets, each with X randomly generated line numbers where
    X is the given sample size.
    """
    sample_set = set()
    sample_size = sample_size * 1000000

    while len(sample_set) < sample_size:
        num = random.randrange(num_alignments)
        sample_set.add(num)

    return sample_set


def find_supporting_alignments(bamfile, outfile, sample_set, num_alignments):
    total = len(sample_set)
    count = 0

    for n, line in enumerate(bamfile.fetch()):
        if n in sample_set:
            count += 1
            sample_set.remove(n)
            outfile.write(line)
            eprint('\r{}% ({}/{})'.format(count*100//total, count, total), end= '')
        if len(sample_set) <= 0:
            eprint(' ')
            break


if __name__ == '__main__':
    # parse command line arguments
    bamfile = pysam.AlignmentFile(sys.argv[1], 'rb')
    eprint('Alignment file: {}'.format(sys.argv[1]))

    sample_size = int(sys.argv[2])
    eprint('Sample size: {:,}M'.format(sample_size))

    if len(sys.argv) >= 4:
        num_alignments = int(sys.argv[3])
    else:
        eprint('Counting alignments...')
        num_alignments = count_alignments(bamfile)
    eprint('Total number of alignments: {:,}'.format(num_alignments))

    # randomly generate sets of line numbers to sample
    fname = 'out{}M.bam'.format(sample_size)
    outfile = pysam.AlignmentFile(fname, 'wb', template=bamfile)

    eprint('Building random sample sets...')
    sample_set = build_sample_set(sample_size, num_alignments)

    eprint('Sampling BAM file...')
    find_supporting_alignments(bamfile, outfile, sample_set, num_alignments)
