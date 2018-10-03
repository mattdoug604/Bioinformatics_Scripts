#!/usr/local/bin/python3
#last updated: 2/3/2018

import sys, random

def sort_lines(lines):
    """Sort GFF3 formatted entries by chromosome, then start position, then end
    position, then strand. NOTE: C. elegans uses roman numerals for chromosomes names.
    """
    numerals = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'X':10, 'MtDNA':11}
    lines = [i.split('\t') for i in lines]

    try:
        sorted_lines = sorted(lines, key = lambda x: (numerals[x[0]], int(x[3]), int(x[4]), x[6]))
    except KeyError:
        sorted_lines = sorted(lines, key = lambda x: (x[0], int(x[3]), int(x[4]), x[6]))
    sorted_lines = ['\t'.join(i) for i in sorted_lines]

    return sorted_lines


if __name__ == '__main__':
    infile = sys.argv[1]
    number = [int(sys.argv[2]) if len(sys.argv) > 2 else 1][0]
    output = []

    with open(infile, 'r') as f:
        lines = [i.strip() for i in f.readlines() if i[0] != '#']
        for line in random.sample(lines, number):
            output.append(line)

    for line in sort_lines(output):
        print(line)
