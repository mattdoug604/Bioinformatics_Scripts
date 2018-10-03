#!/usr/local/bin/python3

# INPUT:
# II:2107356-2107465 (- strand)
# OUTPUT:
# II   .   intron   2107356   2107465   .   -   .   .

import sys

with open(sys.argv[1], 'r') as f:
    print('##gff-version 3')
    for line in f:
        # split the line
        line        = line.strip().replace(',', '')
        chrom, line = line.split(':')
        beg, line   = line.split('-', 1)
        end, line   = line.split(' ', 1)
        strand      = line[1]
        # print in GFF3 format
        print(chrom, '.', 'exon', beg, end, 0, strand, '.', '.', sep='\t')
