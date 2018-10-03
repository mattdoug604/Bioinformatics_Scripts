#!/usr/local/bin/python3

#BED format: chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
#GFF3 format: seqid, source, type, start, end, score, strand, phase, attributes

import sys

# GFF3 header for WormBase C. elegans WS250
#print('track name=' + sys.argv[1].split('.')[0], 'introns description=\"none\"')

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] != '#':
            chrom, source, feat, start, end, score, strand, phase, attr = line
            new_line = (chrom,
                        start,
                        end,
                        feat + str(n),
                        score,
                        strand,
                        start,
                        end,
                        '255,0,0',
                        '1',
                        str(int(end) - int(start) + 1),
                        '0')
            print(new_line, sep='\t')
