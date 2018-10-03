#!/usr/local/bin/python3

#BED format: chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
#GFF3 format: seqid, source, type, start, end, score, strand, phase, attributes

import sys

if __name__ == '__main__':
    n = 0

    with open(sys.argv[1], 'r') as f:
        print('##gff-version 3')
        for line_no, line in enumerate(f):
            i = line.strip().split('\t')
            if len(i) > 1:
                n += 1
                try:
                    new_line = (i[0],
                                '.',
                                'intron',
                                int(i[1]) + int(i[10].split(',')[0]) + 1,
                                int(i[2]) - int(i[10].split(',')[1]),
                                i[4],
                                i[5],
                                '.',
                                'ID=' + str(n))
                except IndexError:
                    print('Incomplete line or unexpected format at line {}: {}').format(line_no, line)
                    sys.exit(1)
                print(*new_line, sep='\t')
