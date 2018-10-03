#!/usr/local/bin/python3

# STAR tab format:
#   0) chromosome
#   1) start (1-based)
#   2) end (1-based)
#   3) strand
#   4) motif
#   5) annotated?
#   6) no. uniquely mapped reads
#   7) no. multi-mapped reads
#   8) max splice overhang
# GFF3 format: seqid, source, type, start, end, score, strand, phase, attributes

import sys

strand_dict = {'0':'.', '1':'+', '2':'-'}

if __name__ == '__main__':
    with open(sys.argv[1], 'r') as f:
        print('##gff-version 3')
        for n, line in enumerate(f, start=1):
            i = line.strip().split('\t')
            try:
                new_line = (i[0],
                            '.',
                            'intron',
                            i[1],
                            i[2],
                            int(i[6]) + int(i[7]),
                            strand_dict[i[3]],
                            '.',
                            'ID='+str(n))
            except IndexError:
                print('Incomplete line or unexpected format at line {}: {}').format(n, line)
                sys.exit(1)
            print(*new_line, sep='\t')
