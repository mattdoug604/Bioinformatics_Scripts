#!/usr/bin/python

import sys

def main():
    switch_out = lambda x: pair2 if x == pair1 else pair1    
    
    infile = open(sys.argv[1], 'r')
    pair1 = open(sys.argv[1].split('/')[-1].split('.')[0]+'_1.fq', 'w')
    pair2 = open(sys.argv[1].split('/')[-1].split('.')[0]+'_2.fq', 'w')    

    outfile = pair1
    count = 0
    while True:
        line = infile.readline()
        if not line:
            break
        else:
            if count >= 4: 
                outfile = switch_out(outfile)
                count = 1
            else:
                count += 1
            outfile.write(line)
    
    infile.close()
    pair1.close()
    pair2.close()

if __name__ == '__main__':
    main()

