#!/usr/local/bin/python3
#Last updated: 31/7/2017

import sys
from collections import defaultdict

def parse_attr(lines):
    """Parse attributes column in format 'tag=value' and return a dictionary in
    the form tag:values
    """
    attr_dict = defaultdict(set)

    for line in lines:
        attrs = line.split('\t')[8].strip().split(';')
        for i in attrs:
            tag, val = i.split('=', 1)
            attr_dict[tag].add(val)

    return attr_dict


def parse_GFF3(file_list):
    """Parse a GFF3 format file and merge entries of the same type and
    position.
    """
    header = ['##gff-version 3']
    entries = defaultdict(list)
    new_entries = []

    # group entries by feature position
    for n, infile in enumerate(file_list):
        header.append('##File {} = {}'.format(n+1, infile))
        with open(infile, 'r') as f:
            for line in f:
                line = line.strip()
                if line[0] != '#':
                    i = line.split('\t')
                    entry = i[0], i[2], int(i[3]), int(i[4]), i[6]
                    entries[entry].append(line)

    # merge GFF3 attributes
    counter = 1
    for pos, lines in entries.items():
        new_line = lines[0].split('\t')[:8]
        # sum coverage; if none is specified, leave as "."
        try:
            new_cov = sum([int(i.split('\t')[5]) for i in lines])
        except ValueError:
            new_cov = '.'
        new_line[5] = new_cov
        # merge attributes
        new_attr = []
        attr_dict = parse_attr(lines)
        attr_dict['ID'] = [str(counter)] # for now, ID will just be a number
        for attr in sorted(attr_dict):
            i = sorted(attr_dict[attr])
            if len(i) > 0:
                j = attr + '=' + ','.join(i)
                new_attr.append(j)
        new_attr = ';'.join(new_attr)
        new_line.append(new_attr)
        new_entries.append(new_line)
        counter += 1

    return header, new_entries


def sort_features(entries):
    """Sort GFF3 formatted entries by chromosome then start position.
    Note: C. elegans uses roman numerals for chromosomes names.
    """
    numerals = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'X':10, 'MtDNA':11}
    try:
        return sorted(entries, key = lambda x: (numerals[x[0]], int(x[3]), int(x[4]), x[6]))
    except KeyError:
        return sorted(entries, key = lambda x: (x[0], int(x[3]), int(x[4]), x[6]))


if __name__ == '__main__':
    file_list = sys.argv[1:] # files specified as arguments
    header, entries = parse_GFF3(file_list)
    entries_sorted = sort_features(entries)

    for i in header:
        print(i)
    for i in entries_sorted:
        print('\t'.join(map(str, i)))
