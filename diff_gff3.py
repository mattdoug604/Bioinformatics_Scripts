#!/home2/mattdoug/bin/python3
#last updated: 22/3/2018

# PURPOSE: Compare multiple GFF3 formatted files and list the features unique
#          to each, and features shared by all.
# USAGE:   diff_gff3.py -i features_1.gff3 features_2.gff3 ... features_N.gff3 -o output_dir -t feature_type

import argparse, errno, operator, os, shutil, sys
from collections import defaultdict

def print_to_log(*args, **kwargs):
    """Print to both stdout and the log file."""
    print(*args, file=sys.stderr, **kwargs)
    print(*args, file=log, **kwargs)


def parse_commandline_arguments():
    dirname = lambda x: os.path.splitext(os.path.basename(x))[0]

    parser = argparse.ArgumentParser(description='Compare features by position in two or more GFF3 files.')
    parser.add_argument('-i',
                        type=str,
                        nargs='+',
                        help='two or more files in GFF3 format')
    parser.add_argument('-o',
                        type=str,
                        nargs='?',
                        help='output directory to create')
    parser.add_argument('-t',
                        type=str,
                        nargs='+',
                        help="[optional] only count features matching 'type' (column 3 in the file)")
    parser.add_argument('-f',
                        action='store_true',
                        help='overwrite the output directory if it already exists (THIS WILL OVERWRITE FILES OF THE SAME NAME)')
    args = parser.parse_args()

    # make sure 2 or more input files are specified
    if args.i == None:
        parser.print_help()
        sys.exit(1)
    elif len(args.i) < 2:
        print('Two or more input files are required! Exiting.')
        sys.exit(1)

    # make sure all the input files exist before starting
    for i in args.i:
        if not os.path.exists(i):
            print('One of more input files could not be found! Exiting.')
            sys.exit(1)

    # default output directory is just <name of this script> + "_output"
    if args.o is None:
        args.o = '{}_output'.format(dirname(sys.argv[0]))

    return args


def file_pairs(files):
    """Return all possible pairs of elements from a list."""
    result = []
    for f1 in range(len(files)):
        for f2 in range(f1+1, len(files)):
            result.append([files[f1], files[f2]])

    return result


def sort_features(pos):
    """Sort GFF3 formatted entries by chromosome then start position.
    Note: C. elegans uses roman numerals for chromosomes names.
    """
    numerals = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'X':10, 'MtDNA':11}
    try:
        return sorted(pos, key = lambda x: (numerals[x[0]], x[1], x[2], x[3]))
    except KeyError:
        return sorted(pos, key = lambda x: (x[0], x[1], x[2], x[3]))


def compare_features(files, target=None):
    """Compare sets of features and find those unique to each set."""
    index = {f:{} for f in files}
    features = {f:set() for f in files}
    unique = {}
    paired = {}
    common = []

    if target is not None:
        print('Searching for features matching: {}'.format(' or '.join(target)))
        target = [i.strip() for i in target]

    for x, f in enumerate(files):
        print_to_log('File #{} = {}'.format(x+1, os.path.abspath(f)))
        with open(f, 'r') as infile:
            for n, line in enumerate(infile):
                line = line.strip()
                if line[0] != '#':
                    i = line.split('\t')
                    f_pos = i[0], int(i[3]), int(i[4]), i[6]
                    f_type = i[2]
                    if target is not None and f_type not in target:
                        continue
                    features[f].add(f_pos)
                    index[f][f_pos] = line # Note: if more than one feature
                                           # shares the same position, only the
                                           # last feature appearing will be
                                           # counted.
        print_to_log("  {:,} features found".format(len(features[f])))

    # list features unique to each file
    for f in files:
        f_unique = features[f].copy()
        for f2 in [i for i in files if i != f]:
            f_unique -= features[f2]
        f_unique = sort_features(f_unique) # sort the features by postiion
        unique[f] = [index[f][i] for i in f_unique] # return the line index for each unique feature

    # get features common to each pair
    if 2 < len(files) < 4:
        for f1, f2 in file_pairs(files):
            pair_common = set.intersection(features[f1], features[f2])
            pair_common = sort_features(pair_common)
            pair_common = [index[f1][i] for i in pair_common] # just use the lines as they appear in the first file
            paired[(f1, f2)] = pair_common

    # get the features common to all files
    common = set.intersection(*[features[f] for f in files])
    common = sort_features(common)
    common = [index[files[0]][i] for i in common] # just use the lines as they appear in the first file

    return unique, paired, common


if __name__ == '__main__':
    args = parse_commandline_arguments()
    files = args.i
    out_path = args.o
    f_type = args.t
    force = args.f

    # Make a directory to output the results to
    if not os.path.isdir(out_path):
        os.makedirs(out_path)
    else:
        if force:
            print("WARNING: Output directory '{}' already exists! Some files may be overwritten...".format(out_path))
        else:
            print("WARNING: Output directory '{}' already exists! Exiting.".format(out_path))
            sys.exit(1)

    log = open(out_path + '/files.log', 'w')
    unique, paired, common = compare_features(files, f_type)

    # output the features common to all files
    print_to_log('{:,} features common to all files'.format(len(common)))
    with open(out_path + '/features_common_to_all.gff3', 'w') as outf:
        print('##gff-version 3', file=outf)
        for line in common:
            print(line, file=outf)

    # output the results common to each pair of files
    if len(paired) > 0:
        for f1, f2 in sorted(paired):
            x1, x2 = files.index(f1)+1, files.index(f2)+1
            lines = paired[(f1, f2)]
            print_to_log('  {:,} features common to files #{} and #{}'.format(len(lines), x1, x2))
            with open(out_path + '/features_common_to_{}_and_{}.gff3'.format(x1, x2), 'w') as outf:
                print('##gff-version 3', file=outf)
                for line in lines:
                    print(line, file=outf)

    # output the features unique to each files
    for x, f in enumerate(files):
        print_to_log("{:,} features unique to file #{}".format(len(unique[f]), x+1))
        with open(out_path + '/features_unique_to_{}.gff3'.format(x+1), 'w') as outf:
            print('##gff-version 3', file=outf)
            for line in unique[f]:
                print(line, file=outf)

    print_to_log('Wrote output to', os.path.abspath(out_path))
    log.close()
