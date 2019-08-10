#!/usr/bin/env python3
'''Replaces drawfrag C program

From each Fasta sequence in a file, draw random substrings of size k covering it c times, returning a new multi-fasta file with labels
'''

__version__ = "0.0.2"
import argparse
import sys
import random

from fasta_functions import fasta_reader, check_acgt


def main_not_commandline(args):
    '''All the main code except for the parser'''
    input_file = open(args.input, 'r')
    output_file = open(args.output, 'w')
    taxid_infile = open(args.taxids, 'r')
    gi2taxid_outfile = open(args.gi2taxid, 'w')
    k = args.size

    if args.seed:
        random.seed(args.seed)

    read_num = 0
    for name, seq in fasta_reader(input_file):
        tlabel = taxid_infile.readline().rstrip('\n')
        firstname = name.split()[0]
        coverage = 0
        desired_coverage = args.coverage * len(seq)
        if len(seq) < k:
            pass
        else:
            try_num = 0
            while coverage < desired_coverage:
                try_num = try_num + 1
                pos = random.randint(0, len(seq) - k)
                sample = seq[pos:pos + k]
                if args.atgc and not check_acgt(sample):
                    pass
                else:
                    coverage = coverage + k
                    read_num = read_num + 1
                    output_file.write(">{}\n".format(read_num))
                    output_file.write("{}\n".format(sample))
                    gi2taxid_outfile.write("{}\t{}\n".format(firstname, tlabel))
                if try_num > 10 * len(seq):
                    break


def main(argv):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__)
    parser.add_argument(
        '--version', action='version',
        version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('-i', '--input', help='file containing input sequences (fasta or fastq) [required]')
    parser.add_argument('-l', '--size', help='size of drawn items [required]', type=int)
    parser.add_argument('-t', '--taxids', help='one-column file containing the taxid of each input sequence] [required]')
    parser.add_argument('-c', '--coverage', help='mean coverage value for drawing fragments [required]', type=float)
    parser.add_argument('-g', '--gi2taxid', help='output gi2taxids file: two-column file containing genome ids and taxids of the drawn fragments [required]')
    parser.add_argument('-s', '--seed', help='value used to initialize the random seed (to use for reproducibility purposes; if not set, will be randomly initialized by Python', type=int)
    parser.add_argument('-o', '--output', help='output sequence file [required]')
    parser.add_argument('--atgc', help='draw fragments made of ATCG only', action='store_true')

    args = parser.parse_args(argv)
    print(args)
    main_not_commandline(args)


if __name__ == "__main__":
    main(sys.argv[1:])
