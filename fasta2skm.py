#!/usr/bin/env python
'''Replaces fasta2skm C program

Generates vowpal_wabbit compatible features based off of Opal LDPC hashing.

    The C program has several bugs, including not accepting single character labels, and not being able to read non-numeric labels from the dico file, though it can generate them.
'''

from __future__ import print_function
__version__ = "0.0.1"
import argparse
import os
import sys
import operator
import itertools

from fasta_functions import fasta_reader, reverse_complement, get_all_substrings

def update_dictionary(labels, dico_file):
    '''Updates the dictionary file converting labels to vwid with an iterator over new labels'''
    label2vwid = {}
    vwidset = set()
    if os.path.isfile(dico_file):
        with open(dico_file, "r") as df:
            for line in df.readlines():
                label, vwid_str = line.strip().split()[:2]
                vwid = int(vwid_str)
                label2vwid[label] = vwid
                vwidset.add(vwid)
    for line in labels:
        label = line.rstrip('\n')
        if label in label2vwid:
            pass
        else:
            if len(vwidset)==0:
                vwid = 1
            else:
                vwid = max(vwidset)+1
            label2vwid[label] = vwid
            vwidset.add(vwid)
    with open(dico_file, "w") as df:
        for label, vwid in label2vwid.items():
            df.write("{}\t{}\n".format(label, vwid))
    return label2vwid

def main_not_commandline(args):
    '''All the main code except for the parser'''
    if not args.input or not args.kmer:
        raise ValueError("fasta2skm requires input and kmer arguments")

    # Updates the dictionary file converting labels to vwid
    if args.taxid:
        if not args.dico:
            raise ValueError("If --taxid is set, then so must be --dico")
        with open(args.taxid, 'r') as taxid_file:
            label2vwid = update_dictionary(taxid_file, args.dico)

    # Reads in the pattern file
    if args.pattern:
        with open(args.pattern, 'r') as pattern_file:
            file_contents = pattern_file.readlines()
    else:
        file_contents = None
    pattern_getters = create_pattern_getters(file_contents, args.kmer)
    
    if args.taxid:
        taxid_file = open(args.taxid, 'r')
        labels = (label2vwid[l.rstrip('\n')] for l in taxid_file)
    else:
        labels = itertools.repeat(1)

    with open(args.input, 'r') as input_file:
        for _, seq in fasta_reader(input_file):
            feature_list = gen_features(pattern_getters, seq, args.kmer)
            if args.reverse:
                feature_list.extend(gen_features(pattern_getters, reverse_complement(seq), args.kmer))
            features = " ".join(feature_list)
            args.output.write('{} | {}\n'.format(labels.next(), features))
    args.output.close()

    if args.taxid:
        taxid_file.close()

def main(argv):
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description=__doc__)
    parser.add_argument('--version', action='version',
            version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('-i', '--input', help='file containing input sequences (fasta or fastq) [required]')
    parser.add_argument('-k', '--kmer', help='k-mer size [required]', type=int)
    parser.add_argument('-t', '--taxid', help='one-column file containing the taxid of each input sequence] (if not specified, instances are labeled as 1)')
    parser.add_argument('-p', '--pattern', help='LDPC code locations. Space delimited. First line specifies number of hashes and row_weight. Following lines specify which locations are selected within the k-mer.')
    parser.add_argument('-d', '--dico', help='''two-column file giving the correspondance between taxids and VW class labels]
        (created if does not exist; possibly updated and overwritten if exists)
        (must be specified if --taxid option is used)
        (useless if --taxid option is not used)''')
    parser.add_argument('-o', '--output', help='output file', type=argparse.FileType('w'), default='-')
    parser.add_argument('-r', '--reverse', help='Take the reverse complements of sequences; fails if non-ACGT sequences provided', action='store_true')
    
    args = parser.parse_args(argv)
    print(args)
    main_not_commandline(args)


def create_pattern_getters(pattern_file_contents, kmer):
    '''Reads in the pattern file'''
    pattern_list = []
    if pattern_file_contents:
        num_hash, row_weight = [int(x) for x in pattern_file_contents[0].split()[:2]]
        assert(kmer%row_weight ==0)
        for i in range(num_hash):
            row = [int(x) for x in pattern_file_contents[1+i].split()]
            assert(len(row)==row_weight)
            pattern_list.append(row)
    else:
        row = [x for x in range(kmer)]
        pattern_list.append(row)
    pattern_getters = [operator.itemgetter(*pl) for pl in pattern_list]
    return pattern_getters

def gen_features(pattern_getters, seq, k):
    '''Generates features from a pattern list and a sequence'''
    kmers = get_all_substrings(seq, k)
    feature_list = ["".join(pat(kmer))+str(i) for kmer in kmers for i, pat in enumerate(pattern_getters)]
    return feature_list

if __name__=="__main__":
    main(sys.argv[1:])
