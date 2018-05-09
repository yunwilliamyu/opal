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

def fasta_reader(f):
    '''Generator expression that returns a fasta sequence
        
        Ignores quality score string of FASTQ file
    '''
    seq = ''
    name = ''
    ignore_line = False
    first_line = True
    while True:
        line = f.readline()
        if line=='':
            if seq=='':
                break
            else:
                yield (name, seq)
                name = ''
                seq = ''
        elif line[0]=='>':
            if first_line:
                first_line = False
                pass
            else:
                yield (name, seq)
            name = line[1:].rstrip('\n')
            seq = ''
            ignore_line = False
        elif line[0]=='+':
            # Ignore quality score strings
            ignore_line = True
        else:
            if ignore_line:
                pass
            else:
                seq = seq + line.rstrip('\n')
    

def main(argv):
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description=__doc__)
    parser.add_argument('--version', action='version',
            version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('-i', '--input', help='file containing input sequences (fasta or fastq) [required]', type=argparse.FileType('r'))
    parser.add_argument('-k', '--kmer', help='k-mer size [required]', type=int)
    parser.add_argument('-t', '--taxid', help='one-column file containing the taxid of each input sequence] (if not specified, instances are labeled as 1)', type=argparse.FileType('r'))
    parser.add_argument('-p', '--pattern', help='LDPC code locations. Space delimited. First line specifies number of hashes and row_weight. Following lines specify which locations are selected within the k-mer.', type=argparse.FileType('r'))
    parser.add_argument('-d', '--dico', help='''two-column file giving the correspondance between taxids and VW class labels]
        (created if does not exist; possibly updated and overwritten if exists)
        (must be specified if --taxid option is used)
        (useless if --taxid option is not used)''')
    parser.add_argument('-o', '--output', help='output file', type=argparse.FileType('w'), default='-')
    
    args = parser.parse_args(argv)
    if not args.input or not args.kmer:
        parser.print_help()
        exit(-1)

    # Updates the dictionary file converting labels to vwid
    label2vwid = {}
    vwidset = set()
    if args.taxid:
        if not args.dico:
            parser.print_help()
            exit(-1)
        if os.path.isfile(args.dico):
            with open(args.dico, "r") as df:
                for line in df.readlines():
                    label, vwid_str = line.strip().split()[:2]
                    vwid = int(vwid_str)
                    label2vwid[label] = vwid
                    vwidset.add(vwid)
        for line in args.taxid:
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
        with open(args.dico, "w") as df:
            for label, vwid in label2vwid.items():
                df.write("{}\t{}\n".format(label, vwid))

    # Reads in the pattern file
    pattern_list = []
    if args.pattern:
        num_hash, row_weight = [int(x) for x in args.pattern.readline().split()[:2]]
        assert(args.kmer%row_weight ==0)
        for _ in range(num_hash):
            row = [int(x) for x in args.pattern.readline().split()]
            assert(len(row)==row_weight)
            pattern_list.append(row)
    else:
        row = [x for x in range(args.kmer)]
        pattern_list.append(row)
    pattern_getters = [operator.itemgetter(*pl) for pl in pattern_list]
    
    args.taxid.seek(0)
    for _, seq in fasta_reader(args.input):
        label = args.taxid.readline().rstrip('\n')
        #feature_list = gen_features(pattern_list, seq, args.kmer)
        feature_list = gen_features2(pattern_getters, seq, args.kmer)
        features = " ".join(feature_list)
        args.output.write('{} | {}\n'.format(label2vwid[label], features))


def get_all_substrings(input_string, k):
    length = len(input_string)
    return [input_string[i:i+k] for i in range(length - k + 1)]

def gen_features2(pattern_getters, seq, k):
    '''Generates featrues from a pattern list and a sequence'''
    kmers = get_all_substrings(seq, k)
    feature_list = ["".join(pat(kmer))+str(i) for kmer in kmers for i, pat in enumerate(pattern_getters)]
    return feature_list

if __name__=="__main__":
    main(sys.argv[1:])
