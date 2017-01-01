#!/usr/bin/env python
'''
Runs the Opal LDPC k-mer hash based metagenomic classifier. Based off the paper
"Low-density locality-sensitive hashing boosts metagenomic binning" by Yunan
Luo, Jianyeng Zeng, Bonnie Berger, and Jian Peng in the conference Recomb 2016.
Journal version yet to appear.

This Python wrapper was written by Yun William Yu <contact@yunwilliamyu.net>,
and is based off an earlier set of prototyping Bash scripts by Yunan Luo.

The implementation of the metagenomic binning is adapted from the source code
of K. Vervier, P. Mahe, M. Tournoud, J.-B. Veyrieras, and J.-P. Vert.
Large-scale Machine Learning for Metagenomics Sequence Classification ,
Technical report HAL-01151453, May, 2015.  This code is included in the util/
directory, with modifications to enable using the Opal Gallagher code based
hashes in util/ldpc.py.

The code from Verview, et al, requires the Genetic Data Analysis Library, which
we have included a copy of under util/ext/ for ease of installation.

This pipeline depends on Python scikit-learn and on Vowpal Wabbit. Vowpal
Wabbit must be properly installed in the system path.
'''
from __future__ import print_function
__version__ = "0.9.0"

import argparse
import os
import sys
import glob
import subprocess
import random
import threading
import pandas as pd
import numpy as np
from sklearn.metrics import precision_score, recall_score
from datetime import datetime

script_loc = os.path.realpath(__file__)
sys.path.append(os.path.join(os.path.dirname(script_loc),'util'))
import ldpc

my_env = os.environ.copy()
my_env["PATH"]=(
        os.path.join(os.path.dirname(script_loc),'util') + ":" + 
        os.path.join(os.path.dirname(script_loc),
                'util','ext','gdl-1.1','GDL','bin') + ":" + 
        os.path.join(os.path.dirname(script_loc),
                'util','ext','gdl-1.1','GDL','include') + ":" +
        my_env.get("PATH", ""))
my_env["LD_LIBRARY_PATH"]=(
        os.path.join(os.path.dirname(script_loc),
                'util','ext','gdl-1.1','GDL','lib') + ":" +
        my_env.get("LD_LIBRARY_PATH", ""))

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def unique_lines(file):
    '''gets number of unique lines in file'''
    seen = set()
    with open(file) as f:
        for line in f:
            seen.add(line)
    return len(seen)

def safe_makedirs(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return 0

def extract_column_two(infile, outfile):
    """cut -f2 infile > outfile"""
    with open(infile, 'r') as inf:
        with open(outfile, 'w') as outf:
            for line in inf:
                parts = line.split()
                if len(parts) > 1:
                    print(parts[1], file=outf)
                else:
                    print('',file=outf)

def vw_class_to_taxid(inputfile, dicofile, outputfile):
    '''Converts vw IDs in a newline delimited list (inputfile) to 
    outputfile using the mapping specified in dicofile'''
    dico = {}
    with open(dicofile, "r") as fin:
        for line in fin:
            txid, vwid = line.strip().split()[:2]
            dico[vwid] = txid
    predout = open(outputfile, "w")
    with open(inputfile, "r") as fin:
        for line in fin:
            predout.write("%s\n"%(dico[str(int(float(line.strip())))]))
    predout.close()

def get_fasta_and_taxid(directory):
    '''finds the 'first' fasta file in directory, and returns a tuple with
    it and the matching named taxid file in the directory if both exist'''
    try:
        fasta = glob.glob(directory + "/*.fasta")[0]
    except IndexError:
        raise RuntimeError("Could not find fasta file in:" + directory)
    taxids = os.path.splitext(fasta)[0] + ".taxid"
    if not os.path.isfile(taxids):
        raise RuntimeError("Could not find matching taxid: " + taxids)
    return [fasta, taxids]

def get_final_model(directory):
    '''gets a 'final' model from a directory. Note, will match the first
    file ending in _final.model'''
    try:
        model = glob.glob(directory + "/*_final.model")[0]
    except IndexError:
        raise RuntimeError("Could not find final model file in: " + directory)
    return model

def evaluate_predictions(reffile, predfile):
    '''Evaluates how good a predicted list is compared to a reference gold standard'''
    with open(predfile, "r") as fin:
        pred = fin.read().splitlines()
    with open(reffile, "r") as fin:
        ref = fin.read().splitlines()

    pred = map(int, pred)
    ref = map(int, ref)
    correct = np.equal(pred, ref)

    perf = pd.DataFrame({"pred":pred, "ref":ref, "correct":correct})
    tmp = perf.groupby("ref")
    species = tmp["correct"].agg(np.mean)
    micro = np.mean(correct)
    macro = np.mean(species)
    median = np.median(species)
    print("micro = {:.2f}".format(micro*100))
    print("macro = {:.2f}".format(macro*100))
    print("median = {:.2f}".format(median*100))
    sys.stdout.flush()


def frag(test_dir, frag_dir, args):
    '''Draws fragments from the fasta file found in test_dir. Note that
    there must be a taxid file of the same basename with matching ids for
    each of the fasta lines.
    
    test_dir (string):  must be a path to a directory with a single fasta
                        and taxid file
    frag_dir (string):  must be a path to an output directory

    Unpacking args:
        frag_length (int):  length of fragments to be drawn
        coverage (float):   fraction of times each location is to be covered
                            by drawn fragments
    '''
    # Unpack args
    frag_length = args.frag_length
    coverage = args.coverage
    # Finish unpacking args

    fasta, taxids = get_fasta_and_taxid(test_dir)
    safe_makedirs(frag_dir)
    fasta_out = os.path.join(frag_dir, "test.fragments.fasta")
    gi2taxid_out = os.path.join(frag_dir, "test.fragments.gi2taxid")
    taxid_out = os.path.join(frag_dir, "test.fragments.taxid")
    starttime = datetime.now()
    print(
    '''================================================
Drawing fragments
{:%Y-%m-%d %H:%M:%S}
'''.format(starttime) + '''
frag_length = {frag_length}
coverage = {coverage}
------------------------------------------------
Fasta input:    {fasta}
taxids input:   {taxids}

Fasta output:   {fasta_out}
gi2taxid output:{gi2taxid_out}
taxids output:  {taxid_out}'''.format(
    frag_length=frag_length, coverage=coverage, fasta=fasta,
    taxids=taxids, fasta_out=fasta_out, gi2taxid_out=gi2taxid_out,
    taxid_out=taxid_out)
    )
    sys.stdout.flush()
    # set seed (for reproducibility)
    seed = 42
    # draw fragments
    subprocess.check_call(["drawfrag",
        "-i", fasta,
        "-t", taxids,
        "-l", str(frag_length),
        "-c", str(coverage),
        "-o", fasta_out,
        "-g", gi2taxid_out,
        "-s", str(seed)],
        env=my_env)

    # extract taxids
    extract_column_two(gi2taxid_out, taxid_out)
    print('''------------------------------------------------
Total wall clock runtime (sec): {}
================================================'''.format(
    (datetime.now() - starttime).total_seconds()))
    sys.stdout.flush()

    return 0

def train(ref_dir, model_dir, args):
    '''Draws fragments from the fasta file found in data_dir. Note that
    there must be a taxid file of the same basename with matching ids for
    each of the fasta lines.
    
    ref_dir (string):   must be a path to a directory with a single fasta
                        and taxid file
    model_dir (string): must be a path to an output directory

    Unpacking args:
        frag_length (int):  length of fragments to be drawn
        coverage (float):   fraction of times each location is to be covered
                            by drawn fragments
        kmer (int):         size of k-mers used
        row_weight (int):   how many positions will be randomly chosen in the
                            contiguous k-mer (k-mer length should be multiple
                            of row_weight)

        num_hash (int):     number of hashing functions
        num_batches (int):  number of times to run vowpal_wabbit
        num_passes (int):   number of passes within vowpal_wabbit
    '''
    # Unpack args
    frag_length = args.frag_length
    coverage = args.coverage
    kmer = args.kmer
    row_weight = args.row_weight
    hierarchical = args.hierarchical_weight # only comes into play if > 0
    num_hash = args.num_hash
    num_batches = args.num_batches
    num_passes = args.num_passes
    bits = args.bits
    lambda1 = args.lambda1
    lambda2 = args.lambda2
    # Finish unpacking args

    fasta, taxids = get_fasta_and_taxid(ref_dir)
    starttime = datetime.now()

    if kmer % row_weight != 0:
        raise ValueError("Row weight [{}] must divide into k-mer length [{}].".format(row_weight, kmer))
    if (hierarchical > 0):
        if kmer % hierarchical != 0:
            raise ValueError("Hierarchy middle level [{}] must divide into k-mer length [{}].".format(hierarchical, kmer))
        if hierarchical % row_weight != 0:
            raise ValueError("Row weight[{}] must divide into middle hierarchical structure weight [{}].".format(row_weight, hierarchical))

    print(
    '''================================================
Training using Opal + vowpal-wabbit
{:%Y-%m-%d %H:%M:%S}
'''.format(starttime) + '''
frag_length = {frag_length}
coverage:       {coverage}
k-mer length:   {kmer}
row weight:     {row_weight}'''.format(
    frag_length=frag_length,
    coverage=coverage,
    kmer=kmer,
    row_weight=row_weight))
    if hierarchical > 0:
        print('''hierarchical:   {}'''.format(hierarchical))
    print('''num hashes:     {num_hash}
num batches:    {num_batches}
num passes:     {num_passes}
------------------------------------------------
Fasta input:    {fasta}
taxids input:   {taxids}
------------------------------------------------'''.format(
    num_hash=num_hash,
    num_batches=num_batches,
    num_passes=num_passes,
    fasta=fasta,
    taxids=taxids)
    )
    sys.stdout.flush()
    num_labels = unique_lines(taxids)
    print("Number labels:  {}".format(num_labels))
    sys.stdout.flush()

    safe_makedirs(model_dir)

    # define output "dictionary" : taxid <--> vw classes
    dico = os.path.join(model_dir, "vw-dico.txt")
    
    # define model prefix
    model_prefix = os.path.join(model_dir, "vw-model")

    # generate LDPC spaced pattern
    pattern_file = os.path.join(model_dir, "patterns.txt")
    ldpc.ldpc_write(k=kmer, t=row_weight, _m=num_hash, d=pattern_file)

    seed = 42
    for i in range(num_batches):
        seed = seed + 1
        batch_prefix = os.path.join(model_dir, "train.batch-{}".format(i))
        fasta_batch = batch_prefix + ".fasta"
        gi2taxid_batch = batch_prefix + ".gi2taxid"
        taxid_batch = batch_prefix + ".taxid"

        # draw fragments
        subprocess.check_call(["drawfrag",
            "-i", fasta,
            "-t", taxids,
            "-l", str(frag_length),
            "-c", str(coverage),
            "-o", fasta_batch,
            "-g", gi2taxid_batch,
            "-s", str(seed)],
            env=my_env)
        # extract taxids
        extract_column_two(gi2taxid_batch, taxid_batch)

        # learn model
        fasta2skm_param_list = ["fasta2skm",
            "-i", fasta_batch,
            "-t", taxid_batch,
            "-k", str(kmer),
            "-d", dico,
            "-p", pattern_file]
        print("Getting training set ...")
        sys.stdout.flush()
        training_list = subprocess.check_output(
                fasta2skm_param_list, env=my_env).splitlines()
        print("Shuffling training set ...")
        sys.stdout.flush()
        random.shuffle(training_list)
        curr_model = model_prefix + "_batch-{}.model".format(i)
        prev_model = model_prefix + "_batch-{}.model".format(i-1) # May not exist if first run
        vw_param_base = ["vw",
            "--random_seed", str(seed),
            "-f", curr_model,
            "--cache_file", batch_prefix + ".cache",
            "--passes", str(num_passes),
            "--save_resume"]
        vw_param_firstrun = [
            "--oaa", str(num_labels),
            "--bit_precision", str(bits),
            "--l1", str(lambda1),
            "--l2", str(lambda2)]
        if i > 0:
            vw_param_list = vw_param_base + ["-i", prev_model]
        else:
            vw_param_list = vw_param_base + vw_param_firstrun
        print(vw_param_list)
        sys.stdout.flush()
        vwps = subprocess.Popen(vw_param_list, env=my_env,
                stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)
        # Cannot use
        #   vwps.communicate(input='\n'.join(training_list))
        # Because I want to display output as it comes, so have to
        # manually to threading. -ywy 2016-12-28
        def vw_pipe_writer():
            vwps.stdin.write('\n'.join(training_list))
            vwps.stdin.close()
        thread = threading.Thread(target=vw_pipe_writer)
        thread.start()
        while vwps.poll() is None:
            l = vwps.stdout.readline()
            sys.stdout.write(l)
            sys.stdout.flush()
        thread.join() # This shouldn't be necessary, but just being safe.

        if i > 0:
            os.remove(prev_model)
        if i == num_batches - 1:
            os.rename(curr_model, model_prefix + "_final.model")
        os.remove(batch_prefix + ".cache")
        os.remove(fasta_batch)
        os.remove(taxid_batch)
        os.remove(gi2taxid_batch)
    print('''------------------------------------------------
Total wall clock runtime (sec): {}
================================================'''.format(
    (datetime.now() - starttime).total_seconds()))
    sys.stdout.flush()
    return 0


def predict(model_dir, test_dir, predict_dir, args):
    '''Draws fragments from the fasta file found in data_dir. Note that
    there must be a taxid file of the same basename with matching ids for
    each of the fasta lines.
    
    ref_dir (string):   must be a path to a directory with a single fasta
                        and taxid file
    model_dir (string): must be a path to a directory with a vw model file
    predict_dir (string):output directory of predictions

    Unpacking args:
        kmer (int):         size of k-mers used

    Returns a tuple with (reffile, predicted_labels_file) for easy input
    into evaluate_predictions.
    '''
    # Unpack args
    kmer = args.kmer
    # Finish unpacking args

    # Don't need to get taxids until eval
    #fasta, taxids = get_fasta_and_taxid(test_dir)
    try:
        fasta = glob.glob(test_dir + "/*.fasta")[0]
    except:
        raise RuntimeError("No fasta file found in: " + test_dir)
    model = get_final_model(model_dir)
    dico = os.path.join(model_dir, "vw-dico.txt")
    pattern_file = os.path.join(model_dir, "patterns.txt")
    starttime = datetime.now()
    print(
    '''================================================
Predicting using Opal + vowpal-wabbit
{:%Y-%m-%d %H:%M:%S}
'''.format(starttime) + '''
k-mer length:   {kmer}
------------------------------------------------
Fasta input:    {fasta}
Model used:     {model}
Dict used:      {dico}
LDPC patterns:  {pattern_file}
------------------------------------------------'''.format(
    kmer=kmer,
    fasta=fasta,
    model=model,
    dico=dico,
    pattern_file=pattern_file)
    )
    sys.stdout.flush()
    safe_makedirs(predict_dir)
    prefix = os.path.join(predict_dir, "test.fragments-db")

    # get vw predictions
    fasta2skm_param_list = ["fasta2skm",
        "-i", fasta,
        "-k", str(kmer),
        "-p", pattern_file]
    vw_param_list = ["vw", "-t",
        "-i", model,
        "-p", prefix + ".preds.vw"]
    ps = subprocess.Popen(fasta2skm_param_list, env=my_env, 
            stdout=subprocess.PIPE)
    vwps = subprocess.Popen(vw_param_list, env=my_env,
            stdin=ps.stdout, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
    while vwps.poll() is None:
        l = vwps.stdout.readline()
        sys.stdout.write(l)
        sys.stdout.flush()

    # Convert back to standard taxonomic IDs instead of IDs
    vw_class_to_taxid(prefix + '.preds.vw', dico, prefix + '.preds.taxid')

    print('''------------------------------------------------
Predicted labels:   {pl}
Total wall clock runtime (sec): {s}
================================================'''.format(
    pl=prefix + '.preds.taxid',
    s=(datetime.now() - starttime).total_seconds()))
    sys.stdout.flush()
    return (prefix + '.preds.taxid')


def parse_extra(parser, namespace):
    namespaces = []
    extra = namespace.extra
    while extra:
        n = parser.parse_args(extra)
        extra = n.extra
        namespaces.append(n)
    return namespaces

class ArgClass:
    '''So I don't have to duplicate argument info'''
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

def main(argv):
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description=__doc__)
    parser.add_argument('--version', action='version',
            version='%(prog)s {version}'.format(version=__version__))

    # Shared arguments
    frag_length_arg = ArgClass("-l", "--frag-length",
            help="length of fragments to be drawn from fasta",
            type=int, default=16)
    kmer_arg = ArgClass("-k", "--kmer", help="length of k-mers used",
            type=int, default=8)
    coverage_arg = ArgClass("-c", "--coverage", help="""number/fraction of
            times each location in a fragment should be covered by a k-mer""",
            type=float, default=1.0)
    hierarchical_arg = ArgClass("--hierarchical-weight",
            help="intermediate organization of positions chosen in the k-mer in row_weight; should be a multiple of row_weight and a divisor of k-mer length if set", type=int, default=-1)
    row_weight_arg = ArgClass("--row-weight", help="""the number of positions
            that will be randomly chosen in the contiguous k-mer; k-mer
            length should be a multiple of row_weight""", type=int, default=4)
    num_hash_arg = ArgClass("--num-hash", help="""number of k-mer hashing
            functions to get features""", type=int, default=1)
    num_batches_arg = ArgClass("--num-batches", help="""Number of times to
            generate a random batch of training data for VW""",
            type=int, default=1)
    num_passes_arg = ArgClass("--num-passes",
            help="Number of VW passes in each training batch",
            type=int, default=1)
    bits_arg = ArgClass("--bits", help="Number of bits used in VW model",
            type=int, default=31)
    lambda1_arg = ArgClass("--lambda1", help="VW model lambda1 training parameter", type=float, default=0.)
    lambda2_arg = ArgClass("--lambda2", help="VW model lambda2 training parameter", type=float, default=0.)


    subparsers = parser.add_subparsers(help="sub-commands", dest="mode")

    parser_frag = subparsers.add_parser("frag", help="Fragment a fasta file into substrings for training/testing",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_frag.add_argument("test_dir", help="Input directory for test data")
    parser_frag.add_argument("frag_dir", help="Output directory for fasta fragments")
    parser_frag.add_argument(*frag_length_arg.args, **frag_length_arg.kwargs)
    parser_frag.add_argument(*coverage_arg.args, **coverage_arg.kwargs)

    parser_train = subparsers.add_parser("train", help="Train a Vowpal Wabbit model using Opal hashes",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_train.add_argument("train_dir", help="Input directory for train data")
    parser_train.add_argument("model_dir", help="Output directory for VW model")
    parser_train.add_argument(*frag_length_arg.args, **frag_length_arg.kwargs)
    parser_train.add_argument(*coverage_arg.args, **coverage_arg.kwargs)
    parser_train.add_argument(*kmer_arg.args, **kmer_arg.kwargs)
    parser_train.add_argument(*num_batches_arg.args, **num_batches_arg.kwargs)
    parser_train.add_argument(*num_passes_arg.args, **num_passes_arg.kwargs)
    parser_train.add_argument(*num_hash_arg.args, **num_hash_arg.kwargs)
    parser_train.add_argument(*row_weight_arg.args, **row_weight_arg.kwargs)
    parser_train.add_argument(*hierarchical_arg.args, **hierarchical_arg.kwargs)
    parser_train.add_argument(*bits_arg.args, **bits_arg.kwargs)
    parser_train.add_argument(*lambda1_arg.args, **lambda1_arg.kwargs)
    parser_train.add_argument(*lambda2_arg.args, **lambda2_arg.kwargs)

    parser_predict = subparsers.add_parser("predict", help="Predict metagenomic classifications given a Opal/VW model",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_predict.add_argument("model_dir", help="Input directory for VW model")
    parser_predict.add_argument("test_dir", help="Input directory for already fragmented test data")
    parser_predict.add_argument("predict_dir", help="Output directory for predictions")
    parser_predict.add_argument(*kmer_arg.args, **kmer_arg.kwargs)

    parser_eval = subparsers.add_parser('eval', help="Evaluate quality of predictions given a reference",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_eval.add_argument("reference_file", help="Gold standard labels")
    parser_eval.add_argument("predicted_labels", help="Predicted labels")

    parser_simulate = subparsers.add_parser('simulate', help=
    '''Run a full pipeline of frag, train, predict, and eval to
determine how good a model is under particular parameter
ranges''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_simulate.add_argument("test_dir", help="Input directory for test data")
    parser_simulate.add_argument("train_dir", help="Input directory for train data")
    parser_simulate.add_argument("out_dir", help="Output directory for all steps")
    parser_simulate.add_argument("--do-not-fragment", help="If set, will use test_dir fasta files as is without fragmenting", action="store_true")
    parser_simulate.add_argument(*frag_length_arg.args, **frag_length_arg.kwargs)
    parser_simulate.add_argument(*coverage_arg.args, **coverage_arg.kwargs)
    parser_simulate.add_argument(*kmer_arg.args, **kmer_arg.kwargs)
    parser_simulate.add_argument(*num_batches_arg.args, **num_batches_arg.kwargs)
    parser_simulate.add_argument(*num_passes_arg.args, **num_passes_arg.kwargs)
    parser_simulate.add_argument(*num_hash_arg.args, **num_hash_arg.kwargs)
    parser_simulate.add_argument(*row_weight_arg.args, **row_weight_arg.kwargs)
    parser_simulate.add_argument(*hierarchical_arg.args, **hierarchical_arg.kwargs)
    parser_simulate.add_argument(*bits_arg.args, **bits_arg.kwargs)
    parser_simulate.add_argument(*lambda1_arg.args, **lambda1_arg.kwargs)
    parser_simulate.add_argument(*lambda2_arg.args, **lambda2_arg.kwargs)

    args = parser.parse_args(argv)

    print(args)
    sys.stdout.flush()

    mode = args.mode
    if (mode == "simulate"):
        fullstarttime = datetime.now()
        print("Full simulation")
        print("{:%Y-%m-%d %H:%M:%S}".format(fullstarttime))
        print("Fragment mode: {}".format(not args.do_not_fragment))
        output_dir = args.out_dir
        frag_dir = os.path.join(output_dir, '1frag')
        model_dir = os.path.join(output_dir, '2model')
        predict_dir = os.path.join(output_dir, '3predict')
        if args.do_not_fragment:
            train(args.train_dir, model_dir, args)
            pf = predict(model_dir, args.test_dir, predict_dir, args)
            _, rf = get_fasta_and_taxid(args.test_dir)
        else:
            frag(args.test_dir, frag_dir, args)
            train(args.train_dir, model_dir, args)
            pf = predict(model_dir, frag_dir, predict_dir, args)
            _, rf = get_fasta_and_taxid(frag_dir)

        print("Evaluation reference file: " + rf)
        sys.stdout.flush()
        evaluate_predictions(rf, pf)
        print("Total full sim wall clock runtime (sec): {}".format(
            (datetime.now() - fullstarttime).total_seconds()))
    elif mode == "frag":
        frag(args.test_dir, args.frag_dir, args)
    elif mode == "train":
        train(args.train_dir, args.model_dir, args)
    elif mode == "predict":
        predict(args.model_dir, args.test_dir, args.predict_dir, args)
    elif mode == "eval":
        evaluate_predictions(args.reference_file, args.predicted_labels)

if __name__ == "__main__":
    main(sys.argv[1:])

