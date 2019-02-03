This code is associated with the following manuscript:
Metagenomic binning through low density hashing. Yunan Luo, Yun William Yu, Jianyang Zeng, Bonnie Berger, and Jian Peng. Submitted for publication.

Further information can be found at http://opal.csail.mit.edu/

0. Requirments
    Vowpal Wabbit >= 8.1.1 (the command-line program must be installed)

    Python 2.7 (this code fails on Python 3)
    Python packages:
        sklearn
        pandas

    This code has been tested on Ubuntu 14.04 and 16.04, running
    under Bash 4.3.11.

    Additionally, while RAM requirements will vary by model size, we recommend
    at least 32GiB of RAM for running with default options.

1. Directory structure
data/: training and testing data should be given a subfolder here
    A1/: a small size example.
        train/: fasta and tax id for training model, 
                10 fasta records of 10 species.
        test/: short fragments and corresponding tax id for testing model.
    example/: a large size example.
        train/: fasta and tax id for training model, 
                1564 fasta records (different strains) of huandreds of species.
        test/: fasta and tax id for testing model.
                Note that it is fasta file of sequence, not short fragments. 
                We need to first simulated test data by randomly drawing fragments 
                from these fasta records (see below).
util/
    drawfrag.py: draw fragments from fasta records.
    fasta2skm.py: construct feature (spaced k-mer profile), and convert to VW input format.
    ldpc.py: generate LSH function using LDPC code.
    fasta_functions.py: parse FASTA files

2. Install and test:
    You must ensure that all dependencies are installed. On Ubuntu, you can do:
        sudo apt-get install vowpal-wabbit python-pip
        pip install pandas, sklearn
    Note that you must install the command-line version of vowpal-wabbit, as Opal depends on the "vw" command.

    Altenrately, the provided  setup script will test your environment for
    dependencies and download some example data files to play with.
    It will fail when it does not find a required dependency, and give
    instructions for how to install them.
        bash SETUP.sh

    To test Opal after setup if you downloaded the example data files using SETUP.sh:
        ./opal.py simulate data/A1/test data/A1/train out_dir
    Note that this may take a very long time to train, so you may change the coverage
    to have a quicker (but horrendously inaccurate test of the pipeline as follows):
        ./opal.py simulate data/A1/test data/A1/train out_dir -c0.1
    (default is -c15, for a coverage of 15)

    The simulate command is equivalent to the following mult-step process with
    an optional -c0.1 for fast but inaccurate results. (default is -c15)
    Create output directory:
        mkdir out_dir

    Fragmentation to prepare the fragmented testing data (not used for training):
        ./opal.py frag data/A1/test out_dir/1frag [-c0.1]
    
    Training based on the training directory fasta:
        ./opal.py train data/A1/train out_dir/2model [-c0.1]

    Predicting examples from the fasta file found in the test directory:
        ./opal.py predict out_dir/2model data/A1/test out_dir/3predict

    Evaulating the accuracy of the predicted assignments:
        ./opal.py eval data/A1/test/A1.test.taxid out_dir/3predict/test.fragments-db.preds.taxid

    As an aside, if you are classifying DNA, you will probably want the "-r" argument on everything
    to also treat reverse complements as the same. We do not use it by default for backwards compatability
    since some users use Opal on protein sequences.

3. Usage: opal.py assumes it lives in the current directory structure, but can be symlinked elsewhere.

Modes:
    (default --optional-arguments such as k-mer length, fragment size,
    hash functions, etc. are set a single batch, and so will use too much
    RAM for extremely large training sets. For training larger data sets,
    be sure to set num-batches and coverage per batch.)
    (Also, use "-r" to enable reverse complements for ACGT genomic sequence
    data. Otherwise, the sequence is treated as simple text.)

    1) ./opal.py frag [--optional-arguments] test_dir frag_dir [-h]

        Looks for a fasta file in test_dir with matching taxid file.
        Randomly draws fragments of length and coverage specified in
        optional-arguments. (use "./opal.py frag -h" for details)

        Outputs these fragments with corresponding taxid into frag_dir.
        
    2) ./opal.py train [--optional-arguments] train_dir model_dir [-h]

        Looks for a fasta file in train_dir with matching taxid file.
        For each batch of training, randomly draw fragments and generate
        feature vectors using Opal LDPC hashes, and trains Vowpal_Wabbit
        One-Against-All classifier against all batches sequentially.

        Outputs the generated classifier model into model_dir.

    3) ./opal.py predict [--optional-arguments] model_dir test_dir predict_dir [-h]

        Looks for a classifier model in model_dir, and a fasta file in
        test_dir. Note that if the fasta file is of a full refence, you
        may need to run ./opal.py frag and use the outputted frag_dir in
        place of test_dir so that you are predicting on small fragments.
        Or, if you have a fasta file of reads, that is correct input too.

        Outputs the predictions in predict_dir as a fasta file with
        corresponding a corresponding taxid file.

    4) ./opal.py eval reference_file predicted_labels [-h]

        Naive evaluation of prediction accuracy.

    5) ./opal.py simulate [--optional-arguments] test_dir train_dir out_dir [-h]

        Runs a full pipeline training on data in train_dir, testing on
        data in test_dir, and outputting everything under out_dir in the
        following directory structure:

        1frag/
            simulated test data (drawn fragments) are saved here.
            (ignored if --do-not-fragment)
        2model/
            classifier will be saved here.
        3predict/
            fragment classifications are saved here.

Contact
    Yunan Luo, luoyunan@gmail.com (original author)
    Yun William Yu, contact@yunwilliamyu.net (author of Python rewrite)

Acknowledgement
    This implementation of Opal is a complete reimplementation based on the ideas and software from the following paper:
    K. Vervier, P. Mahe, M. Tournoud, J.-B. Veyrieras, and J.-P. Vert. Large-scale Machine Learning for Metagenomics Sequence Classification , Technical report HAL-01151453, May, 2015.

