This code is associated with the following manuscript:
Metagenomic binning through low density hashing. Yunan Luo, Yun William Yu, Jianyang Zeng, Bonnie Berger, and Jian Peng. Submitted for publication.

Further information can be found at http://opal.csail.mit.edu/

0. Requirments
    Python 2.7 (this code fails on Python 3)
    Vowpal Wabbit >= 8.3.0 (the command-line program must be installed)
    scikit-learn (sklearn)
    pandas
    scipy

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
    If the dependencies aren't already installed:
        On Ubuntu, run the following:
        sudo apt-get install vowpal-wabbit ("vw" must be on the command line)
        pip install sklearn pandas scipy
    Note that the Python pip package "vowpalwabbit" doesn't actually build
    the command line tool, so it will not work with Opal.

    Then afterwards download some data files for the tests:
        bash SETUP.sh (just downloads some data files to play with)

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

