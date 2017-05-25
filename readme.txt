This code is associated with the following manuscript:
Metagenomic binning through low density hashing. Yunan Luo, Y. William Yu, Jianyang Zeng, Bonnie Berger, and Jian Peng. Submitted for publication.

Upon publication, further information can be found at http://opal.csail.mit.edu/

0. Requirments
    Vowpal Wabbit >= 8.3.0
    scikit-learn

    This code has been tested with GCC 4.8.4 on Ubuntu 14.04 and 16.04, running
    under Bash 4.3.11. There are reports of compilation errors on Mac OS X, so
    YMMV.

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
    ext/: external libararies.
    test/: test drawfrag.c and fasta2skm.c
    drawfrag.c: draw fragments from fasta records.
    fasta2skm.c: construct feature (spaced k-mer profile), and convert to VW input format.
    ldpc.py: generate LSH function using LDPC code.
    kseq.h: parse FASTA files

2. Install and test:
    bash SETUP.sh

   Old manual installation:
    1) install
        $ cd util
        $ sh INSTALL.sh
    2) test
        $ cd util/test
        $ sh test_skm.sh
        Results will be saved in output/

3. Usage: opal.py assumes it lives in the current directory structure, but can be symlinked elsewhere.

Modes:
    (default --optional-arguments such as k-mer length, fragment size,
    hash functions, etc. are set for quick run, and so will give terrible
    predictions. If you're actually going to use for training, be sure to
    set better parameters.)

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
    This implementation of Opal is adapted from the source code of the following paper:
    K. Vervier, P. Mahe, M. Tournoud, J.-B. Veyrieras, and J.-P. Vert. Large-scale Machine Learning for Metagenomics Sequence Classification , Technical report HAL-01151453, May, 2015.

