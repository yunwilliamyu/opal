This code is associated with the following manuscript:
Yunan Luo, Yun William Yu, Jianyang Zeng, Bonnie Berger, Jian Peng; Metagenomic binning through low-density hashing, Bioinformatics, Volume 35, Issue 2, 15 January 2019, Pages 219–226, https://doi.org/10.1093/bioinformatics/bty611

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
    First download and expand Opal to a working directory:
        git clone https://github.com/yunwilliamyu/opal.git
        cd opal
        bash SETUP.sh

    The provided setup script will test your environment for dependencies and
    download some example data files to play with. It will fail when it does
    not find a required dependency, and attempt to give instructions for how to
    install them. On Ubuntu, the following should install all needed
    dependencies:
        sudo apt-get install vowpal-wabbit python-pip
        pip install pandas, sklearn
    Note that you must install the command-line version of vowpal-wabbit, as
    Opal depends on the "vw" command.

    To test Opal after setup if you downloaded the example data files using SETUP.sh:
        ./opal.py simulate data/A1/test data/A1/train out_dir -r
    Note that this may take a very long time to train, so you may change the coverage
    to have a quicker (but horrendously inaccurate) test of the pipeline as follows:
        ./opal.py simulate data/A1/test data/A1/train out_dir -r -c0.1

    For context, the coverage parameter -c is the most important one,
    controlling how much training data the model is given. A coverage parameter
    of 1 implies that on average every position in the original fasta file is
    covered by one training example. Similarly, coverage of 10 implies every
    position is covered by 10 training examples, and coverage of 0.1 implies
    that you expect only 10% of positions to be covered by training examples.
    All else held equal, increasing the coverage parameter will increase the
    accuracy, though at the cost of longer training times. For most cases,
    setting a coverage of 15 (-c15) gives reasonable accuracy, though if you
    wish to do subspecies classification, you may need a much higher coverage
    level, depending on how closely related your strains are. (See paper for
    more details). For *very* high coverage (e.g. 1000), we recommend running
    in multiple batches. Instead of "-c1000", use "-c1 --num-batches=1000".

    As an aside, if you are classifying DNA, you will probably also want the
    "-r" argument on everything to also treat reverse complements as the same.
    We do not use it by default for backwards compatability since some users
    use Opal on protein sequences, for which reverse complements are not
    defined or meaningful.

    The simulate command is equivalent to the following mult-step process:
    Create output directory:
        mkdir out_dir

    Fragmentation to prepare the fragmented testing data (not used for training):
        ./opal.py frag data/A1/test out_dir/1frag [-c0.1]
    
    Training based on the training directory fasta:
        ./opal.py train data/A1/train out_dir/2model [-c0.1] [-r]

    Predicting examples from the fasta file found in the test directory:
        ./opal.py predict out_dir/2model out_dir/1frag out_dir/3predict [-r]

    Evaulating the accuracy of the predicted assignments:
        ./opal.py eval out_dir/1frag/test.fragments.taxid out_dir/3predict/test.fragments-db.preds.taxid

    Note that we fragment the testing data because for the simulation, we
    accept fasta files consisting of entire genomes, so we must fragment it in
    order to get test reads. If you are classifying the Fasta files from an NGS
    pipeline, or benchmarking against a real dataset, you will probably want to
    follow the following pipeline. To duplicate the experiment we ran for Opal
    in Figure 3 of our paper on the A1 dataset, you'll want to do the following:

    We will start with:
        data/A1/train/A1.train.fasta: each Fasta entry in the file is a genome
        data/A1/train/A1.train.taxid: each line corresponds to the taxid of the fasta entries

        data/A1/test/A1.test.fasta: each Fasta entry is a read to be classified
        data/A1/test/A1.test.taxid: each line is the ground-truth taxid. Obviously,
            we only have the ground truth because this is a benchmark.

    Create output directory:
        mkdir out_dir

    Train on the given genomes, drawing fragments of length 200, and setting coverage=15
    (by combining -c1 and --num-batches=15).
        ./opal.py train data/A1/train out_dir/2model -c1 --num-batches=15 -r --frag-length=200

    Predict on the test data:
        ./opal.py predict out_dir/2model data/A1/test out_dir/3predict -r

    Evaluate the predicted assignments:
        ./opal.py eval data/A1/test/A1.test.taxid out_dir/3predict/test.fragments-db.preds.taxid

    For other evaluations (e.g. F1-score), the predicted taxids are just in the file
        out_dir/3predict/test.fragments-db.preds.taxid
    So they can be compared using whatever accuracy measure we desire.

    Note that the results may differ slightly from the paper itself because of
    both updates to the code since running our original experiments and the
    randomness involved in the training process.

    Also, the training process can be extremely slow, and can take many hours.
    On our test machine, the above pipeline takes around 4 hours to complete.
    This is of course why we provide the much faster commands above for playing
    with the code. However, once trained on a set of genomes, classification is
    relatively fast, and in line with the speed of other metagenomic binners.

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

