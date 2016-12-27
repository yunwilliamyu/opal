0. Requirments
    Vowpal Wabbit >= 8.3.0
    scikit-learn

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
output/: All scripts will output to subfolders here corresponding to data/[example]
    [example]/
        1-generate-test-datasets/
            simulated test data (drawn fragents) are saved here.
        2-build-models/
            classifier will be saved here.
        3-make-predictions/
            fragment classifications are saved here.
        TMP/:
            temp files.
src/
    1-generate-test-datasets/
        opal-generate.sh: draw fragments for testing.
    2-build-models/
        opal-train.sh: train a multi-class svm using VW library.
    3-make-predictions/
        opal-predict.sh: classify fragments.
        vw-class-to-taxid.py: map VW labels to tax id.
        eval.py: evaluation of supervised classification
util/
    ext/: external libararies.
    test/: test drawfrag.c and fasta2skm.c
    drawfrag.c: draw fragments from fasta records.
    fasta2skm.c: construct feature (spaced k-mer profile), and convert to VW input format.
    ldpc.py: generate LSH function using LDPC code.
    kseq.h: parse FASTA files

2. Install and test
    1) install
        $ cd util
        $ sh INSTALL.sh
    2) test
        $ cd util/test
        $ sh test_skm.sh
        Results will be saved in output/

3. Usage: runall.sh will automatically do all of the following steps and can be
        modified with command line parameters.

    1) simulate fragments for testing.
        $ cd src/1-generate-test-datasets
        $ sh opal-generate.sh

        Simulated fragments will be saved in output/1-generate-test-datasets.
        If existing test fragments are provided, this step can be ignored, 
        but directory of the existing test fragments should be set manually 
        in step 3). (see below)

    2) train classifier
        $ cd src/2-build-models
        $ bash opal-train.sh

        Important parameters:
        DB: name of the folder of training/test data.
        NBATCHES: number of train batches.
        L: fragment length.
        COVERAGE: average coverage of each position in the sequence.
        K: k-mer size.
        row_weight: how many positions will be randomly chosen in the contiguous k-mer. (k should be multiple of row_weight, e.g., k=32, row_weight=16)
        numHash: number of hashing function
        (The current parameters are set as such for quick example, so the performance is not good.)

    3) make prediction
        $ cd src/3-make-predictions
        $ sh opal-predict.sh
        if existing test fragments are provided, they should be located in the standard location under data/example/test/ and opal-predict.sh should be called with option --use-existing

Contact
    Yunan Luo, luoyunan@gmail.com

Acknowledgement
    This implementation of Opal is adapted from the source code of the following paper:
    K. Vervier, P. Mahe, M. Tournoud, J.-B. Veyrieras, and J.-P. Vert. Large-scale Machine Learning for Metagenomics Sequence Classification , Technical report HAL-01151453, May, 2015.

