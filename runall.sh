#!/bin/bash
# This runs the full pipeline of simulating fragments, training, and then prediction.

DB=subspecies
L=200
COVERAGE=1
K=128
NUMHASH=32
NPASSES=1000
ROWWEIGHT=32
NBATCHES=64

echo ================
echo Generating Test Datasets
date
echo "bash opal-generate.sh -d $DB -l $L -c $COVERAGE 2>&1"
echo ================
cd src/1-generate-test-datasets/src
bash opal-generate.sh -d $DB -l $L -c $COVERAGE 2>&1
cd ../../..

echo ================
echo Building Models
date
echo "bash opal-train.sh -d $DB -l $L -c $COVERAGE --nbatches $NBATCHES --kmer $K --row_weight $ROWWEIGHT --numHash $NUMHASH --npasses $NPASSES 2>&1"
echo ================
cd src/2-build-models/src
rm .cache
bash opal-train.sh -d $DB -l $L -c $COVERAGE --nbatches $NBATCHES --kmer $K --row_weight $ROWWEIGHT --numHash $NUMHASH --npasses $NPASSES 2>&1
cd ../../..

echo ================
echo Making Predictions
date
echo "bash opal-predict.sh -d $DB --nbatches $NBATCHES --kmer $K 2>&1"
echo ================
cd src/3-make-predictions/src
bash opal-predict.sh -d $DB --nbatches $NBATCHES --kmer $K 2>&1
cd ../../..


