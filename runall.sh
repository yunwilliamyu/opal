#!/bin/bash
# This runs the full pipeline of simulating fragments, training, and then prediction.

DB=A1
L=32
COVERAGE=0.1
K=4
NUMHASH=2
NPASSES=1
ROWWEIGHT=1
NBATCHES=1

#DB=subspecies
#L=200
#COVERAGE=1
#K=128
#NUMHASH=32
#NPASSES=1000
#ROWWEIGHT=32
#NBATCHES=64

function show_help {
echo This script assumes that it is run from the directory in which it lives.
echo Usage: $0 -d [database] -l [fragment-length] -c [coverage] --nbatches [NBATCHES] --kmer [KMER_LENGTH] --row_weight [row_weight] --numHash [numHash] --npasses [NPASSES]

}

while :; do
    case $1 in
        -h|-\?|--help)
            show_help
            exit 0
            ;;
        -d|--database)
            if [ -n "$2" ]; then
                DB=$2
                shift
            else
                printf 'ERROR: "--database" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --database=?*)
            DB=${1#*=}
            ;;
        --database=)
            printf 'ERROR: "--database" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        --nbatches)
            if [ -n "$2" ]; then
                NBATCHES=$2
                shift
            else
                printf 'ERROR: "--nbatches" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --nbatches=?*)
            NBATCHES=${1#*=}
            ;;
        --nbatches=)
            printf 'ERROR: "--nbatches" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        -l|--length)
            if [ -n "$2" ]; then
                L=$2
                shift
            else
                printf 'ERROR: "--length" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --length=?*)
            L=${1#*=}
            ;;
        --length=)
            printf 'ERROR: "--length" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        -k|--kmer)
            if [ -n "$2" ]; then
                K=$2
                shift
            else
                printf 'ERROR: "--kmer" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --npasses)
            if [ -n "$2" ]; then
                NPASSES=$2
                shift
            else
                printf 'ERROR: "--npasses" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --npasses=?*)
            NPASSES=${1#*=}
            ;;
        --npasses=)
            printf 'ERROR: "--npasses" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        -c|--coverage)
            if [ -n "$2" ]; then
                COVERAGE=$2
                shift
            else
                printf 'ERROR: "--coverage" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --coverage=?*)
            COVERAGE=${1#*=}
            ;;
        --coverage=)
            printf 'ERROR: "--coverage" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        --kmer=?*)
            K=${1#*=}
            ;;
        --kmer=)
            printf 'ERROR: "--kmer" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        --)
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        --row_weight)
            if [ -n "$2" ]; then
                row_weight=$2
                shift
            else
                printf 'ERROR: "--row_weight" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --row_weight=?*)
            row_weight=${1#*=}
            ;;
        --row_weight=)
            printf 'ERROR: "--row_weight" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        --numHash)
            if [ -n "$2" ]; then
                numHash=$2
                shift
            else
                printf 'ERROR: "--numHash" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --numHash=?*)
            numHash=${1#*=}
            ;;
        --numHash=)
            printf 'ERROR: "--numHash" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        --)
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)
            break
    esac
    shift
done

exec > >(tee -i output/${DB}/`date +%s`.log)

echo ================
echo Invocation
echo $0 $@
echo ================

echo ================
echo Generating Test Datasets
date
echo "bash opal-generate.sh -d $DB -l $L -c $COVERAGE 2>&1"
echo ================
cd src/1-generate-test-datasets
bash opal-generate.sh -d $DB -l $L -c $COVERAGE 2>&1
cd ../..

echo ================
echo Building Models
date
echo "bash opal-train.sh -d $DB -l $L -c $COVERAGE --nbatches $NBATCHES --kmer $K --row_weight $ROWWEIGHT --numHash $NUMHASH --npasses $NPASSES 2>&1"
echo ================
cd src/2-build-models
bash opal-train.sh -d $DB -l $L -c $COVERAGE --nbatches $NBATCHES --kmer $K --row_weight $ROWWEIGHT --numHash $NUMHASH --npasses $NPASSES 2>&1
cd ../..

echo ================
echo Making Predictions
date
echo "bash opal-predict.sh -d $DB --nbatches $NBATCHES --kmer $K 2>&1"
echo ================
cd src/3-make-predictions
bash opal-predict.sh -d $DB --nbatches $NBATCHES --kmer $K 2>&1
cd ../..

echo ================
echo Fine.
date
echo ================
