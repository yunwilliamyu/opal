#!/bin/bash
# This runs the full pipeline of simulating fragments, training, and then prediction.

DB=A1
L=16
COVERAGE=1
K=8
NUMHASH=1
NPASSES=1
ROWWEIGHT=1
NBATCHES=1
fragment_test_set=true

#DB=subspecies
#L=200
#COVERAGE=1
#K=128
#NUMHASH=32
#NPASSES=1000
#ROWWEIGHT=32
#NBATCHES=64

ROOT=.
if hash readlink
then
    ROOT=`readlink -f $ROOT`
fi

function show_help {
echo This script assumes that it is run from the directory in which it lives.
echo Usage: $0 -d [database] -l [fragment-length] -c [coverage] --nbatches [NBATCHES] --kmer [KMER_LENGTH] --row_weight [row_weight] --numHash [numHash] --npasses [NPASSES] [--fragment-test-set]

}

while :; do
    case $1 in
        -h|-\?|--help)
            show_help
            exit 0
            ;;
        --fragment-test-set)
            fragment_test_set=true
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

RUNID=`date +%s`
OUTDIR=$ROOT/output/${DB}/${RUNID}
DATADIR=$ROOT/data/${DB}
mkdir -p $OUTDIR
exec > >(tee -i ${OUTDIR}/run.log) 2>&1
set -o posix > ${OUTDIR}/params.txt

echo ================
echo Invocation
echo $0 $@
echo ================

echo ================
echo Generating Test Datasets
date
echo ================
set -x
cd src/1-generate-test-datasets
bash opal-generate.sh -i $DATADIR -l $L -c $COVERAGE -o $OUTDIR/1-generate-test-datasets 2>&1
cd ../..
set +x

echo ================
echo Building Models
date
echo ================
set -x
cd src/2-build-models
bash opal-train.sh -i $DATADIR -l $L -c $COVERAGE --nbatches $NBATCHES --kmer $K --row_weight $ROWWEIGHT --numHash $NUMHASH --npasses $NPASSES -o $OUTDIR/2-build-models 2>&1
cd ../..
set +x

echo ================
echo Making Predictions
date
echo ================
if [ "$fragment_test_set" = true ]; then
    test_dir=$OUTDIR/1-generate-test-datasets
else
    test_dir=$ROOT/data/$DB/test
fi
set -x
cd src/3-make-predictions
bash opal-predict.sh --nbatches $NBATCHES --kmer $K -m $OUTDIR/2-build-models -t $test_dir -o $OUTDIR/3-make-predictions 2>&1
cd ../..
set +x

echo ================
echo Fine.
date
echo ================
