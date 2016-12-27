#!/bin/bash

ROOT=../..
if hash readlink
then
    ROOT=`readlink -f $ROOT`
fi

export PATH=$ROOT/util/ext/gdl-1.1/GDL/bin:$ROOT/util/ext/gdl-1.1/GDL/include:$PATH
export LD_LIBRARY_PATH=$ROOT/util/ext/gdl-1.1/GDL/lib:$LD_LIBRARY_PATH

# specify path to fasta2skm
fasta2skm=$ROOT/util/fasta2skm

# specify reference database
DB=A1

# specify parameters
NBATCHES=2
K=64

function show_help {
echo This script assumes that it is run from the directory in which it lives.
echo Usage: $0 -d [database] --nbatches [NBATCHES] --kmer [KMER_LENGTH] --outputdir [OUTPUTDIR] --modeldir [MODELDIR] --testdir [TESTDIR]
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
        -k|--kmer)
            if [ -n "$2" ]; then
                K=$2
                shift
            else
                printf 'ERROR: "--kmer" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --kmer=?*)
            K=${1#*=}
            ;;
        --kmer=)
            printf 'ERROR: "--kmer" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        -o|--outputdir)
            if [ -n "$2" ]; then
                outputDir=$2
                shift
            else
                printf 'ERROR: "--outputdir" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --outputdir=?*)
            outputDir=${1#*=}
            ;;
        --outputdir=)
            printf 'ERROR: "--outputdir" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        -m|--modeldir)
            if [ -n "$2" ]; then
                modelDir=$2
                shift
            else
                printf 'ERROR: "--modeldir" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --modeldir=?*)
            modelDir=${1#*=}
            ;;
        --modeldir=)
            printf 'ERROR: "--modeldir" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        -t|--testdir)
            if [ -n "$2" ]; then
                testDir=$2
                shift
            else
                printf 'ERROR: "--testdir" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --testdir=?*)
            testDir=${1#*=}
            ;;
        --testdir=)
            printf 'ERROR: "--testdir" requires a non-empty option argument.\n' >&2
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

if [[ -z $modelDir ]]; then
    modelDir=$ROOT/output/$DB/2-build-models/train_$DB-db
fi
model=$modelDir/vw-model_batch-${NBATCHES}.model
dico=$modelDir/vw-dico.txt

# make predictions fragments 

if [[ -z $testDir ]]; then
    testDir=$ROOT/output/$DB/1-generate-test-datasets
fi

fasta=$testDir/test.fragments.fasta

if [[ -z $outputDir ]]; then
    outputDir=$ROOT/output/$DB/3-make-predictions
fi
mkdir -p $outputDir
prefix=$outputDir/test.fragments.$DB-db

# get vw predictions
$fasta2skm -i $fasta -k $K  -p $modelDir/patterns.txt | vw -t -i $model -p $prefix.preds.vw
# convert vw class to taxid
python vw-class-to-taxid.py ${prefix}.preds.vw $dico ${prefix}.preds.taxid
#python eval.py -d $DB
python eval.py -p ${prefix}.preds.taxid -r $testDir/test.fragments.taxid


