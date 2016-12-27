#!/bin/bash

ROOT=../..
if hash readlink
then
    ROOT=`readlink -f $ROOT`
fi

export PATH=$ROOT/util/ext/gdl-1.1/GDL/bin:$ROOT/util/ext/gdl-1.1/GDL/include:$PATH
export LD_LIBRARY_PATH=$ROOT/util/ext/gdl-1.1/GDL/lib:$LD_LIBRARY_PATH

OPTIND=1 # reset for getopts

# specify path to drawfrag
drawfrag=$ROOT/util/drawfrag
DB="A1"

# specify fragments parameters (length and cover ~ number)
L=64
COVERAGE=1

function show_help {
echo This script assumes that it is run from the directory in which it lives.
echo Usage: $0 -d "[database]" -l "[fragment-length]" -c "[coverage]" -o "[outputdir]"
}

while getopts ":hd:l:c:o:" opt; do
    case $opt in
        h)
            show_help
            exit 0
            ;;
        d)
            DB=$OPTARG
            ;;
        l)
            L=$OPTARG
            ;;
        c)
            COVERAGE=$OPTARG
            ;;
        o)
            OUTPUTDIR=$OPTARG
            ;;
        \?)
            echo "Invalid option -$OPTARG"
            ;;
    esac
done
shift "$((OPTIND-1))"

echo Database=$DB
echo length=$L
echo coverage=$COVERAGE

# specify input data
fasta=$ROOT/data/$DB/test/$DB.test.fasta
taxids=$ROOT/data/$DB/test/$DB.test.taxid

# set seed (for reproducibility)
SEED=42
if [[ -z $OUTPUTDIR ]]; then
    OUTPUTDIR=$ROOT/output/$DB/1-generate-test-datasets
fi
mkdir -p $OUTPUTDIR

# draw fragments
$drawfrag -i $fasta -t $taxids -l $L -c $COVERAGE -o $OUTPUTDIR/test.fragments.fasta -g $OUTPUTDIR/test.fragments.gi2taxid -s $SEED

# extract taxids
cut -f 2 $OUTPUTDIR/test.fragments.gi2taxid > $OUTPUTDIR/test.fragments.taxid
