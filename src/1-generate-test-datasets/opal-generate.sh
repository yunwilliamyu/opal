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
echo Usage: $0 -i "[datadir]" -l "[fragment-length]" -c "[coverage]" -o "[outputdir]"
}

while getopts ":hi:l:c:o:" opt; do
    case $opt in
        h)
            show_help
            exit 0
            ;;
        i)
            DATADIR=$OPTARG
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

echo Data_dir=$DATADIR
echo length=$L
echo coverage=$COVERAGE

# specify input data

testfasta_files="$DATADIR/test/*.fasta"
testfasta_array=( $testfasta_files )
fasta="${testfasta_array[0]}"
testid_files="$DATADIR/test/*.taxid"
testid_array=( $testid_files )
taxids="${testid_array[0]}"

# set seed (for reproducibility)
SEED=42
mkdir -p $OUTPUTDIR

# draw fragments
$drawfrag -i $fasta -t $taxids -l $L -c $COVERAGE -o $OUTPUTDIR/test.fragments.fasta -g $OUTPUTDIR/test.fragments.gi2taxid -s $SEED

# extract taxids
cut -f 2 $OUTPUTDIR/test.fragments.gi2taxid > $OUTPUTDIR/test.fragments.taxid
