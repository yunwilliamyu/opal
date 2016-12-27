#!/bin/bash

ROOT=../..
if hash readlink
then
    ROOT=`readlink -f $ROOT`
fi

export PATH=$ROOT/util/ext/gdl-1.1/GDL/bin:$ROOT/util/ext/gdl-1.1/GDL/include:$PATH
export LD_LIBRARY_PATH=$ROOT/util/ext/gdl-1.1/GDL/lib:$LD_LIBRARY_PATH

DB="A1"

# specify path to drawfrag and fasta2skm (fasta to spaced kmer)
drawfrag=$ROOT/util/drawfrag
fasta2skm=$ROOT/util/fasta2skm

# specify fragments parameters #
NBATCHES=2
COVERAGE=0.1
L=200

# specify kmer size #
K=64

# spaced kmer weight, number of hash functions
row_weight=4
numHash=2

# specify generic VW options #
LAMBDA1=0
LAMBDA2=0
NPASSES=10
BITS=31

function show_help {
echo This script assumes that it is run from the directory in which it lives.
echo Usage: $0 -d [database] -l [fragment-length] -c [coverage] --nbatches [NBATCHES] --kmer [KMER_LENGTH] --row_weight [row_weight] --numHash [numHash] --npasses [NPASSES] -o [OUTPUT_DIR]
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

# specify input data #
fasta=$ROOT/data/$DB/train/$DB.train.fasta
taxids=$ROOT/data/$DB/train/$DB.train.taxid
# extract number of labels
NLABELS=$(sort -u $taxids | wc -l)



# specify output data #
if [[ -z $outputDir ]]; then
    outputDir=$ROOT/output/$DB/2-build-models
fi
mkdir -p $outputDir
# define output "dictionary" : taxid <--> wv classes
dico=$outputDir/vw-dico.txt
# define model prefix
modelPrefix=$outputDir/vw-model


# specify temporary directory (to shuffle input files) #
TMPDIR=$ROOT/output/$DB/TMP
mkdir -p $TMPDIR

# generate LDPC spaced pattern
python $ROOT/util/ldpc.py -k $K -t $row_weight -m $numHash -d $outputDir/patterns.txt


# loop on batches #
SEED=42
for i in $(seq 1 ${NBATCHES})
do
	# modify random seed
        SEED=$(($SEED+$i))

	#  draw fragments
	fastaBatch=$outputDir/train.batch-$i.fasta
	gi2taxidBatch=$outputDir/train.batch-$i.gi2taxid
	taxidBatch=$outputDir/train.batch-$i.taxid
	$drawfrag -i $fasta -t $taxids -l $L -c $COVERAGE -o $fastaBatch -g $gi2taxidBatch -s $SEED
	# extract taxids
	cut -f 2 $gi2taxidBatch > $taxidBatch

	# learn model
	if [[ $i -eq 1 ]]; then
		# first iteration : no previous model to build upon 		
		$fasta2skm -i $fastaBatch -t $taxidBatch -k $K -d $dico -p $outputDir/patterns.txt | awk 'BEGIN{srand($SEED);} {printf "%06d %s\n", rand()*1000000, $0;}' | sort -n -T $TMPDIR | cut -c8- | vw --random_seed $SEED --cache_file ${modelPrefix}_batch.cache --passes $NPASSES -f ${modelPrefix}_batch-${i}.model.sr --oaa $NLABELS -b $BITS --l1 $LAMBDA1 --l2 $LAMBDA2 --save_resume
	else
		# last iteration : build final model (no save-resume mechanism)
		if [[ $i -eq ${NBATCHES} ]]; then
			$fasta2skm -i $fastaBatch -t $taxidBatch -k $K -d $dico -p $outputDir/patterns.txt | awk 'BEGIN{srand($SEED);} {printf "%06d %s\n", rand()*1000000, $0;}' | sort -n -T $TMPDIR | cut -c8- | vw --random_seed $SEED --cache_file ${modelPrefix}_batch.cache -f ${modelPrefix}_batch-${i}.model -i ${modelPrefix}_batch-$(($i - 1)).model.sr
		else
			# futher iterations : save-resume mechanism
			$fasta2skm -i $fastaBatch -t $taxidBatch -k $K -d $dico -p $outputDir/patterns.txt | awk 'BEGIN{srand($SEED);} {printf "%06d %s\n", rand()*1000000, $0;}' | sort -n -T $TMPDIR | cut -c8- | vw --random_seed $SEED --cache_file ${modelPrefix}_batch.cache -f ${modelPrefix}_batch-${i}.model.sr -i ${modelPrefix}_batch-$(($i - 1)).model.sr --save_resume
		fi
			

		# remove old model
		rm ${modelPrefix}_batch-$(($i - 1)).model.sr
	fi

	# remove temp files
	rm  $outputDir/train.batch-$i*
	rm ${modelPrefix}_batch.cache
done

# Fix to rename the first model for future scripts
if [[ $NBATCHES -eq 1 ]]; then
    mv ${modelPrefix}_batch-1.model.sr ${modelPrefix}_batch-1.model
fi

