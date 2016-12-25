export PATH=../../../util/ext/gdl-1.1/GDL/bin:../../../util/ext/gdl-1.1/GDL/include:$PATH
export LD_LIBRARY_PATH=../../../util/ext/gdl-1.1/GDL/lib:$LD_LIBRARY_PATH
#export PATH=/mnt/work/jpeng/Yunan/Metagenomics/code/vowpal_wabbit/vowpalwabbit:$PATH
#export PATH=/mnt/work/jpeng/Yunan/Metagenomics/code/boost_1_60_0/BOOST/include:$PATH
#export LD_LIBRARY_PATH=/mnt/work/jpeng/Yunan/Metagenomics/code/boost_1_60_0/BOOST/lib:$LD_LIBRARY_PATH

# specify path to fasta2skm
fasta2skm=../../../util/fasta2skm

# specify reference database
DB=A1

# specify parameters
NBATCHES=2
K=64

function show_help {
echo This script assumes that it is run from the directory in which it lives.
echo Usage: $0 -d [database] --nbatches [NBATCHES] --kmer [KMER_LENGTH] 
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

modelDir=../../2-build-models/output/train_$DB-db
model=$modelDir/vw-model_batch-${NBATCHES}.model
dico=$modelDir/vw-dico.txt

# make predictions fragments 
# generate test fragments from fasta
fasta=../../1-generate-test-datasets/output_${DB}/test.fragments.fasta
# use existing test fragments
#fasta=../../../data/$DB/test/${DB}.test.fasta
prefix=../output/test.fragments.$DB-db

# get vw predictions
$fasta2skm -i $fasta -k $K  -p $modelDir/patterns.txt | vw -t -i $model -p $prefix.preds.vw
# convert vw class to taxid
python vw-class-to-taxid.py $prefix.preds.vw $dico $prefix.preds.taxid
python eval.py -d $DB


