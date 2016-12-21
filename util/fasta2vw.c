/*
Copyright (c) by respective owners including bioMerieux, and
individual contributors Kevin Vervier and Pierre Mahe. All rights reserved.  Released under a BSD (revised)
license as described in the file LICENSE.
 */
#include <getopt.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>
#include <math.h>

#include <gdl/gdl_common.h>
#include <gdl/gdl_version.h>
#include <gdl/gdl_errno.h>
#include <gdl/gdl_sort_uint.h>
#include <gdl/gdl_io.h>
#include <gdl/gdl_hash.h>
#include <gdl/gdl_list.h>
#include <gdl/gdl_string.h>
#include <gdl/gdl_util.h>

#include <zlib.h>
#include <kseq.h>
KSEQ_INIT(gzFile, gzread);


static gdl_string * PROGRAM = "fasta2vw";

static int help_flag    = 0;
static int verbose_flag = 0;
static int stdout_flag = 1;
static int taxid_flag = 0;

static gdl_string * INPUT = NULL;
static gdl_string * OUTPUT = NULL;
static gdl_string * DICO = NULL;
static gdl_string * K = NULL;
static gdl_string * TAXID = NULL;

static struct option long_options[] =
{
		/* These options set a flag. */
		{"help", no_argument,			&help_flag, 1},
		{"verbose", no_argument,		&verbose_flag, 1},
		/* These options don't set a flag.
         	We distinguish them by their indices. */
		{"input",   required_argument, 0, 'i'},
		{"taxid",   required_argument, 0, 't'},
		{"output",   required_argument, 0, 'o'},
		{"dico", required_argument, 0, 'd'},
		{"kmer_size", required_argument, 0, 'k'},
		{0, 0, 0, 0}
};

static int
parse_argument (int argc, char *argv[])
{
	int c;

	while (1)
	{

		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "i:o:d:k:t:",
				long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
		case 0:
			/* If this option set a flag, do nothing else now. */
			if (long_options[option_index].flag != 0)
				break;
			printf ("option %s", long_options[option_index].name);
			if (optarg)
				printf (" with arg %s", optarg);
			printf ("\n");
			break;
		case 'i':
			INPUT = gdl_string_clone (optarg);
			break;
		case 'o':
			OUTPUT = gdl_string_clone (optarg);
			stdout_flag = 0;
			break;
		case 'd':
			DICO = gdl_string_clone (optarg);
			break;
		case 'k':
			K = gdl_string_clone (optarg);
			break;
		case 't':
			TAXID = gdl_string_clone (optarg);
			taxid_flag = 1;
			break;
		case '?':
			GDL_ERROR_VAL ("Unknown arguments", GDL_EINVAL, -1);
		default:
			GDL_ERROR_VAL ("Bad arguments", GDL_EINVAL, -1);
		}
	}
	return GDL_SUCCESS;
}

static int
check_argument (void)
{
	if (INPUT == 0)
	{
		GDL_ERROR_VAL ("No input file provided",  GDL_FAILURE, 1);
	}
	if ( (TAXID != 0) & (DICO == 0) )
	{
		GDL_ERROR_VAL ("Must provide a dico if option taxid is specified",  GDL_FAILURE, 1);
	}

	if ( (TAXID == NULL) & (DICO != 0) )
	{
		printf ("WARNING : the --dico option useless if taxid not specified \n");
	}
	if(K == 0){
		GDL_ERROR_VAL ("No kmer size provided", GDL_FAILURE, 1);
	}
	return GDL_SUCCESS;
}

static int
help (void)
{
	printf("%s - version 1.0\n", PROGRAM);
   	printf("Copyright (C) 2015 bioMerieux, France\n");
	printf ("\n");
	printf("Extract k-mer profiles of sequences stored in a fasta or fastq file into a Vowpal Wabbit (VW) compliant output file\n");
	printf ("\n");
	printf ("--help\tDisplay a brief help on program usage\n");
	printf ("--verbose\tOutput message on standard output to see what the program is doing\n");
	printf ("\n");
	printf ("--input or -i\t file containing input sequences (fasta or fastq format)\n");
	printf ("--kmer or -k\t k-mer size\n");
	printf ("[--taxid or -t\tone-column file containing the taxid of each input sequence]\n");
	printf ("\t(if not specified, instances are labeled as 1)\n");
	printf ("--[dico or -d\t two-colum file giving the correspondance between taxids and VW class labels]\n");
	printf ("\t(created if does not exist ; possibly updated and overwritten if exists)\n");
	printf ("\t(must be specified if --taxid option is used)\n");
	printf ("\t(useless if --taxid option is not used)\n");
	printf ("[--output or -o\t output file]\n");
	printf ("\t (if not specified, results printed on standard output)\n");
	return GDL_SUCCESS;
}

gdl_string * reverse_transcribe(gdl_string * s){
        // initialize result
        gdl_string * res;
        res = gdl_string_alloc(strlen(s));
        // process sequence
        size_t i;
        for (i = 0; i < strlen(s); i++){
                // reverse
                res[i] = s[strlen(s)-i-1];
                // transcribe
                switch(res[i]){
                        case 'A' : res[i] = 'T'; break;
                        case 'T' : res[i] = 'A'; break;
                        case 'G' : res[i] = 'C'; break;
                        case 'C' : res[i] = 'G'; break;
                        default : break;
                }
        }
        // return
        return(res);
}



int
main (int argc, char *argv[])
{

	int status;

	parse_argument (argc, argv);

	if (help_flag)
	{
		exit (help());
	}

	status = check_argument ();

	if (status == GDL_SUCCESS)
	{

		gzFile f, dico_table;
		kseq_t *seq;
		FILE * stream = 0;
		size_t i, j, l, n, kmer_size, n_kmer_sizes, cptr, nlabels, tx_write, label_write;
		size_t * label_tx, *kmer_sizes;
		gdl_hashtable *taxid_to_label = 0;
		char taxid_name[20];
		FILE * taxid_stream  = 0;
		gdl_string *taxid = NULL;

		//-------------------------//
 		// parse sizes to consider //
		//-------------------------//
		gdl_string ** toks;
		size_t ntok, n1, n2;
		toks = gdl_string_split(K, ",", &ntok); // search for list of values
		if(ntok == 1){
			toks = gdl_string_split(K, "-", &ntok); // search for range of values
			if(ntok == 1){
				n_kmer_sizes = 1;
				kmer_sizes = GDL_CALLOC(size_t, n_kmer_sizes);
				kmer_sizes[0] = atoi(K);
			}else{
				if(ntok != 2){
					GDL_ERROR_VAL ("problem reading kmer sizes",GDL_FAILURE,1);
				}
				n1 = atoi(toks[0]);
				n2 = atoi(toks[1]);
				if(n1 >= n2){
					GDL_ERROR_VAL ("problem reading kmer sizes",GDL_FAILURE,1);
				}
				n_kmer_sizes = n2-n1+1;
				kmer_sizes = GDL_CALLOC(size_t, n_kmer_sizes);
				for(i = 0; i < n_kmer_sizes; i++){
					kmer_sizes[i] = n1 + i;
				}
			}
		}else{
			n_kmer_sizes = ntok;
			kmer_sizes = GDL_CALLOC(size_t, n_kmer_sizes);
			for(i = 0; i < n_kmer_sizes; i++){
				kmer_sizes[i] = atoi(toks[i]);
			}
		}
		if(verbose_flag){
			printf("** considering %lu sizes of kmers : \n", n_kmer_sizes);
			for(i = 0; i < n_kmer_sizes; i++){
				printf("\t- size no %lu = %lu\n", i, kmer_sizes[i]);
			}
		}

		//---------------------//
		// initialize I/O data //
		//---------------------//
		// open input file
		f = gzopen(INPUT, "r");
		seq = kseq_init(f);
		// open taxid file
		if(taxid_flag){
		 	taxid_stream  = gdl_fileopen (TAXID, "r");
			// initialize dictionary hash table (taxid to VW class index)
  			taxid_to_label = gdl_hashtable_alloc (gdl_interface_uint, 0);
			// read dictionary (if exists)
			dico_table = fopen(DICO, "r");
			if(dico_table)
			{
				while(!feof(dico_table))
				{
					// extract (taxid / class) pair
					fscanf(dico_table,"%lu\t%lu",&tx_write,&label_write);
					sprintf(taxid_name,"%lu",tx_write);
					//sprintf(label_name,"%lu",label_write);
					// add to hash hash table (sanity check : check it does not exist already)
					if((label_tx = (size_t*) gdl_hashtable_lookup (taxid_to_label, taxid_name)) == 0)
					{
						label_tx = GDL_CALLOC (size_t, 1);
						*label_tx = (size_t) (label_write);
						gdl_hashtable_add(taxid_to_label, taxid_name, label_tx, 1);
					}
				}
				fclose(dico_table);
			}
		}
		// open output file
		if(!stdout_flag)
		{
			stream = gdl_fileopen (OUTPUT, "w");
		}
		//-----------------------//
  		// process each sequence //
  		//----------------------//
		cptr = 0;
		nlabels = 1;
		while (kseq_read(seq) >= 0)
		{
			// log message //
			// ------------//
			cptr++;
			if( verbose_flag & (cptr>0) & (cptr%1000 == 0)){
				printf("\t- processing sequence no %lu (length = %lu)\n", cptr, seq->seq.l);
			}	
			// extract taxid, convert to class and write to output file //
			//----------------------------------------------------------//
			if(taxid_flag){
				// read from input file
				gdl_getline (&taxid, &n, taxid_stream);
				// check exists in dico and insert if not
				if((label_tx = (size_t*) gdl_hashtable_lookup (taxid_to_label, taxid)) == 0 ){
					label_tx = GDL_CALLOC (size_t, 1);
					//new key
					*label_tx = (size_t) (nlabels);
					nlabels++;
					// add 
					gdl_hashtable_add(taxid_to_label, taxid, label_tx, 1);
				}
				// write to output file
				if(stdout_flag){
					printf("%lu |", *label_tx);
				}else{
                	  		fprintf(stream,"%lu |", *label_tx);
              			}
			}else{
				// write 1 to output file
				if(stdout_flag){
					printf("1 |");
				}else{
                	  		fprintf(stream,"1 |");
              			}
			}
			// output kmers //
			//--------------//
			// loop on kmer sizes
			for(l = 0; l < n_kmer_sizes; l++){
				kmer_size = kmer_sizes[l];
				//test if the sequence length is greater than K
				if(seq->seq.l > kmer_size){
					for(i = 0; i <= seq->seq.l - kmer_size; i++){
						// extract forward kmer
						gdl_string *kmer;
						kmer = 0;
						kmer = gdl_string_alloc (kmer_size);
						strncpy (kmer, seq->seq.s+i, kmer_size);
						// convert to upper case
						for(j = 0; j < strlen(kmer); j++){
							kmer[j] = toupper(kmer[j]);
						}
   						// print foward kmer
						if(!stdout_flag){
							fprintf(stream," %s",kmer);
						}else{
							printf(" %s",kmer);
						}
						// construct reverse kmer
						gdl_string *rt_kmer;
						rt_kmer = 0;
						rt_kmer = reverse_transcribe(kmer);
						// print reverse kmers
						if(!stdout_flag){
							fprintf(stream," %s",rt_kmer);
						}else{
							printf(" %s",rt_kmer);
						}
						// free kmers
						gdl_string_free(kmer);
						gdl_string_free(rt_kmer);
					}
				}else{
					// simply out "forward" sequence, no upper-case transformation
					if(!stdout_flag){
						fprintf(stream," %s",seq->seq.s);
					}else{
						printf(" %s",seq->seq.s);
					}
				}
			}
			// add an empty line between sequences //
			// -----------------------------------//
			if(!stdout_flag){
				fprintf(stream,"\n");
			}else{
				printf("\n");
			}
			// free taxid string //
			gdl_string_free(taxid);
			taxid = 0;
		} 
	
		//-------------------//
		// process I/O files //
		//-------------------//
		// close input file
		kseq_destroy(seq);
		gzclose(f);
		// close output file
		if(!stdout_flag){
			gdl_fileclose (OUTPUT, stream);
		}
		// write hash-table if necessary
		if(taxid_flag){	
			stream = gdl_fileopen(DICO, "w");
			gdl_hashtable_itr * hash_itr = gdl_hashtable_iterator(taxid_to_label);
			int * val = 0;
			do {
				val = (int *) gdl_hashtable_iterator_value(hash_itr);
				fprintf(stream, "%s\t%d\n", gdl_hashtable_iterator_key(hash_itr),*val);
			} while (gdl_hashtable_iterator_next(hash_itr));
			gdl_hashtable_iterator_free(hash_itr);
			gdl_fileclose (DICO, stream);
			// free hash table
			gdl_hashtable_free(taxid_to_label);
		}
		gdl_string_free(DICO);
		gdl_string_free(OUTPUT);
		gdl_string_free(INPUT);
	}
	exit (0);
}
