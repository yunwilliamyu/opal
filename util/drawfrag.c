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


static gdl_string * PROGRAM = "drawfrag";

static int help_flag    = 0;
static int verbose_flag = 0;
static int atgc_flag = 0;

static gdl_string * INPUT = NULL;
static gdl_string * OUTPUT = NULL;
static gdl_string * TAXIDS = NULL;
static gdl_string * GI2TAXID = NULL;
static double COVERAGE = 0;
static int SIZE = 0;
static int SEED = 0;

static struct option long_options[] =
{
		/* These options set a flag. */
		{"help", no_argument,          &help_flag, 1},
		{"verbose", no_argument,       &verbose_flag, 1},
		{"atgc", no_argument,       &atgc_flag, 1},
		/* These options don't set a flag.
         	We distinguish them by their indices. */
		{"input",   required_argument, 0, 'i'},
		{"output",   required_argument, 0, 'o'},
		{"coverage",   required_argument, 0, 'c'},
		{"taxids", required_argument, 0, 't'},
		{"size",   required_argument, 0, 'k'},
		{"gi2taxid", required_argument,0, 'g'},
		{"seed", required_argument,0, 's'},
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

		c = getopt_long (argc, argv, "g:i:l:c:o:s:t:",
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
		case 'g':
			GI2TAXID = gdl_string_clone (optarg);
			break;
		case 'i':
			INPUT = gdl_string_clone (optarg);
			break;
		case 'l':
			SIZE = atoi(optarg);
			break;
		case 'c':
			COVERAGE = atof(optarg);
			break;
		case 'o':
			OUTPUT = gdl_string_clone (optarg);
			break;
		case 's':
			SEED = atoi (optarg);
			break;
		case 't':
			TAXIDS = gdl_string_clone (optarg);
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
		GDL_ERROR_VAL ("No input file provided",
				GDL_FAILURE,
				1);
	}
	if (OUTPUT == 0)
	{
		GDL_ERROR_VAL ("No output file provided",
				GDL_FAILURE,
				1);
	}
	if (GI2TAXID == 0)
	{
		GDL_ERROR_VAL ("No gi2taxid output file provided",
				GDL_FAILURE,
				1);
	}
	if (SIZE == 0)
	{
		GDL_ERROR_VAL ("No fragment size provided",
				GDL_FAILURE,
				1);
	}
        if (COVERAGE == 0)
        {
                GDL_ERROR_VAL ("No coverage value provided",
                                GDL_FAILURE,
                                1);
        }
	return GDL_SUCCESS;
}


static int
help (void)
{
	printf("%s - version 1.0\n", PROGRAM);
   	printf("Copyright (C) 2015 bioMerieux, France\n");
	printf ("\n");
	printf ("Draw random fragments of size k in each sequence of input multi-fasta or multi-fastq file, returning a new multi-fasta file\n");
	printf ("--help\tDisplay a brief help on program usage\n");
	printf ("--verbose\tOutput message on standard output to see what the program is doing\n");
	printf ("\n");
	printf ("--input or -i\tfile containing input sequences (fasta or fastq format)\n");
	printf ("--taxids or -t\tone-column file containing the taxid of each input sequence \n");
	printf ("--size or -l\tsize of drawn fragments \n");
	printf ("--coverage or -c\tmean coverage value for drawing fragments \n");
	printf ("--output or -o\toutput sequence file (fasta format)\n");
	printf ("--gi2taxid or -g\toutput gi2taxids file : two-colum file containing genome ids and taxids of the drawn fragments\n");	
	printf ("[--seed or -s\tvalue (integer) used to initialize the random seed (to use for repeatability purposes - default: seed based on time)]\n");
	printf ("[--atgc \tdraw fragments made of ATCG only]\n");
	printf ("   (warning : may take long (and even fail) if many non ATGC characters)\n");
	return GDL_SUCCESS;
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
		gzFile f;
		kseq_t *seq;
		FILE * stream, * stream_gi2taxid;
		FILE * taxids;
		gdl_string * line = 0;
		gdl_string * tx = 0;
		gdl_string *frag = 0;
		size_t i, j, k, random_pos, nfrag, n_write, to_write, atgc_only;
		int cptr;
		
		if(verbose_flag){
			printf("**** draw random fragments of size %d with mean coverage %f ****\n", SIZE, COVERAGE);
		}

		//------------------------//
		// initialize random seed //
		//------------------------//
		if(SEED == 0){
			srand(time(NULL));
		}else{
			srand(SEED);
		}

		//------------------//
		// open input files //
		// -----------------//
		f = gzopen(INPUT, "r");
		seq = kseq_init(f);
		taxids = gdl_fileopen (TAXIDS, "r");

		//--------------------//
		// open output files //
		//-------------------/
		stream = gdl_fileopen (OUTPUT, "w");
		stream_gi2taxid = gdl_fileopen (GI2TAXID, "w");

		//-----------------------//
		// process each sequence //
		//----------------------//
		cptr = 0; // sequence counter
		while(gdl_getline (&line, &n_write, taxids)!=-1)
		{
			i=j=0;
			// read taxid
			tx = gdl_string_next_token (line, n_write, &i, &j);
			//read corresponding sequence
			kseq_read(seq);
			// compute number of fragments to draw
			nfrag = COVERAGE*seq->seq.l / SIZE;
			// if sequence is shorter than fragment length, skip
			if(seq->seq.l <= SIZE){
				nfrag = 0;
				if(verbose_flag){
					printf("WARNING : skipping sequence %s : length = %lu < %d\n", seq->name.s, seq->seq.l, SIZE);
				}
			}
			// draw fragments
			for(k = 0; k < nfrag; k++){
				to_write = 0;
				while(to_write == 0){
					atgc_only = 1;
					// generate random position on the sequence
					random_pos = rand()%(seq->seq.l-SIZE);
					// write sub-seq into frag and flag non ATGC characters
					frag = gdl_string_alloc(SIZE);
					for(i = random_pos; i <= (random_pos + SIZE); i++){
						switch(toupper(seq->seq.s[i])){
							case 'A' : break;
							case 'T' : break;
							case 'G' : break;
							case 'C' : break;
							default :
								atgc_only = 0;
								break;
						}
						frag[i-random_pos] = toupper(seq->seq.s[i]);
					}
					// check wether to include fragment or not depending on the non ATGC characters
					if(atgc_only == 1){
						to_write = 1;
					}else{
						if(atgc_flag){
							to_write = 0;
						}else{
							to_write = 1;
							if(verbose_flag)
								printf("WARNING : sequence %s : drawing fragment (out of %lu required ones) containing IUPAC character(s)\n", seq->name.s, nfrag);
						}
					}
				}
				// write sequence header
				cptr++;
				fprintf (stream, ">%d\n", cptr);
				//write sequence
				for(j = 0; j < SIZE; j++)
				{
					if (j && j % 80 == 0)
					{
						fprintf (stream,"\n");
					}
					fputc (frag[j], stream);
				}
				fprintf (stream,"\n");
				// write gi2taxid information
				fprintf (stream_gi2taxid, "%s\t%s\n", seq->name.s, tx);
				// free fragment
				gdl_string_free(frag);
				frag = 0;
			}
		}
		//free
		gdl_string_free(tx);
		//free line
		gdl_string_free (line);
		// close input file
		kseq_destroy(seq);
		gzclose(f);
		gdl_fileclose (TAXIDS, taxids);
		//close files
		gdl_fileclose(OUTPUT, stream);
		gdl_fileclose(GI2TAXID, stream_gi2taxid);
		//free args
		gdl_string_free(INPUT);
		gdl_string_free(OUTPUT);
		gdl_string_free(TAXIDS);
	}
	exit(0);
}


