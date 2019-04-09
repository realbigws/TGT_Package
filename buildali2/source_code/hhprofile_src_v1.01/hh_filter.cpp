#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <getopt.h>
#include <unistd.h>
#include <fstream>
#include <string>
#include <string.h>
#include "hhpred_ali.h"
#include "hhpred_util.h"
using namespace std;



//-------- main process --------//
void HH_Filter(string &infile_,string &outfile_,int verbose,int par_M,int par_Mgaps,
	int max_seqid,int coverage,int qid,float qsc,int Ndiff,float par_Neff)
{
  int par_matrix=0;
  int par_nseqdis=MAXSEQ-1;
  Alignment qali(MAXSEQ,MAXRES,par_matrix,verbose);              //Create an alignment 
 
  // Reads in an alignment from par.infile into matrix X[k][l] as ASCII
  char infile[NAMELEN];
  strcpy(infile,infile_.c_str());
  FILE* inf = fopen(infile, "r");
  qali.Read(inf,infile,par_nseqdis);
  fclose(inf);

  // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i] 
  // and store marked sequences in name[k] and seq[k]
  qali.Compress(infile,par_M,par_Mgaps);

  // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
  qali.N_filtered = qali.Filter(max_seqid,coverage,qid,qsc,Ndiff);
  
  // Atune alignment diversity q.Neff with qsc to value Neff_goal
  if (par_Neff>=0.999) qali.FilterNeff(par_Neff,coverage,max_seqid);


  // Write filtered alignment WITH insert states (lower case) to alignment file
  char outfile[NAMELEN];
  strcpy(outfile,outfile_.c_str());
  qali.WriteToFile(outfile); 
}


/////////////////////////////////////////////////////////////////////////////////////
// Exit function
/////////////////////////////////////////////////////////////////////////////////////
void Usage()
{
  printf("\n");
  printf("HHfilter\n");
  printf("Filter an alignment by maximum sequence identity of match states and minimum coverage\n");
  printf("\n");
  printf("Usage: ./hhfilter -i infile -o outfile [options]                  \n");
  printf(" -i <file>      read input file in A3M/A2M or FASTA format                 \n");
  printf(" -o <file>      write to output file in A3M format                         \n");
  printf("\n");
  printf("Options:                                                                  \n");
  printf(" -v <int>       verbose mode: 0:no screen output  [1]:only warings  2: verbose\n");
  printf(" -b    [0,100]  maximum pairwise sequence identity (%%) (def=90)   \n");
  printf(" -d    [0,inf[  filter MSA by selecting most diverse set of sequences, keeping \n");
  printf("                at least this many seqs in each MSA block of length 50 (def=0) \n");
  printf(" -c    [0,100]  minimum coverage with query (%%) (def=0) \n");
  printf(" -p    [0,100]  minimum sequence identity with query (%%) (def=0) \n");
  printf(" -q   [-20,100] minimum score per column with query  (def=-20.0)\n");
  printf(" -n    [0,inf[  target diversity of alignment (default=off)\n");
  printf("\n");         
  printf("Input alignment format:                                                    \n");
  printf(" -m a2m         use A2M/A3M (default): upper case = Match; lower case = Insert;\n");         
  printf("                '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -m first       use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -m [0,100]     use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf("                                                                          \n");
  printf("Example: ./hhfilter -b 50 -i d1mvfd_.a2m -o d1mvfd_.fil.a2m          \n\n");
}

//------ main -------//
int main(int argc,char **argv)
{
	//---- HH_Filter ----//ws_version
	{
		if(argc < 3)
		{
			Usage();
			exit(-1);
		}
		string input_file = "";
		string output_file = "";
		int verbose=1;
		int max_seqid=90;
		int Ndiff=0;
		int coverage=0;
		int qid=0;
		float qsc=-20.0f;
		float par_Neff=0.0f;
		int par_M=1;
		int par_Mgaps=50;

		//---- process argument ----//
		char c = 0;
		extern char* optarg;
		while ((c = getopt(argc, argv, "i:o:v:b:d:c:p:q:n:m:")) != EOF) 
		{
			switch (c) {
			case 'i':
				input_file = optarg;
				break;
			case 'o':
				output_file = optarg;
				break;
			case 'v':
				verbose = atoi(optarg);
				break;
			case 'b':
				max_seqid = atoi(optarg);
				break;
			case 'd':
				Ndiff = atoi(optarg);
				break;
			case 'c':
				coverage = atoi(optarg);
				break;
			case 'p':
				qid = atoi(optarg);
				break;
			case 'q':
				qsc = atof(optarg);
				break;
			case 'n':
				par_Neff = atof(optarg);
				break;
			case 'm':
				if(optarg=="a2m" || optarg=="a3m")par_M=1;
				else if(optarg=="first")par_M=3;
				else
				{
					if (optarg[0]>='0' && optarg[0]<='9')
					{
						par_Mgaps=atoi(optarg); 
						par_M=2;
					}
					else
					{
						Usage();
						exit(-1);
					}
				}
				break;

			default:
				Usage();
				exit(-1);
			}
		}

		//---- check arguments -----//
		if(input_file=="")
		{
			fprintf(stderr,"input file is blank ! \n");
			exit(-1);
		}
		if(output_file=="")
		{
			fprintf(stderr,"output file is blank ! \n");
			exit(-1);
		}
		if(verbose<0 || verbose>2)
		{
			fprintf(stderr," verbose mode: 0:no screen output  [1]:only warings  2: verbose\n");
			exit(-1);
		}
		if(max_seqid<0 || max_seqid>100)
		{
			fprintf(stderr," -b    [0,100]  maximum pairwise sequence identity (%%) (def=90)   \n");
			exit(-1);
		}
		if(Ndiff<0)
		{
		  fprintf(stderr," -d    [0,inf]  filter MSA by selecting most diverse set of sequences, keeping \n");
		  fprintf(stderr,"                at least this many seqs in each MSA block of length 50 (def=0) \n");
		  exit(-1);
		}
		if(coverage<0 || coverage>100)
		{
			fprintf(stderr," -c    [0,100]  minimum coverage with query (%%) (def=0) \n");
			exit(-1);
		}
		if(qid<0 || qid>100)
		{
			fprintf(stderr," -p    [0,100]  minimum sequence identity with query (%%) (def=0) \n");
			exit(-1);
		}
		if(qsc<-20 || qsc>100)
		{
			fprintf(stderr," -q   [-20,100] minimum score per column with query  (def=-20.0)\n");
			exit(-1);
		}
		if(par_Neff<0)
		{
			fprintf(stderr," -n    [0,inf]  target diversity of alignment (default=off)\n");
			exit(-1);
		}

		//----- main process -----//
		HH_Filter(input_file,output_file,verbose,par_M,par_Mgaps,
			max_seqid,coverage,qid,qsc,Ndiff,par_Neff);
		exit(0);
	}
}
