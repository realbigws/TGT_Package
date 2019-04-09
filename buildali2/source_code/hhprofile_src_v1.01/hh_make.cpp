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
#include "hhpred_hmm.h"
#include "hhpred_ali.h"
#include "hhpred_util.h"
using namespace std;



//-------- main process --------//
void HH_Make(string &infile_,string &outfile_, string &inname, 
  int pseudo,int verbose,int par_M,int par_Mgaps,
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

  // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
  HMM* q = new HMM(MAXRES,verbose);          //Create a HMM with maximum of par.maxres match states
  qali.FrequenciesAndTransitions(q);
  if(inname!="")strcpy(q->name,inname.c_str());

  // Calculate pseudo-count
  if(pseudo==0)q->has_pseudocounts=true;
  q->AddTransitionPseudocounts();
  q->PreparePseudocounts();
  q->AddAminoAcidPseudocounts();
//  q->CalculateAminoAcidBackground();
//  q->Log2LinTransitionProbs(1.0);


  // Write filtered alignment WITH insert states (lower case) to alignment file
  char outfile[NAMELEN];
  strcpy(outfile,outfile_.c_str());
  q->WriteToFile(outfile);

  // Final delete
  delete q;
}


/////////////////////////////////////////////////////////////////////////////////////
// Exit function
/////////////////////////////////////////////////////////////////////////////////////
void Usage()
{
  printf("\n");
  printf("HHmake \n");
  printf("Build an HMM from an input alignment in A2M, A3M, or FASTA format,   \n");
  printf("\n");
  printf("Usage: ./hhmake -i file [options]                                       \n");
  printf(" -i <file>     query alignment (A2M, A3M, or FASTA)                     \n");
  printf(" -o <file>     HMM file to be written                                   \n");
  printf("[note]: the first sequence in input file should be the query sequence.  \n");
  printf("\n");
  printf("Output options:                                                              \n");
  printf(" -v <int>      verbose mode: 0:no screen output  [1]:only warings  2: verbose  \n");
  printf(" -a <name>     use this name for HMM (default: use name of first sequence)   \n");
  printf(" -s <0/1>      consider pseudocount for hhm file [0]:no pseidocount  1:consider \n");
  printf("\n");
  printf("Filter query multiple sequence alignment                                     \n");
  printf(" -b    [0,100]  maximum pairwise sequence identity (%%) (def=90)   \n");
  printf(" -d    [0,inf[  filter MSA by selecting most diverse set of sequences, keeping \n");
  printf("                at least this many seqs in each MSA block of length 50 (def=100) \n");
  printf(" -c    [0,100]  minimum coverage with query (%%) (def=0) \n");
  printf(" -p    [0,100]  minimum sequence identity with query (%%) (def=0) \n");
  printf(" -q   [-20,100] minimum score per column with query  (def=-20.0)\n");
  printf(" -n    [1,inf]  target diversity of alignment (default=off)\n");
  printf("\n");
  printf("Input alignment format:                                                    \n");
  printf(" -m a2m        use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf("               '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -m first      use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -m [0,100]    use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf("                                                                          \n");
  printf("Example: ./hhmake -i test.a3m -o test.hhm                                \n\n");

/*
  //-------- we DO NOT allow pseudocount for hhm file !! ------------//
  if (all)
  {
	  printf("Pseudocount (pc) options:                                                        \n");
	  printf(" -pcm  0-2      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",par.pcm);
	  printf("                0: no pseudo counts:    tau = 0                                  \n");
	  printf("                1: constant             tau = a                                  \n");
	  printf("                2: diversity-dependent: tau = a/(1 + ((Neff[i]-1)/b)^c)          \n");
	  printf("                (Neff[i]: number of effective seqs in local MSA around column i) \n");
	  printf("                3: constant diversity pseudocounts                               \n");
	  printf(" -pca  [0,1]    overall pseudocount admixture (def=%-.1f)                        \n",par.pca);
	  printf(" -pcb  [1,inf[  Neff threshold value for -pcm 2 (def=%-.1f)                      \n",par.pcb);
	  printf(" -pcc  [0,3]    extinction exponent c for -pcm 2 (def=%-.1f)                     \n",par.pcc);
	  printf(" -pre_pca [0,1]   PREFILTER pseudocount admixture (def=%-.1f)                    \n",par.pre_pca);
	  printf(" -pre_pcb [1,inf[ PREFILTER threshold for Neff (def=%-.1f)                       \n",par.pre_pcb);
	  printf("\n");
	  printf("Context-specific pseudo-counts:                                                  \n");
	  printf(" -nocontxt      use substitution-matrix instead of context-specific pseudocounts \n");
	  printf(" -contxt <file> context file for computing context-specific pseudocounts (default=%s)\n",par.clusterfile);
	  printf(" -cslib  <file> column state file for fast database prefiltering (default=%s)\n",par.cs_library);
	  printf("\n");
  }
*/
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
		string input_name = "";
		int pseudo=0;
		int verbose=1;
		int max_seqid=90;
		int Ndiff=100;    //-> this parameter is different from hh_filter.cpp
		int coverage=0;
		int qid=0;
		float qsc=-20.0f;
		float par_Neff=0.0f;
		int par_M=1;
		int par_Mgaps=50;

		//---- process argument ----//
		char c = 0;
		extern char* optarg;
		while ((c = getopt(argc, argv, "i:o:v:a:s:b:d:c:p:q:n:m:")) != EOF) 
		{
			switch (c) {
			case 'i':
				input_file = optarg;
				break;
			case 'o':
				output_file = optarg;
				break;
			case 'a':
				input_name = optarg;
				break;
			case 's':
				pseudo = atoi(optarg);
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
		if(pseudo<0 || pseudo>1)
		{
			fprintf(stderr," consider pseudocount for hhm file [0]:no pseidocount  1:consider \n");
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
		HH_Make(input_file,output_file,input_name,pseudo,verbose,par_M,par_Mgaps,
			max_seqid,coverage,qid,qsc,Ndiff,par_Neff);
		exit(0);
	}
}
