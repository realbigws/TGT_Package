#pragma once
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <stdio.h>    // printf
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <time.h>     // clock
#include <malloc.h>
#include "hhpred_util.h"


//========== class HMM =========//
class HMM
{
public:
	int v;
	int maxres;
    HMM(int maxres_=MAXRES,int v_=1);
    ~HMM();
    HMM& operator=(HMM&);

    //--- macro tag ----//
    bool has_pseudocounts;    // set to true if HMM contains pseudocounts

		//----------- main variables ----------//
    int L;                    // length of HMM = number of match states; set in declaration of HMM object
    char* seq;                // residues of stored sequences (first at pos 1!)
    char* name;               // Full name of first sequence of original alignment (NAME field)
    float* Neff_M;            // Neff_M[i] = diversity of subalignment of seqs that have residue in col i
    float* Neff_I;            // Neff_I[i] = diversity of subalignment of seqs that have insert in col i
    float* Neff_D;            // Neff_D[i] = diversity of subalignment of seqs that have delete in col i
    float Neff_HMM;           // average number of Neff over total length of HMM
    float pb[NAA];            // pb[a] = background freq of amino acids
    float pav[NAA];           // pav[a] = average freq of amino acids in HMM (including subst matrix pseudocounts)
    int* l;                   // l[i] = pos. of j'th match state in aligment
    float** p;                // p[i][a] = prob of finding amino acid a in column i WITH OPTIMUM pseudocounts
    float** f;                // f[i][a] = prob of finding amino acid a in column i WITHOUT pseudocounts
    float** g;                // g[i][a] = prob of finding amino acid a in column i WITH pseudocounts
    float** tr;               // tr[i][X2Y] = log2 of transition probabilities M2M M2I M2D I2M I2I D2M D2D  
    float** trr;              // tr in lin space //__130101__//by WS

		//------------ main functions ----------//
    // Read an HMM from a HHsearch .hhm file and return 0 at end of file
    int Read(FILE* dbf);
    int Read_HMM(FILE* dbf,int L);
    int Read_TGT(FILE* dbf);
    int Read_TPL(FILE* dbf);

    // Add transition pseudocounts to HMM
    void AddTransitionPseudocounts(float gapd=0.15, float gape=1.0, float gapf=0.6, float gapg=0.6, float gaph=0.6, float gapi=0.6, float gapb=1.0);

    // Generate an amino acid frequency matrix g[][] with full pseudocount admixture (tau=1)
    void PreparePseudocounts();

    // Add amino acid pseudocounts to HMM: t.p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
    void AddAminoAcidPseudocounts(char pcm=2, float pca=1.0f, float pcb=1.5f, float pcc=1.0f);

    // Calculate amino acid backround frequencies for HMM
    void CalculateAminoAcidBackground();

    // Add no amino acid pseudocounts to HMM: copy  t.p[i][a] = f[i][a]
    void NoAminoAcidPseudocounts() {for(int i=1; i<=L; i++) for(int a=0; a<NAA; a++) p[i][a]=f[i][a];};

    // Factor Null model into HMM t
    void IncludeNullModelInHMM(HMM* q, HMM* t, int columnscore=1);

    // Transform log to lin transition probs
    void Log2LinTransitionProbs(float beta=1.0);

    // Set query columns in His-tags etc to Null model distribution
    void NeutralizeTags();

    // Calculate effective number of sequences using profiles INCLUDING pseudocounts
    float CalcNeff();

    // Write HMM to output file
    void WriteToFile(char* outfile);
    //------------ main functions ----------//over

    // Utility for Read()
    int Warning(FILE* dbf, char line[], char name[])
    {
        if (v) std::cerr<<"\nWARNING: could not read line\n\'"<<line<<"\'\nin HMM "<<name<<"\n";
//        while (fgetline(line,LINELEN,dbf) && !(line[0]=='/' && line[1]=='/'));
        if (line) return 2;  //return status: skip HMM
        return 0;            //return status: end of database file
    }
};
