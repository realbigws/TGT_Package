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
#include "hhpred_hmm.h"
#include "hhpred_util.h"


/////////////////////////////////////////////////////////////////////////////////////
// // Describes an alignment of two profiles. Used as list element in Hits : List<Hit> 
/////////////////////////////////////////////////////////////////////////////////////
class Hit
{
 public:  
 
  //--- macro tag ----//
  char loc;                 // local or global alignment (default: 1)
  float mact;               // Score threshold for MAC alignment in local mode (set to 0.3501 to track user modification)
  float corr;               // Weight of correlations of scores for |i-j|<=4 (default: 0.1f);
  float shift;              // Score offset for match-match states (default: -0.03f);

	//--- data structure ---//
  float score;          // Score of alignment (i.e. of Viterbi path)
  float score_sort;     // score to sort hits in output list (negative means first/best!)
  float score_aass;     // first: just hit.score, then hit.logPval-SSSCORE2NATLOG*hit.score_ss;(negative means best!)
  double Pforward;      // scaled total forward probability : Pforward * Product_{i=1}^{Lq+1}(scale[i])
  
  int L;                // Number of match states in template
  int nsteps;           // index for last step in Viterbi path; (first=1)
  int wsteps;           // length of alignment
  int* wwi;
  int* wwj;
  int* wsi;             // wsi[step] = query match state at step of Viterbi path
  int* wsj;             // wsj[step] = template match state at step of Viterbi path
  char* states;         // state at step of Viterbi path  0: Start  1: M(MM)  2: A(-D)  3: B(IM)  4: C(D-)  5 D(MI)
  float* S;             // S[step] = match-match score contribution at alignment step
  float* P_posterior;   // P_posterior[step] = posterior prob for MM states (otherwise zero)
  int i1;               // First aligned residue in query
  int i2;               // Last aligned residue in query
  int j1;               // First aligned residue in template 
  int j2;               // Last aligned residue in template
  int matched_cols;     // number of matched columns in alignment against query
  int min_overlap;      // Minimum overlap between query and template
  float sum_of_probs;   // sum of probabilities for Maximum ACcuracy alignment (if dssp states defined, only aligned pairs with defined dssp state contribute to sum)

//--- temp data ---//
  //--- viterbi ---//
  float *Si;           // sMM[i][j] = score of best alignment up to indices (i,j) ending in (Match,Match)
  float *sMM;          // sMM[i][j] = score of best alignment up to indices (i,j) ending in (Match,Match)
  float *sGD;          // sGD[i][j] = score of best alignment up to indices (i,j) ending in (Gap,Delete)
  float *sDG;          // sDG[i][j] = score of best alignment up to indices (i,j) ending in (Delete,Gap)
  float *sIM;          // sIM[i][j] = score of best alignment up to indices (i,j) ending in (Ins,Match)
  float *sMI;          // sMI[i][j] = score of best alignment up to indices (i,j) ending in (Match,Ins)
  //---- forward ---//
  double *F_MM_prev;
  double *F_GD_prev;
  double *F_DG_prev;
  double *F_IM_prev;
  double *F_MI_prev;
  double *F_MM_curr;
  double *F_GD_curr;
  double *F_DG_curr;
  double *F_IM_curr;
  double *F_MI_curr;
  //---- backward ---//
  double *B_MM_prev;
  double *B_DG_prev;
  double *B_MI_prev;
  double *B_MM_curr;
  double *B_GD_curr;
  double *B_DG_curr;
  double *B_IM_curr;
  double *B_MI_curr;
  //---- macalign ----//
  double *S_prev;    // scores
  double *S_curr;    // scores
//--- temp data ---//over


  bool realign_around_viterbi;
  bool forward_allocated;
  bool backward_allocated;

	//--------- constructor and allocate ---------//
  // Constructor (only set pointers to NULL)
  Hit();
  ~Hit(){};
  

  // Allocate/delete memory for dynamic programming matrix
  void AllocateTempMatrix(int wsmax_);
  void AllocateBacktrace(int Nq, int Nt);
  void DeleteBacktrace();
  void AllocateBacktraceMatrix(int Nq, int Nt);
  void DeleteBacktraceMatrix(int Nq);
  void AllocateForwardMatrix(int Nq, int Nt);
  void DeleteForwardMatrix(int Nq);
  void AllocateBackwardMatrix(int Nq, int Nt);
  void DeleteBackwardMatrix(int Nq);


	//--------- alignment function ------------//
  // Compare an HMM with overlapping subalignments
  void Viterbi(HMM* q, HMM* t);
  // Compare two HMMs with each other in lin space
  void Forward(HMM* q, HMM* t);
  // Compare two HMMs with each other in lin space
  void Backward(HMM* q, HMM* t);
  // Find maximum accuracy alignment (after running Forward and Backward algorithms)
  void MACAlignment(HMM* q, HMM* t);
  // Trace back alignment of two profiles based on matrices bXX[][]
  void Backtrace(HMM* q, HMM* t);
  // Trace back MAC alignment of two profiles based on matrix bMM[][]
  void BacktraceMAC(HMM* q, HMM* t);


	//--------- score function ------------//
  // Calculate in log2 space the amino acid similarity score between columns i and j of two HMMs (query and template)
  inline float Score(float* qi, float* tj);
  // Calculate in lin space the amino acid similarity score between columns i and j of two HMMs (query and template)
  inline float ProbFwd(float* qi, float* tj);
  // Calculate score for a given alignment
  void ScoreAlignment(HMM* q, HMM* t, int steps);
  // Comparison (used to sort list of hits)
  int operator<(const Hit& hit2)  {return score_sort<hit2.score_sort;}

  
private:
  char state;          // 0: Start/STOP state  1: MM state  2: GD state (-D)  3: IM state  4: DG state (D-)  5 MI state
  char** bMM;          // (backtracing) bMM[i][j] = STOP:start of alignment  MM:prev was MM  GD:prev was GD etc
  char** bGD;          // (backtracing) bMM[i][j] = STOP:start of alignment  MM:prev was MM  SAME:prev was GD
  char** bDG;          // (backtracing)
  char** bIM;          // (backtracing)
  char** bMI;          // (backtracing)
  char** cell_off;     // cell_off[i][j]=1 means this cell will get score -infinity

  double** F_MM;        // Forward matrices 
  double** B_MM;        // Backward matrices
  double* scale;        // 

  void InitializeForAlignment(HMM* q, HMM* t);
};


double Pvalue(double x, double a[]);
double Pvalue(float x, float lamda, float mu);
double logPvalue(float x, float lamda, float mu);
double logPvalue(float x, double a[]);
double Probab(Hit& hit);




