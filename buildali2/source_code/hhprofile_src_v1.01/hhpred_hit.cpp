#include "hhpred_hit.h"


//------------- definition ------------//
#define CALCULATE_MAX6(max, var1, var2, var3, var4, var5, var6, varb)	\
  if (var1>var2) { max=var1; varb=STOP;}				\
  else           { max=var2; varb=MM;};					\
  if (var3>max)  { max=var3; varb=GD;};					\
  if (var4>max)  { max=var4; varb=IM;};					\
  if (var5>max)  { max=var5; varb=DG;};					\
  if (var6>max)  { max=var6; varb=MI;}; 

#define CALCULATE_MAX4(max, var1, var2, var3, var4, varb)	\
  if (var1>var2) { max=var1; varb=STOP;}			\
  else           { max=var2; varb=MM;};				\
  if (var3>max)  { max=var3; varb=MI;};				\
  if (var4>max)  { max=var4; varb=IM;}; 

// Generate random number in [0,1[
#define frand() ((float) rand()/(RAND_MAX+1.0))


// Function declarations
inline float max2(const float& xMM, const float& xX, char& b); 
inline int pickprob2(const double& xMM, const double& xX, const int& state); 
inline int pickprob3_GD(const double& xMM, const double& xDG, const double& xGD); 
inline int pickprob3_IM(const double& xMM, const double& xMI, const double& xIM); 
inline int pickprob6(const double& x0, const double& xMM, const double& xGD, const double& xIM, const double& xDG, const double& xMI); 
inline int pickmax2(const double& xMM, const double& xX, const int& state); 
inline int pickmax3_GD(const double& xMM, const double& xDG, const double& xGD); 
inline int pickmax3_IM(const double& xMM, const double& xMI, const double& xIM); 
inline int pickmax6(const double& x0, const double& xMM, const double& xGD, const double& xIM, const double& xDG, const double& xMI); 
inline double Pvalue(double x, double a[]);
inline double Pvalue(float x, float lamda, float mu);
inline double logPvalue(float x, float lamda, float mu);
inline double logPvalue(float x, double a[]);
inline double Probab(Hit& hit);

/////////////////////////////////////////////////////////////////////////////////////
//// Constructor
/////////////////////////////////////////////////////////////////////////////////////
Hit::Hit()
{
	//-- macro parameter --//
	mact=0.3501f;
	corr=0.1f;
	shift=-0.03f;

	//-- data structure --//
  bMM = bGD = bDG = bIM = bMI = NULL;
  cell_off = NULL;
  wsi = wsj = NULL;
  wwi = wwj = NULL;
  states = NULL;
  S = P_posterior = NULL;
  B_MM=NULL;
  F_MM=NULL;
  scale = NULL;
  sum_of_probs=0.0; 
  realign_around_viterbi=false;
}

//-------- allocate temp matrix ------//
void Hit::AllocateTempMatrix(int wsmax_)
{
  //--- viterbi ---//
  Si=new float[wsmax_];           // sMM[i][j] = score of best alignment up to indices (i,j) ending in (Match,Match)
  sMM=new float[wsmax_];          // sMM[i][j] = score of best alignment up to indices (i,j) ending in (Match,Match)
  sGD=new float[wsmax_];          // sGD[i][j] = score of best alignment up to indices (i,j) ending in (Gap,Delete)
  sDG=new float[wsmax_];          // sDG[i][j] = score of best alignment up to indices (i,j) ending in (Delete,Gap)
  sIM=new float[wsmax_];          // sIM[i][j] = score of best alignment up to indices (i,j) ending in (Ins,Match)
  sMI=new float[wsmax_];          // sMI[i][j] = score of best alignment up to indices (i,j) ending in (Match,Ins)

  //---- forward ---//
  F_MM_prev=new double[wsmax_];
  F_GD_prev=new double[wsmax_];
  F_DG_prev=new double[wsmax_];
  F_IM_prev=new double[wsmax_];
  F_MI_prev=new double[wsmax_];
  F_MM_curr=new double[wsmax_];
  F_GD_curr=new double[wsmax_];
  F_DG_curr=new double[wsmax_];
  F_IM_curr=new double[wsmax_];
  F_MI_curr=new double[wsmax_];

  //---- backward ---//
  B_MM_prev=new double[wsmax_];
  B_DG_prev=new double[wsmax_];
  B_MI_prev=new double[wsmax_];
  B_MM_curr=new double[wsmax_];
  B_GD_curr=new double[wsmax_];
  B_DG_curr=new double[wsmax_];
  B_IM_curr=new double[wsmax_];
  B_MI_curr=new double[wsmax_];

  //---- macalign ----//
  S_prev=new double[wsmax_];    // scores
  S_curr=new double[wsmax_];    // scores
}


//--------- allocate/delete backtrace matrix --------//
void Hit::AllocateBacktrace(int Nq, int Nt)
{
  wsi = new( int[Nq+Nt]);
  wsj = new( int[Nq+Nt]);
  wwi = new( int[Nq+Nt]);
  wwj = new( int[Nq+Nt]);
  states  = new( char[Nq+Nt]);
  S = new( float[Nq+Nt]);
  P_posterior = new( float[Nq+Nt]);
}
void Hit::DeleteBacktrace(void)
{
  delete [] wsi;
  delete [] wsj;
  delete [] wwi;
  delete [] wwj;
  delete [] states;
  delete [] S;
  delete [] P_posterior;
  wsi=wsj=NULL;
  wwi=wwj=NULL;
  states=NULL;
  S=P_posterior=NULL;
}


/////////////////////////////////////////////////////////////////////////////////////
//// Allocate/delete memory for dynamic programming matrix
/////////////////////////////////////////////////////////////////////////////////////
void Hit::AllocateBacktraceMatrix(int Nq, int Nt)
{
  int i;
  bMM=new(char*[Nq]);
  bMI=new(char*[Nq]);
  bIM=new(char*[Nq]);
  bDG=new(char*[Nq]);
  bGD=new(char*[Nq]);
  cell_off=new(char*[Nq]);
  for (i=0; i<Nq; ++i) 
    {
      bMM[i]=new(char[Nt]);
      bMI[i]=new(char[Nt]);
      bIM[i]=new(char[Nt]);
      bGD[i]=new(char[Nt]);
      bDG[i]=new(char[Nt]);
      cell_off[i]=new(char[Nt]);
      if (!bMM[i] || !bMI[i] || !bIM[i] || !bGD[i] || !bDG[i] || !cell_off[i]) 
	{
	  fprintf(stderr,"Error: out of memory while allocating row %i (out of %i) for dynamic programming matrices \n",i+1,Nq);
	  fprintf(stderr,"Please decrease your memory requirements to the available memory using option -maxmem <GBs>\n");
	  fprintf(stderr,"You may want to check and increase your stack size limit (Linux: ulimit -a)\n");
	  exit(3);
	} 
    }
}
void Hit::DeleteBacktraceMatrix(int Nq)
{
  int i;
  for (i=0; i<Nq; ++i) 
    {
      delete[] bMM[i];
      delete[] bMI[i];
      delete[] bIM[i];
      delete[] bGD[i];
      delete[] bDG[i];
      delete[] cell_off[i];
    }
  delete[] bMM;
  delete[] bMI;
  delete[] bIM;
  delete[] bDG;
  delete[] bGD;
  delete[] cell_off;
  bMM = bGD = bDG = bIM = bMI = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Allocate/delete memory for Forward dynamic programming matrix
/////////////////////////////////////////////////////////////////////////////////////
void Hit::AllocateForwardMatrix(int Nq, int Nt)
{
  F_MM=new(double*[Nq]);
  scale=new(double[Nq+3]); // need Nq+3?
  for (int i=0; i<Nq; ++i) 
    {
      F_MM[i] = new(double[Nt]);
      if (!F_MM[i]) 
	{
	  fprintf(stderr,"Error: out of memory while allocating row %i (out of %i) for dynamic programming matrices \n",i+1,Nq);
	  fprintf(stderr,"Please decrease your memory requirements to the available memory using option -maxmem <GBs>\n");
	  fprintf(stderr,"You may want to check and increase your stack size limit (Linux: ulimit -a)\n");
	  exit(3);
	} 
      for (int j=0; j<Nt; ++j) 
	F_MM[i][j]=0.0; // This might be time-consuming! Is it necessary???? JS
	
    }
  forward_allocated = true;
}
void Hit::DeleteForwardMatrix(int Nq)
{
  for (int i=0; i<Nq; ++i) 
    {
      delete[] F_MM[i];
    }
  delete[] F_MM;
  delete[] scale;
  F_MM=NULL;
  forward_allocated = false;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Allocate/delete memory for Backward dynamic programming matrix (DO ONLY AFTER FORWARD MATRIX HAS BEEN ALLOCATED)
/////////////////////////////////////////////////////////////////////////////////////
void Hit::AllocateBackwardMatrix(int Nq, int Nt)
{
  B_MM=new(double*[Nq]);
  for (int i=0; i<Nq; ++i) 
    {
      B_MM[i] = new(double[Nt]);
      if (!B_MM[i]) 
	{
	  fprintf(stderr,"Error: out of memory while allocating row %i (out of %i) for dynamic programming matrices \n",i+1,Nq);
	  fprintf(stderr,"Please decrease your memory requirements to the available memory using option -maxmem <GBs>\n");
	  fprintf(stderr,"You may want to check and increase your stack size limit (Linux: ulimit -a)\n");
	  exit(3);
	} 
      for (int j=0; j<Nt; ++j) 
	B_MM[i][j]=0.0;   // This might be time-consuming! Is it necessary???? JS
    }
  backward_allocated = true;
}
void Hit::DeleteBackwardMatrix(int Nq)
{
  for (int i=0; i<Nq; ++i) 
    {
      delete[] B_MM[i];
    }
  delete[] B_MM;
  B_MM=NULL;
  backward_allocated = false;
}


/////////////////////////////////////////////////////////////////////////////////////
// Compare HMMs with one another and look for sub-optimal alignments that share no pair with previous ones
// The function is called with q and t
/////////////////////////////////////////////////////////////////////////////////////
void Hit::Viterbi(HMM* q, HMM* t)
{
  
  // Linear topology of query (and template) HMM:
  // 1. The HMM HMM has L+2 columns. Columns 1 to L contain 
  //    a match state, a delete state and an insert state each.
  // 2. The Start state is M0, the virtual match state in column i=0 (j=0). (Therefore X[k][0]=ANY)
  //    This column has only a match state and it has only a transitions to the next match state.
  // 3. The End state is M(L+1), the virtual match state in column i=L+1.(j=L+1) (Therefore X[k][L+1]=ANY)
  //    Column L has no transitions to the delete state: tr[L][M2D]=tr[L][D2D]=0.
  // 4. Transitions I->D and D->I are ignored, since they do not appear in PsiBlast alignments 
  //    (as long as the gap opening penalty d is higher than the best match score S(a,b)). 
  
  // Pairwise alignment of two HMMs:
  // 1. Pair-states for the alignment of two HMMs are 
  //    MM (Q:Match T:Match) , GD (Q:Gap T:Delete), IM (Q:Insert T:Match),  DG (Q:Delelte, T:Match) , MI (Q:Match T:Insert) 
  // 2. Transitions are allowed only between the MM-state and each of the four other states.
  
  // Saving space:
  // The best score ending in pair state XY sXY[i][j] is calculated from left to right (j=1->t->L) 
  // and top to bottom (i=1->q->L). To save space, only the last row of scores calculated is kept in memory.
  // (The backtracing matrices are kept entirely in memory [O(t->L*q->L)]).
  // When the calculation has proceeded up to the point where the scores for cell (i,j) are caculated,
  //    sXY[i-1][j'] = sXY[j']   for j'>=j (A below)  
  //    sXY[i][j']   = sXY[j']   for j'<j  (B below)
  //    sXY[i-1][j-1]= sXY_i_1_j_1         (C below) 
  //    sXY[i][j]    = sXY_i_j             (D below)
  //                   j-1   
  //                     j
  // i-1:               CAAAAAAAAAAAAAAAAAA
  //  i :   BBBBBBBBBBBBBD
  
  
  // Variable declarations
//  float Si[MAXRES];           // sMM[i][j] = score of best alignment up to indices (i,j) ending in (Match,Match) 
//  float sMM[MAXRES];          // sMM[i][j] = score of best alignment up to indices (i,j) ending in (Match,Match) 
//  float sGD[MAXRES];          // sGD[i][j] = score of best alignment up to indices (i,j) ending in (Gap,Delete) 
//  float sDG[MAXRES];          // sDG[i][j] = score of best alignment up to indices (i,j) ending in (Delete,Gap)
//  float sIM[MAXRES];          // sIM[i][j] = score of best alignment up to indices (i,j) ending in (Ins,Match)
//  float sMI[MAXRES];          // sMI[i][j] = score of best alignment up to indices (i,j) ending in (Match,Ins) 
  float smin=(loc? 0:-FLT_MAX);  //used to distinguish between SW and NW algorithms in maximization         
  int i,j;      //query and template match state indices
  float sMM_i_j=0,sMI_i_j,sIM_i_j,sGD_i_j,sDG_i_j;
  float sMM_i_1_j_1,sMI_i_1_j_1,sIM_i_1_j_1,sGD_i_1_j_1,sDG_i_1_j_1;
  int jmin,jmax;

  // Reset crossed out cells?
  InitializeForAlignment(q,t);

  // Initialization of top row, i.e. cells (0,j)
  for (j=0; j<=t->L; ++j) 
    {
      sMM[j] = 0;
      sIM[j] = sMI[j] = sDG[j] = sGD[j] = -FLT_MAX; 
    }
  score=-INT_MAX; i2=j2=0; bMM[0][0]=STOP;

  // Viterbi algorithm
  for (i=1; i<=q->L; ++i) // Loop through query positions i
    {
      
	{
	  // If q is compared to t, exclude regions where overlap of q with t < min_overlap residues
	  jmin=imax( 1, i+min_overlap-q->L);  // Lq-i+j>=Ovlap => j>=i+Ovlap-Lq => jmin=max{1, i+Ovlap-Lq} 
	  jmax=imin(t->L,i-min_overlap+t->L);  // Lt-j+i>=Ovlap => j<=i-Ovlap+Lt => jmax=min{Lt,i-Ovlap+Lt} 
	}      

      // Initialize cells
      if (jmin==1) 
	{
	  sMM_i_1_j_1 = -(i-1)*0;  // initialize at (i-1,0)   // egq -> default:0  // no charge for end gaps as default
	  sMM[0] = -i*0;           // initialize at (i,0)  // egq -> default:0
	  sIM_i_1_j_1 = sMI_i_1_j_1 = sDG_i_1_j_1 = sGD_i_1_j_1 = -FLT_MAX; // initialize at (i-1,jmin-1)
	} 
      else 
	{
	  // Initialize at (i-1,jmin-1) if lower left triagonal is excluded due to min overlap
	  sMM_i_1_j_1 = sMM[jmin-1];     // initialize at (i-1,jmin-1)
	  sIM_i_1_j_1 = sIM[jmin-1];     // initialize at (i-1,jmin-1)
	  sMI_i_1_j_1 = sMI[jmin-1];     // initialize at (i-1,jmin-1)
	  sDG_i_1_j_1 = sDG[jmin-1];     // initialize at (i-1,jmin-1)
	  sGD_i_1_j_1 = sGD[jmin-1];     // initialize at (i-1,jmin-1)
	  sMM[jmin-1] = -FLT_MAX;        // initialize at (i,jmin-1)
	}
      if (jmax<t->L) // initialize at (i-1,jmmax) if upper right triagonal is excluded due to min overlap
	sMM[jmax] = sIM[jmax] = sMI[jmax] = sDG[jmax] = sGD[jmax] = -FLT_MAX; 
      sIM[jmin-1] = sMI[jmin-1] = sDG[jmin-1] = sGD[jmin-1] = -FLT_MAX; // initialize at (i,jmin-1)

      // Precalculate amino acid profile-profile scores
      for (j=jmin; j<=jmax; ++j) 
	Si[j] = Score(q->p[i],t->p[j]);      

      for (j=jmin; j<=jmax; ++j) // Loop through template positions j
	{

	  if (cell_off[i][j])
	    {
	      sMM_i_1_j_1 = sMM[j]; // sMM_i_1_j_1 (for j->j+1) = sMM(i-1,(j+1)-1) = sMM[j] 
	      sGD_i_1_j_1 = sGD[j];
	      sIM_i_1_j_1 = sIM[j];
	      sDG_i_1_j_1 = sDG[j];
	      sMI_i_1_j_1 = sMI[j];
	      sMM[j]=sMI[j]=sIM[j]=sDG[j]=sGD[j]=-FLT_MAX; // sMM[j] = sMM(i,j) is cell_off
	    }
	  else 
	    {
	      // Recursion relations

	      CALCULATE_MAX6( sMM_i_j,
			      smin,
			      sMM_i_1_j_1 + q->tr[i-1][M2M] + t->tr[j-1][M2M], 
			      sGD_i_1_j_1 + q->tr[i-1][M2M] + t->tr[j-1][D2M],
			      sIM_i_1_j_1 + q->tr[i-1][I2M] + t->tr[j-1][M2M],
			      sDG_i_1_j_1 + q->tr[i-1][D2M] + t->tr[j-1][M2M],
			      sMI_i_1_j_1 + q->tr[i-1][M2M] + t->tr[j-1][I2M],
			      bMM[i][j]
			      );

 	      sMM_i_j += Si[j] + shift ; 
	      

	      sGD_i_j = max2
		(
		 sMM[j-1] + t->tr[j-1][M2D], // MM->GD gap opening in query 
		 sGD[j-1] + t->tr[j-1][D2D], // GD->GD gap extension in query 
		 bGD[i][j]
		 );
	      sIM_i_j = max2
		(
		 sMM[j-1] + q->tr[i][M2I] + t->tr[j-1][M2M] ,
		 sIM[j-1] + q->tr[i][I2I] + t->tr[j-1][M2M], // IM->IM gap extension in query 
		 bIM[i][j]
		 );
	      sDG_i_j = max2
		(
		 sMM[j] + q->tr[i-1][M2D],
		 sDG[j] + q->tr[i-1][D2D], //gap extension (DD) in query
		 bDG[i][j]
		 );
	      sMI_i_j = max2
		(
		 sMM[j] + q->tr[i-1][M2M] + t->tr[j][M2I], // MM->MI gap opening M2I in template 
		 sMI[j] + q->tr[i-1][M2M] + t->tr[j][I2I], // MI->MI gap extension I2I in template 
		 bMI[i][j]
		 );

	      sMM_i_1_j_1 = sMM[j];
	      sGD_i_1_j_1 = sGD[j];
	      sIM_i_1_j_1 = sIM[j];
	      sDG_i_1_j_1 = sDG[j];
	      sMI_i_1_j_1 = sMI[j];
	      sMM[j] = sMM_i_j;
	      sGD[j] = sGD_i_j;
	      sIM[j] = sIM_i_j;
	      sDG[j] = sDG_i_j;
	      sMI[j] = sMI_i_j;
	      
	      // Find maximum score; global alignment: maxize only over last row and last column
	      if(sMM_i_j>score && (loc || i==q->L)) { i2=i; j2=j; score=sMM_i_j; }

	    } // end if 

	} //end for j
      
      // if global alignment: look for best cell in last column
      if (!loc && sMM_i_j>score) { i2=i; j2=jmax; score=sMM_i_j; }
      
    } // end for i

  state=MM; // state with maximum score is MM state

  //   printf("Template=%-12.12s  i=%-4i j=%-4i score=%6.3f\n",t->name,i2,j2,score);

  return;
}



/////////////////////////////////////////////////////////////////////////////////////
// Compare two HMMs with Forward Algorithm in lin-space (~ 2x faster than in log-space)
/////////////////////////////////////////////////////////////////////////////////////
void Hit::Forward(HMM* q, HMM* t)
{

  // Variable declarations
  int i,j;      // query and template match state indices
  double pmin=(loc? 1.0: 0.0);    // used to distinguish between SW and NW algorithms in maximization         
  double Cshift = pow(2.0,shift);   // score offset transformed into factor in lin-space
  double Pmax_i;                        // maximum of F_MM in row i
  double scale_prod=1.0;                // Prod_i=1^i (scale[i])
  int jmin; 

//  double F_MM_prev[MAXRES];
//  double F_GD_prev[MAXRES];
//  double F_DG_prev[MAXRES];
//  double F_IM_prev[MAXRES];
//  double F_MI_prev[MAXRES];

//  double F_MM_curr[MAXRES];
//  double F_GD_curr[MAXRES];
//  double F_DG_curr[MAXRES];
//  double F_IM_curr[MAXRES];
//  double F_MI_curr[MAXRES];

  // First alignment of this pair of HMMs?
    {
      q->trr[0][M2D] = q->trr[0][M2I] = 0.0;
      q->trr[0][I2M] = q->trr[0][I2I] = 0.0;
      q->trr[0][D2M] = q->trr[0][D2D] = 0.0;
      t->trr[0][M2M] = 1.0;
      t->trr[0][M2D] = t->trr[0][M2I] = 0.0;
      t->trr[0][I2M] = t->trr[0][I2I] = 0.0;
      t->trr[0][D2M] = t->trr[0][D2D] = 0.0;
      q->trr[q->L][M2M] = 1.0;
      q->trr[q->L][M2D] = q->trr[q->L][M2I] = 0.0;
      q->trr[q->L][I2M] = q->trr[q->L][I2I] = 0.0;
      q->trr[q->L][D2M] = 1.0;
      q->trr[q->L][D2D] = 0.0;
      t->trr[t->L][M2M] = 1.0;
      t->trr[t->L][M2D] = t->trr[t->L][M2I] = 0.0;
      t->trr[t->L][I2M] = t->trr[t->L][I2I] = 0.0;
      t->trr[t->L][D2M] = 1.0;
      t->trr[t->L][D2D] = 0.0;
      InitializeForAlignment(q,t);
    }	

  if (realign_around_viterbi)
    {
      int step;
      
      // Cross out regions
      for (i=1; i<=q->L; ++i) 
	for (j=1; j<=t->L; ++j) 
	  if (!((i < i1 && j < j1) || (i > i2 && j > j2)))
	    cell_off[i][j]=1;   
      
      // Clear Viterbi path
      for (step=nsteps; step>=1; step--)
	{
	  int path_width=40;
	  for (i=imax(1,abs(wsi[step])-path_width); i <= imin(q->L,abs(wsi[step])+path_width); ++i)
	    cell_off[i][abs(wsj[step])]=0;
	  for (j=imax(1,abs(wsj[step])-path_width); j <= imin(t->L,abs(wsj[step])+path_width); ++j)
	    cell_off[abs(wsi[step])][j]=0;
	}
    }


  for (j=1; j<=t->L; ++j) 
    F_MM_curr[j] = 0.0;

  F_MM_curr[0] = 0.0;
  F_IM_curr[0] = 0.0;
  F_GD_curr[0] = 0.0;
  for (j=1; j<=t->L; ++j)
  {
    if (cell_off[1][j])
      F_MM_curr[j] = F_MI_curr[j] = F_DG_curr[j] = F_IM_curr[j] = F_GD_curr[j] = 0.0;
    else
    {
      F_MM_curr[j] = ProbFwd(q->p[1],t->p[j]) * Cshift  ;
      F_MI_curr[j] = F_DG_curr[j] = 0.0;
      F_IM_curr[j] = F_MM_curr[j-1] * q->trr[1][M2I] * t->trr[j-1][M2M] + F_IM_curr[j-1] * q->trr[1][I2I] * t->trr[j-1][M2M];
      F_GD_curr[j] = F_MM_curr[j-1] * t->trr[j-1][M2D]                + F_GD_curr[j-1] * t->trr[j-1][D2D];
    }
  }

  for (int j = 0; j <= t->L; j++)
  {
    F_MM[0][j] = F_MM_prev[j];
    F_MM[1][j] = F_MM_curr[j];

    F_MM_prev[j] = F_MM_curr[j];
    F_MI_prev[j] = F_MI_curr[j];
    F_IM_prev[j] = F_IM_curr[j];
    F_DG_prev[j] = F_DG_curr[j];
    F_GD_prev[j] = F_GD_curr[j];
  }


  scale[0]=scale[1]=scale[2]=1.0;

  // Forward algorithm
  for (i=2; i<=q->L; ++i) // Loop through query positions i
    {
      jmin=1;

      if (scale_prod<DBL_MIN*100) scale_prod = 0.0; else scale_prod *= scale[i];

      // Initialize cells at (i,0)
      if (cell_off[i][jmin]) 
	F_MM_curr[jmin] = F_MI_curr[jmin] = F_DG_curr[jmin] = F_IM_curr[jmin] = F_GD_curr[jmin] = 0.0;
      else 
      {
	F_MM_curr[jmin] = scale_prod * ProbFwd(q->p[i],t->p[jmin]) * Cshift ;
	F_IM_curr[jmin] = F_GD_curr[jmin] = 0.0; 
	F_MI_curr[jmin] = scale[i] * (F_MM_prev[jmin] * q->trr[i-1][M2M] * t->trr[jmin][M2I] + F_MI_prev[jmin] * q->trr[i-1][M2M] * t->trr[jmin][I2I]);
	F_DG_curr[jmin] = scale[i] * (F_MM_prev[jmin] * q->trr[i-1][M2D]                   + F_DG_prev[jmin] * q->trr[i-1][D2D]);
      }

      /* copy back */
      F_MM[i][jmin] = F_MM_curr[jmin];

      Pmax_i=0;
 
      for (j=jmin+1; j<=t->L; ++j) // Loop through template positions j
      {
	// Recursion relations

	if (cell_off[i][j]) 
	  F_MM_curr[j] = F_MI_curr[j] = F_DG_curr[j] = F_IM_curr[j] = F_GD_curr[j] = 0.0;
	else
	{
	  F_MM_curr[j] = ProbFwd(q->p[i],t->p[j]) * Cshift  * scale[i] *
	    ( pmin
	      + F_MM_prev[j-1] * q->trr[i-1][M2M] * t->trr[j-1][M2M] // BB -> MM (BB = Begin/Begin, for local alignment)
	      + F_GD_prev[j-1] * q->trr[i-1][M2M] * t->trr[j-1][D2M] // GD -> MM
	      + F_IM_prev[j-1] * q->trr[i-1][I2M] * t->trr[j-1][M2M] // IM -> MM
	      + F_DG_prev[j-1] * q->trr[i-1][D2M] * t->trr[j-1][M2M] // DG -> MM
	      + F_MI_prev[j-1] * q->trr[i-1][M2M] * t->trr[j-1][I2M] // MI -> MM
	    );
	  F_GD_curr[j] = 
	    (   F_MM_curr[j-1] * t->trr[j-1][M2D]                    // GD -> MM
	      + F_GD_curr[j-1] * t->trr[j-1][D2D]                    // GD -> GD
	    );
	  F_IM_curr[j] = 
	    (   F_MM_curr[j-1] * q->trr[i][M2I] * t->trr[j-1][M2M]     // MM -> IM
	      + F_IM_curr[j-1] * q->trr[i][I2I] * t->trr[j-1][M2M]     // IM -> IM
	    );
	  F_DG_curr[j] = scale[i] * 
	    (   F_MM_prev[j] * q->trr[i-1][M2D]                    // DG -> MM
	      + F_DG_prev[j] * q->trr[i-1][D2D]                    // DG -> DG
	    ) ;
	  F_MI_curr[j] = scale[i] * 
	    (   F_MM_prev[j] * q->trr[i-1][M2M] * t->trr[j][M2I]     // MI -> MM 
	      + F_MI_prev[j] * q->trr[i-1][M2M] * t->trr[j][I2I]     // MI -> MI
	    );

	  if(F_MM_curr[j]>Pmax_i)
	    Pmax_i=F_MM_curr[j];

	} // end else  	  

      } //end for j

      /* F_MM_prev = F_MM_curr */
      for (int jj = 0; jj <= t->L; jj++)
      {
	F_MM_prev[jj] = F_MM_curr[jj];
	F_MI_prev[jj] = F_MI_curr[jj];
	F_IM_prev[jj] = F_IM_curr[jj];
	F_DG_prev[jj] = F_DG_curr[jj];
	F_GD_prev[jj] = F_GD_curr[jj];

	/* and fill matrix because it is reused elsewhere */
	F_MM[i][jj] = F_MM_curr[jj];
      }

      pmin *= scale[i];
      if (pmin<DBL_MIN*100) pmin = 0.0;
      scale[i+1] = 1.0/(Pmax_i+1.0);
//    scale[i+1] = 1.0;   // to debug scaling

     
    } // end for i
  
  // Calculate P_forward * Product_{i=1}^{Lq+1}(scale[i])
  if (loc) 
    {
      Pforward = 1.0; // alignment contains no residues (see Mueckstein, Stadler et al.)
      for (i=1; i<=q->L; ++i) // Loop through query positions i
	{
	  jmin=1;
	  for (j=jmin; j<=t->L; ++j) // Loop through template positions j
	    Pforward += F_MM[i][j];
	  Pforward *= scale[i+1];
	}
    }
  else  // global alignment
    {
      Pforward = 0.0;
      for (i=1; i<q->L; ++i) Pforward = (Pforward + F_MM[i][t->L]) * scale[i+1];
      for (j=1; j<=t->L; ++j) Pforward += F_MM[q->L][j];
      Pforward *= scale[q->L+1];
    }

  // Calculate log2(P_forward)
  score = log2(Pforward)-10.0f;
  for (i=1; i<=q->L+1; ++i) score -= log2(scale[i]);
  //   state = MM;
  
  if (loc) 
    {
	score=score-log(t->L*q->L)/LAMDA+14.; // +14.0 to get approx same mean as for -global
    }

  return;
}



/////////////////////////////////////////////////////////////////////////////////////
// Compare two HMMs with Backward Algorithm (in lin-space, 2x faster), for use in MAC alignment 
/////////////////////////////////////////////////////////////////////////////////////
void Hit::Backward(HMM* q, HMM* t)
{

  // Variable declarations
  int i,j;      // query and template match state indices
  double pmin;  // this is the scaled 1 in the SW algorithm that represents a starting alignment         
  double Cshift = pow(2.0,shift);   // score offset transformed into factor in lin-space
  double scale_prod=scale[q->L+1];
  int jmin;
  
//  double B_MM_prev[MAXRES];
//  double B_DG_prev[MAXRES];
//  double B_MI_prev[MAXRES];

//  double B_MM_curr[MAXRES];
//  double B_GD_curr[MAXRES];
//  double B_DG_curr[MAXRES];
//  double B_IM_curr[MAXRES];
//  double B_MI_curr[MAXRES];

  // Initialization of top row, i.e. cells (0,j)
  for (int j=t->L; j>=1; j--) 
  {
    if (cell_off[q->L][j]) 
      B_MM[q->L][j] = B_MM_prev[j] = 0.0;
    else 
      B_MM[q->L][j] = B_MM_prev[j] = scale[q->L+1];
    B_MI_prev[j] = B_DG_prev[j] = 0.0;
  }
  if (loc) pmin = scale[q->L+1]; else pmin = 0.0; // transform pmin (for local alignment) to scale of present (i'th) row 

  // Backward algorithm
  for (i=q->L-1; i>=1; i--) // Loop through query positions i
    {
      jmin=1;

      // Initialize cells at (i,t->L+1)
      scale_prod *= scale[i+1];
      if (scale_prod<DBL_MIN*100) scale_prod = 0.0;
      if (cell_off[i][t->L]) 
	B_MM[i][t->L] = B_MM_curr[t->L] = 0.0;  
      else 
	B_MM[i][t->L] = B_MM_curr[t->L] = scale_prod;
      B_IM_curr[t->L] = B_MI_curr[t->L] = B_DG_curr[t->L] = B_GD_curr[t->L] = 0.0; 
      pmin *= scale[i+1]; // transform pmin (for local alignment) to scale of present (i'th) row 
      if (pmin<DBL_MIN*100) pmin = 0.0;
 
      for (j=t->L-1; j>=jmin; j--) // Loop through template positions j
	{
	  // Recursion relations
	  //	      printf("S[%i][%i]=%4.1f  ",i,j,Score(q->p[i],t->p[j]));

	  if (cell_off[i][j]) 
	    B_MM_curr[j] = B_GD_curr[j] = B_IM_curr[j] = B_DG_curr[j] = B_MI_curr[j] = 0.0;  
	  else 
	    {
	      double pmatch = B_MM_prev[j+1] * ProbFwd(q->p[i+1],t->p[j+1]) * Cshift * scale[i+1];
	      B_MM_curr[j] =  
		(
		 + pmin                                                    // MM -> EE (End/End, for local alignment)
		 + pmatch       * q->trr[i][M2M] * t->trr[j][M2M]              // MM -> MM
		 + B_GD_curr[j+1]                * t->trr[j][M2D]              // MM -> GD (q->tr[i][M2M] is already contained in GD->MM)
		 + B_IM_curr[j+1] * q->trr[i][M2I] * t->trr[j][M2M]              // MM -> IM
		 + B_DG_prev[j] * q->trr[i][M2D]                * scale[i+1] // MM -> DG (t->tr[j][M2M] is already contained in DG->MM)
		 + B_MI_prev[j] * q->trr[i][M2M] * t->trr[j][M2I] * scale[i+1] // MM -> MI
		 );
	      B_GD_curr[j] = 
		(
		 + pmatch       * q->trr[i][M2M] * t->trr[j][D2M]              // GD -> MM 
		 + B_GD_curr[j+1]                * t->trr[j][D2D]              // DG -> DG   
		 );
	      B_IM_curr[j] = 
		(
		 + pmatch       * q->trr[i][I2M] * t->trr[j][M2M]              // IM -> MM
		 + B_IM_curr[j+1] * q->trr[i][I2I] * t->trr[j][M2M]              // IM -> IM
		 );
	      B_DG_curr[j] =  
		(
		 + pmatch       * q->trr[i][D2M] * t->trr[j][M2M]              // DG -> MM
		 + B_DG_prev[j] * q->trr[i][D2D]                * scale[i+1] // DG -> DG
		 //   	         + B_GD[i][j+1] * q->tr[i][D2M] * t->tr[j][M2D]              // DG -> GD
		 );
	      B_MI_curr[j] = 
		(
		 + pmatch       * q->trr[i][M2M] * t->trr[j][I2M]              // MI -> MM       
		 + B_MI_prev[j] * q->trr[i][M2M] * t->trr[j][I2I] * scale[i+1] // MI -> MI
		 // 	         + B_IM[i][j+1] * q->tr[i][M2I] * t->tr[j][I2M]              // MI -> IM    
		 );

	    } // end else	      

	  /* Copy back to matrix */
	  B_MM[i][j] = B_MM_curr[j];
	} //end for j

      for(int jj = 0; jj <= t->L; jj++)
      {
	B_MM_prev[jj] = B_MM_curr[jj];
	B_DG_prev[jj] = B_DG_curr[jj];
	B_MI_prev[jj] = B_MI_curr[jj];
      }
    } // end for i
  
  // Calculate Posterior matrix and overwrite Backward matrix with it
  for (i=1; i<=q->L; ++i) 
    for (j=1; j<=t->L; ++j) 
      B_MM[i][j] *= F_MM[i][j]/Pforward;

  return;
}



/////////////////////////////////////////////////////////////////////////////////////
// Maximum Accuracy alignment 
/////////////////////////////////////////////////////////////////////////////////////
void Hit::MACAlignment(HMM* q, HMM* t)
{
  // Use Forward and Backward matrices to find that alignment which 
  // maximizes the expected number of correctly aligned pairs of residues (mact=0)
  // or, more generally, which maximizes the expectation value of the number of 
  // correctly aligned pairs minus (mact x number of aligned pairs)
  // "Correctly aligned" can be based on posterior probabilities calculated with
  // a local or a global version of the Forward-Backward algorithm.

  int i,j;           // query and template match state indices
  int jmin,jmax;     // range of dynamic programming for j
//  double S_prev[MAXRES];    // scores
//  double S_curr[MAXRES];    // scores
  double score_MAC;   // score of the best MAC alignment

  // Initialization of top row, i.e. cells (0,j)
  for (j=0; j<=t->L; ++j)
    S_prev[j] = 0.0;
  score_MAC=-INT_MAX; i2=j2=0; bMM[0][0]=STOP;

  // Dynamic programming 
  for (i=1; i<=q->L; ++i) // Loop through query positions i
    {

	{
	  // If q is compared to t, exclude regions where overlap of q with t < min_overlap residues
	  jmin=imax( 1, i+min_overlap-q->L);  // Lq-i+j>=Ovlap => j>=i+Ovlap-Lq => jmin=max{1, i+Ovlap-Lq} 
	  jmax=imin(t->L,i-min_overlap+t->L);  // Lt-j+i>=Ovlap => j<=i-Ovlap+Lt => jmax=min{Lt,i-Ovlap+Lt} 
	}      

      // Initialize cells
      S_curr[jmin-1] = 0.0;
      if (jmax<t->L) S_prev[jmax] = 0.0; // initialize at (i-1,jmax) if upper right triagonal is excluded due to min overlap
      
      for (j=jmin; j<=jmax; ++j) // Loop through template positions j
	{

	  if (cell_off[i][j]) 
	    {
	      S_curr[j] = -FLT_MIN;
	      bMM[i][j] = STOP;
	      //	      if (i>135 && i<140) 
	      // 		printf("Cell off   i=%i  j=%i b:%i\n",i,j,bMM[i][j]);
	    }
	  else 
	    {
	      // Recursion
	     
	      // NOT the state before the first MM state)
	      CALCULATE_MAX4(
			     S_curr[j],
			     B_MM[i][j] - mact,  // STOP signifies the first MM state, NOT the state before the first MM state (as in Viterbi)
			     S_prev[j-1] + B_MM[i][j] - mact, // B_MM[i][j] contains posterior probability
			     S_prev[j] - 0.5*mact,  // gap penalty prevents alignments such as this: XX--xxXX
			     S_curr[j-1] - 0.5*mact,  //                                               YYyy--YY  
			     bMM[i][j]   // backtracing matrix
			     );

	      //if (i>36 && i<40 && j>2200 && j<2230) 
	      //printf("i=%i  j=%i  S[i][j]=%8.3f  MM:%7.3f  MI:%7.3f  IM:%7.3f  b:%i\n",i,j,S[i][j],S[i-1][j-1]+B_MM[i][j]-mact,S[i-1][j],S[i][j-1],bMM[i][j]);
	      
	      // Find maximum score; global alignment: maximize only over last row and last column
	      if(S_curr[j]>score_MAC && (loc || i==q->L)) { i2=i; j2=j; score_MAC=S_curr[j]; }	      
	      
	    } // end if 
	  
	} //end for j
      
	  // if global alignment: look for best cell in last column
      if (!loc && S_curr[jmax]>score_MAC) { i2=i; j2=jmax; score_MAC=S_curr[jmax]; }
      
      for (j=0; j<=t->L; ++j)
	S_prev[j] = S_curr[j];
    } // end for i

  return;
}


/////////////////////////////////////////////////////////////////////////////////////
// Trace back Viterbi alignment of two profiles based on matrices bXX[][]
/////////////////////////////////////////////////////////////////////////////////////
void Hit::Backtrace(HMM* q, HMM* t)
{
  // Trace back trough the matrices bXY[i][j] until first match state is found (STOP-state)

  int step;      // counts steps in path through 5-layered dynamic programming matrix
  int i,j;       // query and template match state indices
  int wsval1;    // +1 for match, -1 for not-match
  int wsval2;

  L = t->L;
  
  // Make sure that backtracing stops when t:M1 or q:M1 is reached (Start state), e.g. sMM[i][1], or sIM[i][1] (M:MM, B:IM)
  for (i=0; i<=q->L; ++i) bMM[i][1]=bGD[i][1]=bIM[i][1] = STOP;
  for (j=1; j<=t->L; ++j) bMM[1][j]=bDG[1][j]=bMI[1][j] = STOP;
  

  // Back-tracing loop
  matched_cols=0; // for each MACTH (or STOP) state matched_col is incremented by 1
  step=0;         // steps through the matrix correspond to alignment columns (from 1 to nsteps)
  i=i2; j=j2;     // last aligned pair is (i2,j2)
  wsval1=1; wsval2=1;


  while (state)   // while (state!=STOP)  because STOP=0
    {
      step++;
      states[step] = state;
      wsi[step] = i;
      wsj[step] = j;
      // Exclude cells in direct neighbourhood from all further alignments
      for (int ii=imax(i-2,1); ii<=imin(i+2,q->L); ++ii)
	cell_off[ii][j]=1;     
      for (int jj=imax(j-2,1); jj<=imin(j+2,t->L); ++jj)
	cell_off[i][jj]=1;     
      
      switch (state)
	{
	case MM: // current state is MM, previous state is bMM[i][j]
	  matched_cols++; 
	  state = bMM[i--][j--];
	  wsval1=1;
	  wsval2=1;
	  break;	      
	case GD: // current state is GD
		wsval1=-1;wsval2=1;
	  switch (bGD[i][j--])
	    {
	    case STOP:  state = STOP; break; // current state does not have predecessor
	    case MM:    state = MM;   break; // previous state is Match state
	    }                               // default: previous state is same state (GD)
	  break;	      
	case IM: 
		wsval1=-1;wsval2=1;
	  switch (bIM[i][j--]) 
	    {
	    case STOP:  state = STOP; break; // current state does not have predecessor
	    case MM:    state = MM;   break; // previous state is Match state
	    }                               // default: previous state is same state (IM)
	  break;	      
	case DG:
		wsval1=1;wsval2=-1;
	  switch (bDG[i--][j])
	    {
	    case STOP:  state = STOP; break; // current state does not have predecessor
	    case MM:    state = MM;   break; // previous state is Match state
	    }                               // default: previous state is same state (DG)
	  break;	      
	case MI:
		wsval1=1;wsval2=-1;
	  switch (bMI[i--][j])
	    {
	    case STOP:  state = STOP; break; // current state does not have predecessor
	    case MM:    state = MM;   break; // previous state is Match state
	    }                               // default: previous state is same state (MI)
	  break;
	default:
	  fprintf(stderr,"Error: unallowed state value %i occurred during backtracing at step %i, (i,j)=(%i,%i)\n",state,step,i,j);
	  exit(-1);
	} //end switch (state)
//	wwi[step] = (i+1)*wsval1;
//	wwj[step] = (j+1)*wsval2;
    } //end while (state)
 
  i1 = abs(wsi[step]);
  j1 = abs(wsj[step]);
  states[step] = MM;  // first state (STOP state) is set to MM state
  nsteps=step; 


  // Add contribution from secondary structure score, record score along alignment,
  // and record template consensus sequence in master-slave-alignment to query sequence
  wsteps=0;
  for (step=1; step<=nsteps; step++)
    {
      switch(states[step])
	{
	case MM: 
	  i = wsi[step];
	  j = wsj[step];
	  S[step] = Score(q->p[i],t->p[j]);
	wwi[wsteps]=i;
	wwj[wsteps]=j;
	wsteps++;
	  break;
	case MI: //if gap in template  
	case DG:   
	default: //if gap in T or Q
	  S[step]=0.0f;
	  break;
	}
    }
 
  // Add contribution from correlation of neighboring columns to score
  float Scorr=0;
  if (nsteps)
    {
      for (step=2; step<=nsteps; step++) Scorr+=S[step]*S[step-1];
      for (step=3; step<=nsteps; step++) Scorr+=S[step]*S[step-2];
      for (step=4; step<=nsteps; step++) Scorr+=S[step]*S[step-3];
      for (step=5; step<=nsteps; step++) Scorr+=S[step]*S[step-4];
      score+=corr*Scorr;
    }
  
  // Set score, P-value etc.
  score_sort = score_aass = -score;

  return;
}


/////////////////////////////////////////////////////////////////////////////////////
// Trace back alignment of two profiles based on matrices bXX[][]
/////////////////////////////////////////////////////////////////////////////////////
void Hit::BacktraceMAC(HMM* q, HMM* t)
{
  // Trace back trough the matrix b[i][j] until STOP state is found

  char** b=bMM;  // define alias for backtracing matrix
  int step;      // counts steps in path through 5-layered dynamic programming matrix
  int i,j;       // query and template match state indices
  int wsval1;    // +1 for match, -1 for not-match
  int wsval2;

  L = t->L;
  
  // Make sure that backtracing stops when t:M1 or q:M1 is reached (Start state), e.g. sMM[i][1], or sIM[i][1] (M:MM, B:IM)
  for (i=0; i<=q->L; ++i) b[i][1] = STOP;
  for (j=1; j<=t->L; ++j) b[1][j] = STOP;
  
  // Back-tracing loop
  // In contrast to the Viterbi-Backtracing, STOP signifies the first Match-Match state, NOT the state before the first MM state
  matched_cols=1; // for each MACTH (or STOP) state matched_col is incremented by 1
  state=MM;       // lowest state with maximum score must be match-match state 
  step=0;         // steps through the matrix correspond to alignment columns (from 1 to nsteps)
  i=i2; j=j2;     // last aligned pair is (i2,j2)


  wsval1=1; wsval2=1;
  if (b[i][j] != MM) {
    step = 1;
    wsi[step] = i;
    wsj[step] = j;
    state = STOP;
  } else {
    while (state!=STOP) 
      {
	step++;
	states[step] = state = b[i][j];
	wsi[step] = i;
	wsj[step] = j;
//	wwi[step] = i*wsval1;
//	wwj[step] = j*wsval2;
	// Exclude cells in direct neighbourhood from all further alignments
	for (int ii=imax(i-2,1); ii<=imin(i+2,q->L); ++ii)
	  cell_off[ii][j]=1;     
	for (int jj=imax(j-2,1); jj<=imin(j+2,t->L); ++jj)
	  cell_off[i][jj]=1;     
	if (state==MM) matched_cols++; 
	
	switch (state)
	  {
	  case MM: i--; j--; wsval1=1;wsval2=1; break;
	  case IM: j--; wsval1=-1;wsval2=1; break;
	  case MI: i--; wsval1=1;wsval2=-1; break;
	  case STOP: break;
	  default:
	    fprintf(stderr,"Errorin: unallowed state value %i occurred during backtracing at step %i, (i,j)=(%i,%i)\n",state,step,i,j);
	    exit(-1);
	  } //end switch (state)
//	wwi[step] = (i+1)*wsval1;
//	wwj[step] = (j+1)*wsval2;
      } //end while (state)
  }  
  i1 = abs(wsi[step]);
  j1 = abs(wsj[step]);
  states[step] = MM;  // first state (STOP state) is set to MM state
  nsteps=step; 


  // Add contribution from secondary structure score, record score along alignment,
  // and record template consensus sequence in master-slave-alignment to query sequence
  sum_of_probs=0.0;       // number of identical residues in query and template sequence
  //   printf("Hit=%s\n",name); /////////////////////////////////////////////////////////////
  wsteps=0;
  for (step=1; step<=nsteps; step++)
    {
      switch(states[step])
	{
	case MM: 
	  i = wsi[step];
	  j = wsj[step];
	wwi[wsteps]=i;
	wwj[wsteps]=j;
	wsteps++;
	  S[step] = Score(q->p[i],t->p[j]);
	  P_posterior[step] = B_MM[wsi[step]][wsj[step]];
	  // Add probability to sum of probs if no dssp states given or dssp states exist and state is resolved in 3D structure
	  sum_of_probs += P_posterior[step]; 
	  // 	  printf("j=%-3i  dssp=%1i  P=%4.2f  sum=%6.2f\n",j,t->ss_dssp[j],P_posterior[step],sum_of_probs); //////////////////////////
	  break;
	case MI: //if gap in template  
	case DG:   
	default: //if gap in T or Q
	  S[step] = P_posterior[step] = 0.0;
	  break;
	}
    }

  // Add contribution from correlation of neighboring columns to score
  float Scorr=0;
  if (nsteps)
    {
      for (step=1; step<=nsteps-1; step++) Scorr+=S[step]*S[step+1];
      for (step=1; step<=nsteps-2; step++) Scorr+=S[step]*S[step+2];
      for (step=1; step<=nsteps-3; step++) Scorr+=S[step]*S[step+3];
      for (step=1; step<=nsteps-4; step++) Scorr+=S[step]*S[step+4];
      score+=corr*Scorr;
    }
  
  // Set score, P-value etc.
  score_sort = score_aass = -score;

  return;
}



/////////////////////////////////////////////////////////////////////////////////////
//// Functions that calculate probabilities
/////////////////////////////////////////////////////////////////////////////////////
void Hit::InitializeForAlignment(HMM* q, HMM* t)
{
  int i,j;

  // SS scoring during (ssm2>0) or after (ssm1>0) alignment? Query SS known or Template SS known?
  int ssm=0;
  switch (ssm) 
    {
    case 0:
      break;
    }


    {
	// Compare two different HMMs Q and T
	  // Activate all cells in dynamic programming matrix
	  for (i=1; i<=q->L; ++i) 
	    for (j=1; j<=t->L; ++j) 
	      cell_off[i][j]=0;   // no other cells crossed out yet
	    
	  // Cross out cells that are excluded by the minimum-overlap criterion
	    min_overlap = imin(60, (int)(0.333f*imin(q->L,t->L))+1); // automatic minimum overlap
	  
	  for (i=0; i<min_overlap; ++i) 
	    for (j=i-min_overlap+t->L+1; j<=t->L; ++j) // Lt-j+i>=Ovlap => j<=i-Ovlap+Lt => jmax=min{Lt,i-Ovlap+Lt} 
	      cell_off[i][j]=1;
	  for (i=q->L-min_overlap+1; i<=q->L; ++i) 
	    for (j=1; j<i+min_overlap-q->L; ++j)      // Lq-i+j>=Ovlap => j>=i+Ovlap-Lq => jmin=max{1, i+Ovlap-Lq} 
	      cell_off[i][j]=1;
    } 
}
	

//Calculate score between columns i and j of two HMMs (query and template)
inline float Hit::Score(float* qi, float* tj)
{
  return fast_log2(ProbFwd(qi,tj));
}

// Calculate score between columns i and j of two HMMs (query and template)
inline float Hit::ProbFwd(float* qi, float* tj)
{
	return fast_dot_product_single2(qi,tj,20);
//	return ScalarProd20(qi,tj); //
}


/////////////////////////////////////////////////////////////////////////////////////
//// Function for Viterbi()
/////////////////////////////////////////////////////////////////////////////////////
inline float max2(const float& xMM, const float& xX, char& b) 
{
  if (xMM>xX) { b=MM; return xMM;} else { b=SAME;  return xX;}
}


/////////////////////////////////////////////////////////////////////////////////////
// Functions for StochasticBacktrace()
/////////////////////////////////////////////////////////////////////////////////////

inline int pickprob2(const double& xMM, const double& xX, const int& state) 
{
  if ( (xMM+xX)*frand() < xMM) return MM; else return state; 
}
inline int pickprob3_GD(const double& xMM, const double& xDG, const double& xGD) 
{
  double x = (xMM+xDG+xGD)*frand();
  if ( x<xMM) return MM; 
  else if ( x<xMM+xDG) return DG; 
  else return GD;
}
inline int pickprob3_IM(const double& xMM, const double& xMI, const double& xIM) 
{
  double x = (xMM+xMI+xIM)*frand();
  if ( x<xMM) return MM; 
  else if ( x<xMM+xMI) return MI; 
  else return IM;
}
inline int pickprob6(const double& x0, const double& xMM, const double& xGD, const double& xIM, const double& xDG, const double& xMI) 
{
  double x = (x0+xMM+xGD+xIM+xDG+xMI)*frand();
  x-=xMM; if (x<0) return MM; 
  x-=x0;  if (x<0) return STOP; 
  x-=xGD; if (x<0) return GD;
  x-=xIM; if (x<0) return IM;
  if (x < xDG) return DG; else return MI;
}

inline int pickmax2(const double& xMM, const double& xX, const int& state) 
{
  if (xMM > xX) return MM; else return state; 
}
inline int pickmax3_GD(const double& xMM, const double& xDG, const double& xGD) 
{
  char state;
  double x;
  if ( xMM>xDG) {state=MM; x=xMM;} 
  else          {state=DG; x=xDG;}
  if ( xGD>x)   {state=GD; x=xGD;}
  return state;
}inline int pickmax3_IM(const double& xMM, const double& xMI, const double& xIM) 
{
  char state;
  double x;
  if ( xMM>xMI) {state=MM; x=xMM;}
  else          {state=MI; x=xMI;}
  if ( xIM>x)   {state=IM; x=xIM;}
  return state;
}
inline int pickmax6(const double& x0, const double& xMM, const double& xGD, const double& xIM, const double& xDG, const double& xMI) 
{
  char state;
  double x;
  if ( x0 >xMM) {state=STOP; x=x0;} 
  else          {state=MM; x=xMM;}
  if ( xGD>x)   {state=GD; x=xGD;}
  if ( xIM>x)   {state=IM; x=xIM;}
  if ( xDG>x)   {state=DG; x=xDG;}
  if ( xMI>x)   {state=MI; x=xMI;}
  return state;
}



/////////////////////////////////////////////////////////////////////////////////////
//// Functions that calculate P-values and probabilities 
/////////////////////////////////////////////////////////////////////////////////////


//// Evaluate the CUMULATIVE extreme value distribution at point x
//// p(s)ds = lamda * exp{ -exp[-lamda*(s-mu)] - lamda*(s-mu) } ds = exp( -exp(-x) - x) dx = p(x) dx
//// => P(s>S) = integral_-inf^inf {p(x) dx}  = 1 - exp{ -exp[-lamda*(S-mu)] }
inline double Pvalue(double x, double a[])
{
  //a[0]=lamda, a[1]=mu
  double h = a[0]*(x-a[1]);
  return (h>10)? exp(-h) : double(1.0)-exp( -exp(-h));
}

inline double Pvalue(float x, float lamda, float mu)
{
  double h = lamda*(x-mu);
  return (h>10)? exp(-h) : (double(1.0)-exp( -exp(-h)));
}

inline double logPvalue(float x, float lamda, float mu)
{
  double h = lamda*(x-mu);
  return (h>10)? -h : (h<-2.5)? -exp(-exp(-h)): log( ( double(1.0) - exp(-exp(-h)) ) );
}

inline double logPvalue(float x, double a[])
{
  double h = a[0]*(x-a[1]);
  return (h>10)? -h : (h<-2.5)? -exp(-exp(-h)): log( ( double(1.0) - exp(-exp(-h)) ) );
}

// Calculate probability of true positive : p_TP(score)/( p_TP(score)+p_FP(score) )
// TP: same superfamily OR MAXSUB score >=0.1
inline double Probab(Hit& hit)
{
  double s=-hit.score_aass;
  double t;
  if (s>200) return 100.0; 
  if (hit.loc) 
    {
	{
	  // local no SS
	  const double a=sqrt(4000.0);
	  const double b=2.0*2.5;
	  const double c=sqrt(0.15);
	  const double d=2.0*34.0;
	  t = a*exp(-s/b) + c*exp(-s/d);
	}
    }
  else
    {
	{
	  // global no SS
	  const double a=sqrt(6000.0);
	  const double b=2.0*2.5;
	  const double c=sqrt(0.10);
	  const double d=2.0*37.0;
	  t = a*exp(-s/b) + c*exp(-s/d);
	}
    }

  return 100.0/(1.0+t*t);
}

