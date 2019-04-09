#pragma once
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <cstdio>     // printf
#include <stdlib.h>   // exit
#include <time.h>     // clock
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <string.h>   // strcmp, strstr
#include <cassert>
#ifdef HH_SSE3
#include <emmintrin.h>
#include <pmmintrin.h>
#endif
using namespace std;


#define NAA 20
#define ANY 20
#define GAP 21
#define ENDGAP 22
#define NTRANS 7
#define MAXSEQ 65535     // max number of sequences in input alignment (must be <~30000 on cluster nodes??)
#define MAXRES 15002     // max number of states in HMM; must be <= LINELEN
#define MAXCOL 32765     // max number of columns in sequence/MSA input files; must be <= LINELEN and >= maxres
#define LINELEN 524288   // max length of line read in from input files; must be >= MAXCOL
#define DESCLEN 32765    // max length of sequence description (longname)
#define NAMELEN 512      // max length of file names etc., defined in limits.h  ( EXTERN const int NAMELEN=(PATH_MAX>512? PATH_MAX:512); )
#define HMMSCALE 1000 
#define SELFEXCL 3
#define LAMDA 0.388

//------- transition utility -------//
enum transitions {M2M,M2I,M2D,I2M,I2I,D2M,D2D};
enum pair_states {STOP=0,SAME=1,GD=2,IM=3,DG=4,MI=5,MS=6,ML=7,SM=8,LM=9,MM=10};




/////////////////////////////////////////////////////////////////////////////////////
// Transforms the one-letter amino acid code into an integer between 0 and 22
/////////////////////////////////////////////////////////////////////////////////////
inline char aa2i(char c)
{
  //A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
  if (c>='a' && c<='z') c+='A'-'a';
  switch (c)
    {
    case 'A': return 0;
    case 'R': return 1;
    case 'N': return 2;
    case 'D': return 3;
    case 'C': return 4;
    case 'Q': return 5;
    case 'E': return 6;
    case 'G': return 7;
    case 'H': return 8;
    case 'I': return 9;
    case 'L': return 10;
    case 'K': return 11;
    case 'M': return 12;
    case 'F': return 13;
    case 'P': return 14;
    case 'S': return 15;
    case 'T': return 16;
    case 'W': return 17;
    case 'Y': return 18;
    case 'V': return 19;
    case 'X': return ANY;
    case 'J': return ANY;
    case 'O': return ANY;
    case 'U': return 4;  //Selenocystein -> Cystein
    case 'B': return 3;  //D (or N)
    case 'Z': return 6;  //E (or Q)
    case '-': return GAP;
    case '.': return GAP;
    case '_': return GAP;
    }
  if (c>=0 && c<=32) return -1; // white space and control characters
  return -2;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms integers between 0 and 22 into the one-letter amino acid code
/////////////////////////////////////////////////////////////////////////////////////
inline char i2aa(char c)
{
  //A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
  switch (c)
    {
    case 0: return 'A';
    case 1: return 'R';
    case 2: return 'N';
    case 3: return 'D';
    case 4: return 'C';
    case 5: return 'Q';
    case 6: return 'E';
    case 7: return 'G';
    case 8: return 'H';
    case 9: return 'I';
    case 10: return 'L';
    case 11: return 'K';
    case 12: return 'M';
    case 13: return 'F';
    case 14: return 'P';
    case 15: return 'S';
    case 16: return 'T';
    case 17: return 'W';
    case 18: return 'Y';
    case 19: return 'V';
    case ANY: return 'X';
    case GAP: return '-';
    case ENDGAP: return '-';
    }
  return '?';
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms the dssp/psipred secondary structure code into an integer number
/////////////////////////////////////////////////////////////////////////////////////
inline char ss2i(char c)
{
  //- H E C S T G B
  if (c>='a' && c<='z') c+='A'-'a';
  switch (c)
    {
    case '.': return 0;
    case '-': return 0;
    case 'X': return 0;
    case 'H': return 1;
    case 'E': return 2;
    case 'C': return 3;
    case '~': return 3;
    case 'S': return 4;
    case 'T': return 5;
    case 'G': return 6;
    case 'B': return 7;
    case 'I': return 3;
    case ' ': return -1;
    case '\t': return -1;
    case '\n': return -1;
    }
  return -2;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms integers between 0 and 8 into the dssp/psipred secondary structure code
/////////////////////////////////////////////////////////////////////////////////////
inline char i2ss(int c)
{
  //- H E C S T G B
  switch (c)
    {
    case 0: return '-';
    case 1: return 'H';
    case 2: return 'E';
    case 3: return 'C';
    case 4: return 'S';
    case 5: return 'T';
    case 6: return 'G';
    case 7: return 'B';
    case 8: return 'I';
    }
  return '?';
}


/////////////////////////////////////////////////////////////////////////////////////
// Transforms the solvend accessiblity code into an integer number
/////////////////////////////////////////////////////////////////////////////////////
inline char sa2i(char c)
{
  //- A B C D E
  if (c>='a' && c<='z') c+='A'-'a';
  switch (c)
    {
    case '.': return 0;
    case '-': return 0;
    case 'A': return 1;
    case 'B': return 2;
    case 'C': return 3;
    case 'D': return 4;
    case 'E': return 5;
    case 'F': return 6;
    case ' ': return -1;
    case '\t': return -1;
    case '\n': return -1;
    }
  return -2;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms integers between 0 and 5 into the solvent accessibility code
/////////////////////////////////////////////////////////////////////////////////////
inline char i2sa(int c)
{
  //- H E C S T G B
  switch (c)
    {
    case 0: return '-';
    case 1: return 'A';
    case 2: return 'B';
    case 3: return 'C';
    case 4: return 'D';
    case 5: return 'E';
    case 6: return 'F';
    }
  return '?';
}


/////////////////////////////////////////////////////////////////////////////////////
// Transforms alternative secondary structure symbols into symbols
/////////////////////////////////////////////////////////////////////////////////////
inline char ss2ss(char c)
{
  //- H E C S T G B
  switch (c)
    {
    case '~': return 'C';
    case 'I': return 'C';
    case 'i': return 'c';
    case 'H':
    case 'E':
    case 'C':
    case 'S':
    case 'T':
    case 'G':
    case 'B':
    case 'h':
    case 'e':
    case 'c':
    case 's':
    case 't':
    case 'g':
    case 'b':
    case '.':
      return c;
    }
  return '-';
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms confidence values of psipred into internal code
/////////////////////////////////////////////////////////////////////////////////////
inline char cf2i(char c)
{
  switch (c)
    {
    case '-': return 0;
    case '.': return 0;
    case '0': return 1;
    case '1': return 2;
    case '2': return 3;
    case '3': return 4;
    case '4': return 5;
    case '5': return 6;
    case '6': return 7;
    case '7': return 8;
    case '8': return 9;
    case '9': return 10;
    }
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms internal representation of psipred confidence values into printable chars
/////////////////////////////////////////////////////////////////////////////////////
inline char i2cf(char c)
{
  switch (c)
    {
    case 0: return '-';
    case 1: return '0';
    case 2: return '1';
    case 3: return '2';
    case 4: return '3';
    case 5: return '4';
    case 6: return '5';
    case 7: return '6';
    case 8: return '7';
    case 9: return '8';
    case 10: return '9';
    }
  return '-';
}


//------- output --------//
inline void fout(FILE* outf, int d)
{
  if (d>=99999) fprintf(outf,"*\t"); else fprintf(outf,"%i\t",d);
  return;
}



/////////////////////////////////////////////////////////////////////////////////////
// Transform a character to lower case and '.' to '-' and vice versa
/////////////////////////////////////////////////////////////////////////////////////
inline char MatchChr(char c)  {return ((c>='a' && c<='z')? c-'a'+'A' : (c=='.'? '-':c) );}
inline char InsertChr(char c) {return ((c>='A' && c<='Z')? c+'a'-'A' : ((c>='0' && c<='9') || c=='-')? '.':c );}
inline int  WordChr(char c) {return (int)((c>='A' && c<='Z') || (c>='a' && c<='z'));}



/////////////////////////////////////////////////////////////////////////////////////
// Arithmetics
/////////////////////////////////////////////////////////////////////////////////////

//// max and min
inline double dmax(double x, double y) { return (x>y? x : y);}
inline double dmin(double x, double y) { return (x<y? x : y);}
inline int imax(int x, int y) { return (x>y? x : y);}
inline int imin(int x, int y) { return (x<y? x : y);}
inline int iabs(int x) { return (x>=0? x : -x);}

// Rounding up, rounding down and rounding to nearest integer
inline int iceil(double x)  {return int(ceil(x));}
inline int ifloor(double x) {return int(floor(x));}
inline int iround(double x) {return int(floor(x+0.5));}

//// Generalized mean: d=0: sqrt(x*y)  d=1: (x+y)/2  d->-inf: min(x,y)  d->+inf: max(x,y)
inline double fmean(double x, double y, double d) { return pow( (pow(x,d)+pow(y,d))/2 ,1./d);}
inline float frand() { return rand()/(RAND_MAX+1.0); }

// log base 2
inline float log2(float x)  {return (x<=0? (float)(-100000):1.442695041*log(x));}
inline float log10(float x) {return (x<=0? (float)(-100000):0.434294481*log(x));}

//------- fast function -------//
extern float flog2(float x);
extern float fast_log2(float x);
extern float fast_log_gamma(float x);
extern float fpow2(float x);

//------- normalization -------//
extern double normalize_to_one(double* array, size_t length, const float* def_array=NULL);
extern float normalize_to_one(float* array, size_t length, const float* def_array=NULL);
extern float NormalizeTo1(float* array, int length, float* def_array=NULL);
extern float NormalizeToX(float* array, int length, float x, float* def_array=NULL);
extern char* sprintg(float val, int w);

//------- string function -----//
extern int strmcpy(char* dest, const char* source, size_t maxlen);
extern int strtoi(const char*& ptr);
extern int strtoi_(const char*& ptr, int deflt=INT_MAX);
extern int strint(char*& ptr);
extern int strinta(char*& ptr, int deflt=99999);
extern float strflt(char*& ptr);
extern float strflta(char*& ptr, float deflt=99999);
extern int chomp(char str[]);
extern char* fgetline(char str[], const int maxlen, FILE* file);
extern char *substr(char* substr, char* str, int a, int b);
extern char* strscn(char* str);
extern char* strscn_ws(char* str);
extern const char* strscn_c(const char* str);
extern char* strscn_(char* str);
extern char* strscn(char* str, const char c);
extern char* strscn_(char* str, const char c);
extern char* strcut(char* str);
extern char* strcut_(char* str);
extern char* strcut(char* str, const char c);
extern char* strcut_(char* str, const char c);
extern char* strcut(char* str, const char* substr);
extern char* strcut_(char* str, const char* substr);
extern char* strwrd(char* str, char* ptr);
extern char* strwrd(char* str, char* ptr, int maxlen);
extern char* strwrd(char* str, char* ptr, const char c);
extern int strtr(char* str, const char oldchars[], const char newchars[]);
extern int strtrd(char* str, const char chars[]);
extern int strtrd(char* str, char char1, char char2);
extern int strcount(char* str, char char1, char char2);
extern char* uprstr(char* str);
extern char* lwrstr(char* str);
extern char uprchr(char chr);
extern char lwrchr(char chr);
extern char* strsubst(char* str, const char str1[], const char str2[]);

//------ path related ------//
//extern void ElapsedTimeSinceFirstCall(const char str[]);
//extern void ElapsedTimeSinceLastCall(const char str[]);
extern char* RemovePath(char outname[], char filename[]);
extern char* RemoveExtension(char outname[], char filename[]);
extern char* RemovePathAndExtension(char outname[], char filename[]);
extern char* Extension(char extension[], char filename[]);
extern char* Pathname(char pathname[], char filename[]);

//------- qsort -------//
extern void swapi(int k[], int i, int j);
extern void QSortInt(int v[], int k[], int left, int right, int up=+1);
extern void QSortFloat(float v[], int k[], int left, int right, int up=+1);
extern void QSortDouble(double v[], int k[], int left, int right, int up=+1);

//------ others -----//
extern float fast_dot_product_single2(float* qi, float* tj,int num);
extern float ScalarProd20(float* qi, float* tj);


//------ errors -----//
extern int FormatError(const char infile[], const char details[]="");
extern int OpenFileError(const char outfile[]);
extern int MemoryError(const char arrayname[]);
extern int NoMemoryError(const char arrayname[]);
extern int SyntaxError(const char details[]="");
extern int InternalError(const char errstr[]);


