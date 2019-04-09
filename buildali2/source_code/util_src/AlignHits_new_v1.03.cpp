#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <getopt.h>
using namespace std;



//=================== upper and lower case ====================//
//----------upper_case-----------//
void toUpperCase(char *buffer)
{
	for(int i=0;i<(int)strlen(buffer);i++)
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
void toUpperCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++)
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
//----------lower_case-----------//
void toLowerCase(char *buffer)
{
	for(int i=0;i<(int)strlen(buffer);i++)
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}
void toLowerCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++)
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}

//----- get upper case -----//
int getUpperCase(char *buffer)
{
	int count=0;
	for(int i=0;i<(int)strlen(buffer);i++)
	if(buffer[i]>=65 && buffer[i]<=90) count++;
	return count;
}
int getUpperCase(string &buffer)
{
	int count=0;
	for(int i=0;i<(int)buffer.length();i++)
	if(buffer[i]>=65 && buffer[i]<=90) count++;
	return count;
}
//----- get lower case -----//
int getLowerCase(char *buffer)
{
	int count=0;
	for(int i=0;i<(int)strlen(buffer);i++)
	if(buffer[i]>=97 && buffer[i]<=122) count++;
	return count;
}
int getLowerCase(string &buffer)
{
	int count=0;
	for(int i=0;i<(int)buffer.length();i++)
	if(buffer[i]>=97 && buffer[i]<=122) count++;
	return count;
}

//--------- base_name -----------//__110830__//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}



//=============== public data ================//
// internal         A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
const int ch2psi[]={1, 4, 3, 4, 5, 6, 7, 8, 9,21,10,11,12,13,21,14,15,16,17,18, 3,19,20,21,22, 5};
//                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 
// psi              ?  A  B  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  X  Y  Z  ?  ?
//                  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
const int ch2i[]=  {0, 3, 4, 3, 6,13, 7, 8, 9,20,11,10,12, 2,20,14, 5, 1,15,16, 4,19,17,20,18, 6};


// The Gonnet matrix is in units of 10*log10()
// 0    1    2    3    4    5    6    7    8    9   10   11    12  13   14   15   16   17   18   19   20 
// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V    X    
const double GONNET[] = {
 2.4,-0.6,-0.3,-0.3, 0.5,-0.2, 0.0, 0.5,-0.8,-0.8,-1.2,-0.4,-0.7,-2.3, 0.3, 1.1, 0.6,-3.6,-2.2, 0.1, 0.0, // A
-0.6, 4.7, 0.3,-0.3,-2.2, 1.5, 0.4,-1.0, 0.6,-2.4,-2.2, 2.7,-1.7,-3.2,-0.9,-0.2,-0.2,-1.6,-1.8,-2.0, 0.0, // R
-0.3, 0.3, 3.8, 2.2,-1.8, 0.7, 0.9, 0.4, 1.2,-2.8,-3.0, 0.8,-2.2,-3.1,-0.9, 0.9, 0.5,-3.6,-1.4,-2.2, 0.0, // N
-0.3,-0.3, 2.2, 4.7,-3.2, 0.9, 2.7, 0.1, 0.4,-3.8,-4.0, 0.5,-3.0,-4.5,-0.7, 0.5, 0.0,-5.2,-2.8,-2.9, 0.0, // D
 0.5,-2.2,-1.8,-3.2,11.5,-2.4,-3.0,-2.0,-1.3,-1.1,-1.5,-2.8,-0.9,-0.8,-3.1, 0.1,-0.5,-1.0,-0.5, 0.0, 0.0, // C
-0.2, 1.5, 0.7, 0.9,-2.4, 2.7, 1.7,-1.0, 1.2,-1.9,-1.6, 1.5,-1.0,-2.6,-0.2, 0.2, 0.0,-2.7,-1.7,-1.5, 0.0, // Q
 0.0, 0.4, 0.9, 2.7,-3.0, 1.7, 3.6,-0.8, 0.4,-2.7,-2.8, 1.2,-2.0,-3.9,-0.5, 0.2,-0.1,-4.3,-2.7,-1.9, 0.0, // E
 0.5,-1.0, 0.4, 0.1,-2.0,-1.0,-0.8, 6.6,-1.4,-4.5,-4.4,-1.1,-3.5,-5.2,-1.6, 0.4,-1.1,-4.0,-4.0,-3.3, 0.0, // G
-0.8, 0.6, 1.2, 0.4,-1.3, 1.2, 0.4,-1.4, 6.0,-2.2,-1.9, 0.6,-1.3,-0.1,-1.1,-0.2,-0.3,-0.8,-2.2,-2.0, 0.0, // H
-0.8,-2.4,-2.8,-3.8,-1.1,-1.9,-2.7,-4.5,-2.2, 4.0, 2.8,-2.1, 2.5, 1.0,-2.6,-1.8,-0.6,-1.8,-0.7, 3.1, 0.0, // I
-1.2,-2.2,-3.0,-4.0,-1.5,-1.6,-2.8,-4.4,-1.9, 2.8, 4.0,-2.1, 2.8, 2.0,-2.3,-2.1,-1.3,-0.7, 0.0, 1.8, 0.0, // L
-0.4, 2.7, 0.8, 0.5,-2.8, 1.5, 1.2,-1.1, 0.6,-2.1,-2.1, 3.2,-1.4,-3.3,-0.6, 0.1, 0.1,-3.5,-2.1,-1.7, 0.0, // K
-0.7,-1.7,-2.2,-3.0,-0.9,-1.0,-2.0,-3.5,-1.3, 2.5, 2.8,-1.4, 4.3, 1.6,-2.4,-1.4,-0.6,-1.0,-0.2, 1.6, 0.0, // M
-2.3,-3.2,-3.1,-4.5,-0.8,-2.6,-3.9,-5.2,-0.1, 1.0, 2.0,-3.3, 1.6, 7.0,-3.8,-2.8,-2.2, 3.6, 5.1, 0.1, 0.0, // F
 0.3,-0.9,-0.9,-0.7,-3.1,-0.2,-0.5,-1.6,-1.1,-2.6,-2.3,-0.6,-2.4,-3.8, 7.6, 0.4, 0.1,-5.0,-3.1,-1.8, 0.0, // P
 1.1,-0.2, 0.9, 0.5, 0.1, 0.2, 0.2, 0.4,-0.2,-1.8,-2.1, 0.1,-1.4,-2.8, 0.4, 2.2, 1.5,-3.3,-1.9,-1.0, 0.0, // S
 0.6,-0.2, 0.5, 0.0,-0.5, 0.0,-0.1,-1.1,-0.3,-0.6,-1.3, 0.1,-0.6,-2.2, 0.1, 1.5, 2.5,-3.5,-1.9, 0.0, 0.0, // T
-3.6,-1.6,-3.6,-5.2,-1.0,-2.7,-4.3,-4.0,-0.8,-1.8,-0.7,-3.5,-1.0, 3.6,-5.0,-3.3,-3.5,14.2, 4.1,-2.6, 0.0, // W
-2.2,-1.8,-1.4,-2.8,-0.5,-1.7,-2.7,-4.0,-2.2,-0.7, 0.0,-2.1,-0.2, 5.1,-3.1,-1.9,-1.9, 4.1, 7.8,-1.1, 0.0, // Y
 0.1,-2.0,-2.2,-2.9, 0.0,-1.5,-1.9,-3.3,-2.0, 3.1, 1.8,-1.7, 1.6, 0.1,-1.8,-1.0, 0.0,-2.6,-1.1, 3.4, 0.0, // V
 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  // X	  
};

// The BLOSUM62 matrix
// A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X
const int BLOSUM62[] = {
  4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0, 0,
 -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1,
 -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,-1,
 -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,-1,
  0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-2,
 -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,-1,
 -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,-1,
  0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,
 -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,-1,
 -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-1,
 -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-1,
 -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,-1,
 -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-1,
 -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-1,
 -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,
  1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0,
  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0, 0,
 -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-2,
 -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-1,
  0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-1,
  0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1
};

//---- public data structure ----//
double *GONNET_CUR; //-> 21*21
double *BLOSUM_CUR; //-> 21*21
double *PSSM_MTX;   //-> L*26 where L is canonical query length
int *QUERY_MATCH;  //-> length is L. 0 for not match; 1 for match
int QUERY_LENGTH;  //-> L, canonical query length
string QUERY_SEQ;  //-> query sequence from file, might contain lower case.
                   //-> the length might not be L
vector <string> INPUT_ALI; //-> input alignment from a3m. upper case for match, lower case for insert
                           //-> the first sequence is canonical query, all in upper case.
vector <string> INPUT_NAM; //-> input name from a3m
vector <int> INPUT_LEN;    //-> input sequence's length
vector <int> PRUNE_LIST;   //-> final prune list. 0 for remove; 1 for retain
double GAP_PENALTY=11.0/3.0;        //-> gap opening penalty in bits (for BLOSUM62: 11 bits/3)
double EXTEND_PENALTY=1.0/3.0;      //-> gap extension penalty in bits (for BLOSUM62: 1 bits/3)



//---------- data structure init ---------//
void Data_Structure_Init(int input_len)
{
	int AAnum=20;
	GONNET_CUR=new double[(AAnum+1)*(AAnum+1)];
	BLOSUM_CUR=new double[(AAnum+1)*(AAnum+1)];
	PSSM_MTX=new double[input_len*26];
	QUERY_MATCH=new int[input_len];
	//rescale gonnet
	for(int i=0;i<=20;i++)for(int j=0;j<=20;j++)GONNET_CUR[i*21+j]=0.3322*GONNET[i*21+j];
	//init QUERY_MATCH
	for(int i=9;i<input_len;i++)QUERY_MATCH[i]=1; //default is all match
}

//------ read fasta sequence -----//
int Read_FASTA_SEQRES(string &infile,string &seqres,int skip=1) //->from .fasta file
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such file! %s \n",infile.c_str());
		return -1;
	}
	//skip
	int i;
	for(i=0;i<skip;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file format bad! %s \n",infile.c_str());
			return -1;
		}
	}
	//process
	temp="";
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		temp+=buf;
	}
	seqres=temp;
	//return
	return (int)seqres.length();
}

//---------- dynamic programming ----------//
int WWW_Advance_Align_Dyna_Prog_Double(int n1,int n2,const vector<double> &score,
								   double GAP_OPEN1,double GAP_EXT1,double GAP_OPEN2,double GAP_EXT2,
								   double GAP_HEAD1,double GAP_TAIL1,double GAP_HEAD2,double GAP_TAIL2,
								   vector<pair<int,int> > & alignment,double &ali_sco)
{
	int i,j;
	//input
	int m = n1 + 1;  // +1 to account for the extra row,col in
	int n = n2 + 1;  // the DP matrices corresponding to gaps
	int DP_maximal=n;
	int IN_maximal=n2;
	//const value
	const int _H_  = 0;
	const int _S_  = 1;
	const int _V_  = 2;

	//create D and M
	vector <int> D[3];      // the path (directions) matrix
	vector <double> M[3];   // the current scores (values) matrix
	//resize(m,n)
	for (i = 0; i < 3; ++i) 
	{
		D[i].resize(m*n);
		M[i].resize(m*n);
	}
	//init()
	double WS_MIN=-1000000;
	D[_S_][0*DP_maximal+ 0] = -1;
	D[_H_][0*DP_maximal+ 0] = -1;
	D[_V_][0*DP_maximal+ 0] = -1;
	M[_S_][0*DP_maximal+ 0] = 0;
	M[_H_][0*DP_maximal+ 0] = WS_MIN;
	M[_V_][0*DP_maximal+ 0] = WS_MIN;
	for (i = 1; i < m; i++) 
	{
		D[_S_][i*DP_maximal+ 0] = _V_;
		D[_H_][i*DP_maximal+ 0] = _V_;
		D[_V_][i*DP_maximal+ 0] = _V_;
		M[_S_][i*DP_maximal+ 0] = WS_MIN;
		M[_H_][i*DP_maximal+ 0] = WS_MIN;
		M[_V_][i*DP_maximal+ 0] = i*GAP_HEAD1; //-(Params::GAP_OPEN + (i-1)*Params::GAP_EXT);
	}
	for (j = 1; j < n; j++) 
	{
		D[_S_][0*DP_maximal+ j] = _H_;
		D[_H_][0*DP_maximal+ j] = _H_;
		D[_V_][0*DP_maximal+ j] = _H_;
		M[_S_][0*DP_maximal+ j] = WS_MIN;
		M[_H_][0*DP_maximal+ j] = j*GAP_HEAD2; //-(Params::GAP_OPEN + (j-1)*Params::GAP_EXT);
		M[_V_][0*DP_maximal+ j] = WS_MIN;
	}
	//fill(firstSeq, secondSeq, distFunc);
	double gap_open;
	double gap_ext;
	double v1,v2,v3;
	double dist;
	for (i = 1; i < m; i++) 
	{
		for (j = 1; j < n; j++) 
		{
			//condition upper
			if(j==n-1)
			{
				gap_open=GAP_TAIL1;
				gap_ext=GAP_TAIL1;
			}
			else
			{
				gap_open=GAP_OPEN1;
				gap_ext=GAP_EXT1;
			}
			v1 = M[_V_][(i-1)*DP_maximal+ j] + gap_ext;
			v2 = M[_S_][(i-1)*DP_maximal+ j] + gap_open;
			v3 = M[_H_][(i-1)*DP_maximal+ j] + gap_open;
			M[_V_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_V_][i*DP_maximal+ j] == v1) D[_V_][i*DP_maximal+ j] = _V_;
			else if(M[_V_][i*DP_maximal+ j] == v2) D[_V_][i*DP_maximal+ j] = _S_;
			else D[_V_][i*DP_maximal+ j] = _H_;
			//condition left
			if(i==m-1)
			{
				gap_open=GAP_TAIL2;
				gap_ext=GAP_TAIL2;
			}
			else
			{
				gap_open=GAP_OPEN2;
				gap_ext=GAP_EXT2;
			}
			v1 = M[_H_][i*DP_maximal+ j-1] + gap_ext;
			v2 = M[_S_][i*DP_maximal+ j-1] + gap_open;
			v3 = M[_V_][i*DP_maximal+ j-1] + gap_open;
			M[_H_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_H_][i*DP_maximal+ j] == v1) D[_H_][i*DP_maximal+ j] = _H_;
			else if(M[_H_][i*DP_maximal+ j] == v2) D[_H_][i*DP_maximal+ j] = _S_;
			else D[_H_][i*DP_maximal+ j] = _V_;
			//condition diag
			dist = score.at((i-1)*IN_maximal+ j-1);  //Params::K - distFunc(firstSeq[i-1], secondSeq[j-1]);
			v1 = M[_V_][(i-1)*DP_maximal+ j-1] + dist;
			v2 = M[_H_][(i-1)*DP_maximal+ j-1] + dist;
			v3 = M[_S_][(i-1)*DP_maximal+ j-1] + dist;
			M[_S_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_S_][i*DP_maximal+ j] == v3) D[_S_][i*DP_maximal+ j] = _S_;
			else if (M[_S_][i*DP_maximal+ j] == v1) D[_S_][i*DP_maximal+ j] = _V_;
			else D[_S_][i*DP_maximal+ j] = _H_;
		}
	}
	//build(ali, firstSeq, secondSeq, distFunc);
	i = m-1;
	j = n-1;
	v1=M[_V_][i*DP_maximal+ j];
	v2=M[_H_][i*DP_maximal+ j];
	v3=M[_S_][i*DP_maximal+ j];
	double maximal = std::max(v1, std::max(v2, v3));
	int k = -1;
	if(v3==maximal)k = _S_;
	else if(v2==maximal)k = _H_;
	else k = _V_;
	//trace_back
	alignment.clear();
	int count = 0;
	int matches = 0;
	int cur_case=k;
	int pre_case;
	for(;;)
	{
		if(i==0||j==0)break;
		pre_case=D[cur_case][i*DP_maximal+ j];
		switch (cur_case)
		{
			case _S_:
				alignment.push_back(pair<int,int>(i,j)); 
				i--;
				j--;
				++matches;
				break;
			case _V_:
				alignment.push_back(pair<int,int>(i,-j)); 
				i--;
				break;
			case _H_:
				alignment.push_back(pair<int,int>(-i,j)); 
				j--;
				break;
			default:
				cout << "ERROR!! -> advance_global: invalid direction D[" << k << "](" << i << ", " << j << ") = " 
				<< D[k][i*DP_maximal+ j] << endl;
				exit(-1);
		}
		cur_case=pre_case;
		count++;
	}
	while (j> 0) alignment.push_back(pair<int,int>(-i,j)),j--;
	while (i> 0) alignment.push_back(pair<int,int>(i,0)), i--;
	reverse(alignment.begin(), alignment.end());
	ali_sco=maximal;
	return matches;
}

//---- get mapping -----//
//[note]: seqres should be longer than ami_, e.g., SEQRES vs ATOM
int Seqres_DynaProg(string &seqres,string &ami_,int *mapping)
{
	//--[0]check
	int len=(int)seqres.length();
	int totnum=(int)ami_.length();

	//--[1]dynamic_programming
	int i,j;
	int head=0;
	int n1=len;    //SEQRES
	int n2=totnum; //ATOM
	vector <double> score;
	score.resize(len*totnum);
	for(i=0;i<n1;i++)
	{
		for(j=0;j<n2;j++)
		{
			if(seqres[i]==ami_[j+head])score.at(i*n2+j)=10;
			else
			{
				if(seqres[i]=='X'||seqres[i]=='Z'||seqres[i]=='.')score.at(i*n2+j)=0;
				else if(ami_[j+head]=='X'||ami_[j+head]=='Z'||ami_[j+head]=='.')score.at(i*n2+j)=0;
				else score.at(i*n2+j)=-15;
			}
		}
	}
	double sco;
	int matchs;
	vector<pair<int,int> > WWW_alignment;
	matchs=WWW_Advance_Align_Dyna_Prog_Double(n1,n2,score,-11,-1,-110,-10,0,0,-110,-110,
		WWW_alignment,sco);
	int lcmp=(int)WWW_alignment.size();
	
	//extract
	for(i=0;i<len;i++)mapping[i]=-1; //default: NO
	int first,second;
	int retv=1;
	for(i=0;i<lcmp;i++)
	{
		first=WWW_alignment[i].first;
		second=WWW_alignment[i].second;
		if(first<=0)
		{
			if(second>0)
			{
				retv=-1;
				continue;
			}
		}
		if(first>0 && second>0)
		{
			mapping[first-1]=second-1;
		}
	}
	return retv;
}

//---- calculate match states ------//
//[note]: file_query should be equal or longer than a3m_query
int Calculate_Match_States(string &file_query,string &a3m_query,int *match)
{
	int l1=(int)file_query.length();
	int l2=(int)a3m_query.length();
	if(l1<l2)
	{
		fprintf(stderr,"ERROR!! file_query is shorter than a3m_query! [%d < %d]\n",l1,l2);
		return -1;
	}
	//transfer to upper case
	string file_query_=file_query;
	toUpperCase(file_query_);
	//do alignment
	int *mapping=new int[l1];
	int retv=Seqres_DynaProg(file_query_,a3m_query,mapping);
	if(retv!=1)
	{
		fprintf(stderr,"WARNING!! Seqres_DynaProg not return correctly ! \n");
	}
	//final assign
	int i;
	int count=0;
	for(i=0;i<l2;i++)match[i]=0;
	for(i=0;i<l1;i++)
	{
		int pos=mapping[i];
		if(pos!=-1)
		{
			char c=file_query[i];
			if(c>=65 && c<=90) //upper case
			{
				match[pos]=1;  //for match state
				count++;
			}
		}
	}
	//return match number
	delete [] mapping;
	return count;
}

//----- load PSSM file ------//
//-> load XXXX.mtx file
//    # Read in PSSM
//    # The PSSM file *.mtx contains one line for each column, beginning with line 15.
//    # The columns of these lines give the log-odds 3*100*log2(p(i,a)/f(a)) in the following order:
//    # ? A ? C D E F G H I K L M N P Q R S T V W X Y ? ? ? 
//-> [note]: return output is L*26
int Load_PSSM_File(string &infile,int canonical_length,double *output_fin)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",infile.c_str());
		return -1;
	}
	//skip
	for(int i=0;i<14;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file %s format error!\n",infile.c_str());
			return -1;
		}
	}
	//process
	int count=0;
	vector < vector <double> > output;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		//record
		vector <int> wsrec;
		int tmp_int;
		for(int i=0;i<26;i++)
		{
			if( ! (www>>tmp_int) )
			{
				fprintf(stderr,"file %s format error!\n",infile.c_str());
				return -1;
			}
			wsrec.push_back(tmp_int);
		}
		//calculate
		vector <double> wscalc;
		double tmp_double;
		for(int i=0;i<26;i++)
		{
			tmp_double=0.01*wsrec[ch2psi[i]]/3.0;
			wscalc.push_back(tmp_double);
		}
		//final assign
		output.push_back(wscalc);
		count++;
	}
	//final judge
	if(count!=canonical_length)
	{
		fprintf(stderr,"length not equal! [%d!=%d]\n",count,canonical_length);
		return -1;
	}
	//final calc
	for(int i=0;i<count;i++)for(int j=0;j<26;j++)output_fin[i*26+j]=output[i][j];
	return 1;
}


//-------- read in MSA in a3m format (i.e., normal FASTA with upper/lower) ------------//
//[note]: we set the first sequence as the query sequence,
//        that is to say, all the following sequences should be longer than the first
int Multi_FASTA_Input(string &multi_fasta,vector <string> &nam_list,vector <string> &fasta_list)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(multi_fasta.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",multi_fasta.c_str());
		return -1;
	}
	//load
	int relfirst=1;
	int firstlen;
	int first=1;
	int count=0;
	int number=0;
	string name;
	string seq;
	nam_list.clear();
	fasta_list.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf.length()>=1 && buf[0]=='>')
		{
			name=buf.substr(1,buf.length()-1);
			nam_list.push_back(name);
			count++;
			if(first!=1)
			{
				fasta_list.push_back(seq);
				number++;
				if(number==1)
				{
					firstlen=(int)seq.length();
				}
				else
				{
					int lowlen=getLowerCase(seq);
					int curlen=(int)seq.length()-lowlen;
					if(curlen!=firstlen)
					{
						fprintf(stderr,"length not equal at %s, [%d!=%d] \n",buf.c_str(),curlen,firstlen);
						return -1;
					}
				}
			}
			first=0;
			seq="";
		}
		else
		{
			if(first!=1)seq+=buf;
		}
	}
	//final
	if(first!=1)
	{
		fasta_list.push_back(seq);
		number++;
		if(number==1)
		{
			firstlen=(int)seq.length();
		}
		else
		{
			int lowlen=getLowerCase(seq);
			int curlen=(int)seq.length()-lowlen;
			if(curlen!=firstlen)
			{
				fprintf(stderr,"length not equal at %s, [%d!=%d] \n",buf.c_str(),curlen,firstlen);
				return -1;
			}
		}
	}
	//check
	if(number!=count)
	{
		fprintf(stderr,"num %d != count %d \n",number,count);
		return -1;
	}
	return count;
}


//================= Vice process =============//
int Get_Length(string &input)
{
	int i;
	int len=(int)input.length();
	int count=0;
	for(i=0;i<len;i++)
	{
		if(input[i]!='-')count++;
	}
	return count;
}

//============ remove gap-insert-gap ==========//
void Remove_GapInsGap(string &in,string &out)
{
	int i,k;
	int totlen=(int)in.length();
	int first=1;
	int start=-1;
	int len=-1;
	string tmp_in=in;
	for(i=0;i<totlen;i++)
	{
		if(tmp_in[i]>='a' && tmp_in[i]<='z') //insert
		{
			if(first==1)
			{
				first=0;
				start=i;
				len=1;
			}
			else
			{
				len++;
			}
		}
		else
		{
			if(first==0)
			{
				first=1;
				//judge
				if(start==0) //begin with insert
				{
					if(tmp_in[i]=='-')
					{
						for(k=0;k<len;k++)tmp_in[start+k]='@';
					}
				}
				else
				{
					if(tmp_in[start-1]=='-' && tmp_in[i]=='-')
					{
						for(k=0;k<len;k++)tmp_in[start+k]='@';
					}
				}
			}
		}
	}
	//terminal check
	if(first==0)
	{
		//judge
		if(start>0) //begin with insert
		{
			if(tmp_in[start-1]=='-')
			{
				for(k=0;k<len;k++)tmp_in[start+k]='@';
			}
		}
	}


	//---- head remove ----//
	int ii=-1;
	for(i=0;i<totlen;i++)
	{
		if(tmp_in[i]>='A' && tmp_in[i]<='Z')
		{
			ii=i;
			break;
		}
	}
	for(i=0;i<ii;i++)
	{
		if(tmp_in[i]>='a' && tmp_in[i]<='z')tmp_in[i]='@';
	}

	//---- tail remove ----//
	int jj=-1;
	for(i=totlen-1;i>=0;i--)
	{
		if(tmp_in[i]>='A' && tmp_in[i]<='Z')
		{
			jj=i;
			break;
		}
	}
	for(i=totlen-1;i>jj;i--)
	{
		if(tmp_in[i]>='a' && tmp_in[i]<='z')tmp_in[i]='@';
	}

	//final assign
	out="";
	for(i=0;i<totlen;i++)
	{
		if(tmp_in[i]!='@')out+=tmp_in[i];
	}
}


//-> Transfer pairwise alignment from a3m to normal fasta
//[note]: normal fasta means the alignment with equal length
int Transfer_Pairwise_Alignment(string &query_in,string &templa_in,string &query_out,string &templa_out)
{
	query_out="";
	templa_out="";
	int i;
	int l1=(int)query_in.length();
	int l2=(int)templa_in.length();
	int count=0;
	for(i=0;i<l2;i++)
	{
		//---- check ----//start
		if(count>l1)
		{
			fprintf(stderr,"impossible here !!! \n");
			return -1;
		}
		//---- check ----//over
		char c=templa_in[i];
		if(c>=97 && c<=122)  // lower case
		{
			query_out+='-';
			templa_out+=c;
		}
		else
		{
			query_out+=query_in[count];
			templa_out+=c;
			count++;
		}
	}
	//---- check ----//start
	if(count!=l1)
	{
		fprintf(stderr,"impossible here !!! \n");
		return -1;
	}
	//---- check ----//over
	query_out+='\0';
	templa_out+='\0';
	return l2;
}

//--------- get qstart and qend from normal fasta -------//
//[note]: 1-> the alignment of query and templa should be in the same length
//        2-> Xstart and Xend starts from 0
int Get_Qstart_Qend(string &query_in,string &templa_in,int L,int l1,int l2,
	int &qstart,int &qend,int &tstart,int &tend,
	string &query_out,string &templa_out,int &start,int &end)
{
	int i;
	start=-1;
	end=-1;
	//head init
	qstart=-1;
	tstart=-1;
	//head process
	for(i=0;i<L;i++)
	{
		char c1=query_in[i];
		char c2=templa_in[i];
		if( (c1>=65 && c1<=90) && (c2>=65 && c2<=90) ) //both upper case
		{
			qstart++;
			tstart++;
			break;
		}
		else
		{
			if( c1!='-' )qstart++;
			if( c2!='-' )tstart++;
		}
	}
	start=i;
	//tail init
	qend=l1;
	tend=l2;
	//tail process
	for(i=L-1;i>=0;i--)
	{
		char c1=query_in[i];
		char c2=templa_in[i];
		if( (c1>=65 && c1<=90) && (c2>=65 && c2<=90) ) //both upper case
		{
			qend--;
			tend--;
			break;
		}
		else
		{
			if( c1!='-' )qend--;
			if( c2!='-' )tend--;
		}
	}
	end=i;
	//final process
	query_out=query_in.substr(start,end-start+1);
	templa_out=templa_in.substr(start,end-start+1);
	return (end-start+1);
}

//---------- PruneHSP (With PSSM or BLOSUM) ----------//
//[note]: 1-> the alignment of query and templa should be in the same length of L',
//            i.e., after Get_Qstart_Qend() to extract the middle aligned region.
//        2-> qstart and qend starts from 0, maximun is qlen-1
int PruneHSP(string &query,string &templa,int L,
	int *match,int qstart,int qend,int qlen,int bg,double bs,double bl,
	double *input_matrix,double GAP,double EXTEND,int PSSM_or_BLOSUM)
{
	double smin;   // minimum score at current position $i
	double score;  // actual score at current position $i
	int i1;       // last pruned residue of HSP on  left side
	int i2;       // last pruned residue of HSP on right side
	int gap;      // gap 0: no gap currently open   1:gap opened
	int i;        // next column to read from pairwise alignment 
	int j;        // next position to read from query (sequence or profile)

	//[1] Count gaps in template that are aligned with match residues to the left/right of HSP
	int gapsleft=0;
	for(i=qstart-1;i>=0;i--)
	{
		if(match[i]==0)break;
		gapsleft++;
	}
	int gapsright=0;
	for(i=qend+1;i<qlen;i++)
	{
		if(match[i]==0)break;
		gapsright++;
	}

	//[2] determine the parameter
	double bleft,bright;
	if(gapsleft>=bg)bleft=bs;
	else bleft=bl;
	if(gapsright>=bg)bright=bs;
	else bright=bl;
	int proc_range=50;
//	if(PSSM_or_BLOSUM==1)proc_range=50;
//	else proc_range=20;

	//[3] Calculate scores at the end with input_matrix (i.e., either PSSM or BLOSUM)
	//[3-1] left side
	i1=-1;         // last pruned residue of HSP on left side
	if(bleft>-9)
	{
		gap=0;
		smin=0.0;
		score=0.0;
		i=0;      //-> i is the alignment position
		j=qstart; //-> j is the query position
		while(i<L && i<i1+proc_range)
		{
			smin+=bleft;
			if(query[i]=='-' || templa[i]=='-')  // gap in query or template sequence
			{
				if(!gap)
				{
					score-=GAP;
					gap=1;
				}
				else
				{
					score-=EXTEND;
				}
			}
			else      // match state
			{
				if(PSSM_or_BLOSUM==1)score+=input_matrix[ j*26 + (templa[i]-65) ];       //use PSSM
				else score+=input_matrix[ ch2i[query[i]-65]*21 + ch2i[templa[i]-65] ]; //use BLOSUM or GONNET
				gap=0;
			}
			//judge gap
			if (query[i]!='-')
			{
				j++;
			}
			//judge score
			if(score<smin)
			{
				i1=i;
				smin=0;
				score=0;
			}
			//increase alignment position
			i++;
		}
	}
	//[3-2] right side
	i2=L;          // last pruned residue of HSP on right side
	if(bright>-9)
	{
		gap=0;
		smin=0.0;
		score=0.0;
		i=L-1;      //-> i is the alignment position
		j=qend;     //-> j is the query position
		while(i>i1 && i>i2-proc_range)
		{
			smin+=bright;
			if(query[i]=='-' || templa[i]=='-')  // gap in query or template sequence
			{
				if(!gap)
				{
					score-=GAP;
					gap=1;
				}
				else
				{
					score-=EXTEND;
				}
			}
			else      // match state
			{
				if(PSSM_or_BLOSUM==1)score+=input_matrix[ j*26 + (templa[i]-65) ];       //use PSSM
				else score+=input_matrix[ ch2i[query[i]-65]*21 + ch2i[templa[i]-65] ]; //use BLOSUM or GONNET
				gap=0;
			}
			//judge gap
			if (query[i]!='-')
			{
				j--;
			}
			//judge score
			if(score<smin)
			{
				i2=i;
				smin=0;
				score=0;
			}
			//decrease alignment position
			i--;
		}
	}

	//[3-3] delete aligned positions on template
	for (i=0; i<=i1; i++)
	{
		if(templa[i]>='A' && templa[i]<='Z')templa[i]='-';
	}
	for (i=i2; i<L; i++)
	{
		if(templa[i]>='A' && templa[i]<='Z')templa[i]='-';
	}

	//return retained number of unpruned residues
	return i2-i1-1;
}


//---------- CheckScorePerColumn (With PSSM or BLOSUM) ----------//
//[note]: 1-> the alignment of query and templa should be in the same length of L
//        2-> qstart and qend starts from 0, maximun is qlen-1
int CheckScorePerColumn(string &query,string &templa,int L,
	int *match,int qstart,int qend,int qlen,double sc_thrshd,double score_min,
	double *input_matrix,double GAP,double EXTEND,int PSSM_or_BLOSUM)
{
	int i;        // next column to read from pairwise alignment 
	int j;        // next position to read from query (sequence or profile)
	double score=0.0;
	int gap=0;
	j=qstart;
	for(i=0;i<L;i++)
	{
		if(match[j])
		{
			if(query[i]=='-' || templa[i]=='-')  // gap in query or template sequence
			{
				if(!gap)
				{
					score-=GAP;
					gap=1;
				}
				else
				{
					score-=EXTEND;
				}
			}
			else      // match state
			{
				if(PSSM_or_BLOSUM==1)score+=input_matrix[ j*26 + (templa[i]-65) ];       //use PSSM
				else score+=input_matrix[ ch2i[query[i]-65]*21 + ch2i[templa[i]-65] ]; //use BLOSUM or GONNET
				gap=0;
			}
		}
		else
		{
			gap=0;
		}
		//judge gap
		if (query[i]!='-')
		{
			j++;
		}
	}

	//final judge
	int len=qend-qstart+1;
	double score_col=1.0*score/(len+0.1);
	if (score_col<sc_thrshd) return 0;
	if (score<score_min) return 0;
	return 1;
}


//=============== main process ==============//
//[note]: the query is from the first sequence in A3M input file
//        match is the MATCH state in query sequence
//-> need to prune identical sequence name, and bad length before
int Delete_Or_Not(string &query,string &templa,int query_len,int templa_len,int *match,int match_len,
	double seqid_thres,double coverage_thres,double sc_thrshd,double score_min,int bg,double bs,double bl,
	double *input_matrix,double GAP,double EXTEND,int PSSM_or_BLOSUM)
{
	//[1] transfer query and templa
	string query_cur,query_fin;
	string templa_cur,templa_fin;
	int L_cur=Transfer_Pairwise_Alignment(query,templa,query_cur,templa_cur);
	if(L_cur<=0)return -1;
	int qstart,qend,tstart,tend,start,end;
	int L_fin=Get_Qstart_Qend(query_cur,templa_cur,L_cur,query_len,templa_len,
		qstart,qend,tstart,tend,query_fin,templa_fin,start,end);
	if(L_fin<=0)return -1;
	//[2] prune HSP
	int retv;
	if( bs>-9 || bl>-9 )
	{
		retv=PruneHSP(query_fin,templa_fin,L_fin,match,
			qstart,qend,query_len,bg,bs,bl,
		input_matrix,GAP,EXTEND,PSSM_or_BLOSUM);
		if(retv==0)return -1;
	}
	//[3] calculate seqid and coverage
	int i;
	int j;
	int len=0;
	int qid=0;
	j=qstart;
	for(i=0;i<L_fin;i++)
	{
		if(query_fin[i]!='-')
		{
			if(templa_fin[i]!='-' && match[j]) // count only non-gap template residues in match columns!
			{
				len++;
				if(query_fin[i]==templa_fin[i])qid++;
			}
			j++;   // $j = next position in query
		}
	}
	if (len==0) return -1;
	if (100.0*qid/len<seqid_thres) return -1;
	if (100.0*len/match_len<coverage_thres) return -1;

	//[4] judge this template sequence
	if(sc_thrshd>-9 || score_min>0) 
	{
		retv=CheckScorePerColumn(query_fin,templa_fin,L_fin,match,
			qstart,qend,query_len,sc_thrshd,score_min,
			input_matrix,GAP,EXTEND,PSSM_or_BLOSUM);
		if(retv==0)return -1;
	}
	//[5] success
	string head_str=templa.substr(0,start);
	string tail_str=templa.substr(end+1,templa.length()-end-1);
	string templa_final=head_str+templa_fin+tail_str;
	templa=templa_final;
	return 1;
}

//-------- AlignHit -------//
void AlignHit(string &infile,string &outfile,string &seqfile,string &mtxfile,int best,
	double seqid_thres,double coverage_thres,double sc_thrshd,double score_min,int bg,double bs,double bl)
{
	int i;
	int retv;
	//load A3M infile
	vector <string> nam_list;
	vector <string> fasta_list;
	retv=Multi_FASTA_Input(infile,nam_list,fasta_list);
	if(retv<=0)exit(-1);
	if(retv<=2)
	{
//		fprintf(stderr,"WARNING !! infile %s contains %d sequences !! \n",infile.c_str(),retv);
		FILE *fp=fopen(outfile.c_str(),"wb");
		for(i=0;i<retv;i++)
		{
			fprintf(fp,">%s\n",nam_list[i].c_str());
			fprintf(fp,"%s\n",fasta_list[i].c_str());
		}
		fclose(fp);
		return;
	}
	//init
	map<string, int > ws_mapping;      //M, mapping the PDB's name
	map<string, int >::iterator iter;
	ws_mapping.clear();
	int totnum=retv;
	string query_seq=fasta_list[0];
	int query_len=(int)(fasta_list[0].length());
	Data_Structure_Init(query_len);
	vector <int> valid_list;
	valid_list.resize(totnum);
	vector <int> length_list;
	length_list.resize(totnum);
	for(i=0;i<totnum;i++)length_list[i]=Get_Length(fasta_list[i]);
	//load seqfile
	int match_len=query_len;
	int *match=QUERY_MATCH;
	if(seqfile!="")
	{
		string seqfile_query;
		retv=Read_FASTA_SEQRES(seqfile,seqfile_query,1); //->from .fasta file
		if(retv==-1)
		{
			fprintf(stderr,"WARNING !! seqfile %s bad, use canonical all match states !! \n",seqfile.c_str());
		}
		else
		{
			retv=Calculate_Match_States(seqfile_query,query_seq,QUERY_MATCH);
			if(retv==-1)fprintf(stderr,"WARNING !! seqfile %s bad, use canonical all match states !! \n",seqfile.c_str());
			else match_len=retv;			
		}
		//assign QUERY_SEQ
		QUERY_SEQ=query_seq;
		for(i=0;i<query_len;i++)
		{
			if(match[i]==0)QUERY_SEQ[i]=QUERY_SEQ[i]-'A'+'a';
		}
	}
	//load PSSM file
	int PSSM_or_BLOSUM=0;    //default: use BLOSUM to filter head and tail
	double *input_matrix=GONNET_CUR;
	double GAP=GAP_PENALTY;
	double EXTEND=EXTEND_PENALTY;
	if(mtxfile!="")
	{
		retv=Load_PSSM_File(mtxfile,query_len,PSSM_MTX);
		if(retv!=1)
		{
			fprintf(stderr,"WARNING !! mtxfile %s bad, use GONNET matrix !! \n",mtxfile.c_str());
		}
		else
		{
			PSSM_or_BLOSUM=1;    //now we use PSSM to filter head and tail
			input_matrix=PSSM_MTX;
		}
	}

	//---- check each template ------//
	int map_count=0;
	valid_list[0]=1;
	for(i=1;i<totnum;i++)
	{
		//check for redundant template sequence
		iter = ws_mapping.find(nam_list[i]);
		if(iter == ws_mapping.end()) //record
		{
			map_count++;
			ws_mapping.insert(map < string, int >::value_type(nam_list[i], map_count));
		}
		else   // find redundant template sequence
		{
			if(best==1) //only record one redundant template sequence
			{
				valid_list[i]=0;
				continue;
			}
		}
		//check for detailed property of the current template sequence
		retv=Delete_Or_Not(query_seq,fasta_list[i],query_len,length_list[i],match,match_len,
			seqid_thres,coverage_thres,sc_thrshd,score_min,bg,bs,bl,
			input_matrix,GAP,EXTEND,PSSM_or_BLOSUM);
		if(retv!=1)valid_list[i]=0;
		else valid_list[i]=1;
	}

	//----- output ------//
	FILE *fp=fopen(outfile.c_str(),"wb");
	if(seqfile!="")
	{
		fprintf(fp,">%s\n",nam_list[0].c_str());
		fprintf(fp,"%s\n",QUERY_SEQ.c_str());
	}
	for(i=1;i<totnum;i++)
	{
		if(valid_list[i]==1)
		{
			string out_str;
			Remove_GapInsGap(fasta_list[i],out_str);
			fprintf(fp,">%s\n",nam_list[i].c_str());
			fprintf(fp,"%s\n",out_str.c_str());
		}
	}
	fclose(fp);
}


/////////////////////////////////////////////////////////////////////////////////////
// Exit function
/////////////////////////////////////////////////////////////////////////////////////
void Usage()
{
	printf("Version: 1.03 \n");
	printf("AlignHits \n");
	printf("Input BLAST A3M MSA, output pruned A3M MSA. \n");
	printf("\n");
	printf("Usage: ./AlignHits -i input.a3m -o output.a3m -q input.seq                 \n");
	printf(" -i  <file>    : query MSA alignment (in a3m)                              \n");            //-> [ input file ]
	printf(" -j  <file>    : insert a2m-formatted query sequence into output alignment;\n");            //-> [ -q   file ]
	printf("                 upper/lower case determines match/insert columns\n");
	printf(" -o  <file>    : output MSA alignment (in a3m)                             \n");            //-> [ output file ]
	printf("[note]: the first sequence in input file should be the query sequence.     \n");
	printf("\n");
	printf("Options for main thresholds\n");
	printf(" -m percent    : minimum sequence identity to query in (default=0) percent \n");            //-> [ -qid percent ]
	printf("                 (seq-id = # identities in match columns / # hit residues in match\n");
	printf("                 columns)\n");
	printf(" -c percent    : minimum coverage in (default=0) percent      \n");                         //-> [ -cov coverage ]
	printf("Options for HSP prune thresholds (with query or -B alignment) \n");
	printf(" -g   int      : below this number of end gaps the lenient HSP pruning score is used,\n");  //-> [ -bg  int ]
	printf("               : above the strict score is employed (default=30)\n");
	printf(" -l   double   : lenient HSP pruning: min per-residue score in bits \n");                   //-> [ -bl  float ]
	printf("                 at ends of HSPs. Used when number of endgaps at the one end < bg    \n");
	printf("                 (see -bg) (default=-10)\n");
	printf(" -s   double   : strict HSP pruning: like -b, but used when number of endgaps >= bg  \n");  //-> [ -bs  float ]
	printf("                 (default=-10)\n");
	printf("Options for match column thresholds (with query or -P alignment)     \n");
	printf(" -p   double   : maximum p-value of HSP IN MATCH COLUMNS (default=1) \n");                  //-> [ -p   p-value (with -z)]
	printf(" -q   double   : minimum score per column in bits (default=-10)      \n");                  //-> [ -qsc value  (with -z)]
	printf("\n");
	printf("Other options:\n");
	printf(" -b   0/1      : extract only the best HSP per sequence (default=0)  \n");                  //-> [ -best ]
	printf(" -z   file     : read calculated PSSM matrix file, in .mtx format    \n");                  //-> [ -B/P   file ]
	printf("\n");
	printf("Examples: \n");
	printf("./AlignHits -i 1enh.a3m -o 1enh.a3m_out\n");
	printf("\n");
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
		string seq_file = "";
		string output_file = "";
		string pssm_file = "";
		int qid=0;
		int coverage=0;
		int bg=30;
		double bl=-10.0f;
		double bs=-10.0f;
		double pval=1.0f;
		double sc_thrshd=-10.0f;
		int best=0;

		//---- process argument ----//
		char c = 0;
		extern char* optarg;
		while ((c = getopt(argc, argv, "i:j:o:m:c:g:l:s:p:q:b:z:")) != EOF) 
		{
			switch (c) {
			case 'i':
				input_file = optarg;
				break;
			case 'j':
				seq_file = optarg;
				break;
			case 'o':
				output_file = optarg;
				break;
			case 'z':
				pssm_file = optarg;
				break;
			case 'm':
				qid = atoi(optarg);
				break;
			case 'c':
				coverage = atoi(optarg);
				break;
			case 'g':
				bg = atoi(optarg);
				break;
			case 'l':
				bl = atof(optarg);
				break;
			case 's':
				bs = atof(optarg);
				break;
			case 'p':
				pval = atof(optarg);
				break;
			case 'q':
				sc_thrshd = atof(optarg);
				break;
			case 'b':
				best = atoi(optarg);
				break;

			//--- default ----//
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
		if(qid<0 || qid>100)
		{
			fprintf(stderr,"  -m percent    : minimum sequence identity to query in (default=0) percent \n");
			exit(-1);
		}
		if(coverage<0 || coverage>100)
		{
			fprintf(stderr,"  -c percent    : minimum coverage in (default=0) percent \n");
			exit(-1);
		}
		if(best<0 || best>1)
		{
			fprintf(stderr,"  -b   0/1      : extract only the best HSP per sequence (default=0)\n");
			exit(-1);
		}

		//----- main process -----//
		double score_min=-3.0*log(fabs(pval))/log(2.0);
		AlignHit(input_file,output_file,seq_file,pssm_file,best,
			qid,coverage,sc_thrshd,score_min,bg,bs,bl);
		exit(0);
	}
}
