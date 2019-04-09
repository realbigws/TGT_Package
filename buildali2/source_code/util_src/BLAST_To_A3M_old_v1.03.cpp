#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <vector>
using namespace std;


//============= BLAST_To_A3M ===============//
//[purpose]: input pairwise BLAST result, output A3M

//each processing unit is as follows,
/*
> gi|159185141|ref|NP_355240.2| pseudoazurin [Agrobacterium tumefaciens 
str. C58]
  pseudoazurin [Agrobacterium sp. ATCC 31749]
  pseudoazurin [imported] - Agrobacterium tumefaciens (strain C58, 
Dupont)
  pseudoazurin [Agrobacterium tumefaciens str. C58]
  pseudoazurin [Agrobacterium sp. ATCC 31749]
Length=150

 Score =   146 bits (369),  Expect = 1e-42, Method: Compositional matrix adjust.
 Identities = 72/120 (60%), Positives = 92/120 (77%), Gaps = 2/120 (2%)

Query  3    IEVHMLNKGAEG-AMVFEPAYIKANPGDTVTFIPVDKGHNVESIKDMIPEGAEKFKSKIN  61
            IEV MLNKG++G AMVFEPA +KA  GD +TF+PVDKGH+  ++KDMIPEG  +FK K+N
Sbjct  27   IEVKMLNKGSDGQAMVFEPATVKAAVGDVITFVPVDKGHDAAAVKDMIPEGVAEFKGKMN  86

Query  62   ENYVLTVTQPGAYLVKCTPHYAMGMIALIAVGD-SPANLDQIVSAKKPKIVQERLEKVIA  120
            E   +TV + GAY+VKCTPH  MGMIAL+ VGD +PANLD + + K PK  ++RL + IA
Sbjct  87   EAVKITVEKEGAYVVKCTPHLGMGMIALVVVGDATPANLDVVKNGKLPKKARDRLNEEIA  146

*/

//--- check digit ---//
int Check_Digit(string &in_str)
{
	int i;
	int len=(int)in_str.length();
	int success=1;
	for(i=0;i<len;i++)
	{
		if( !  (in_str[i]>='0' && in_str[i]<='9') )
		{
			success=0;
			break;
		}
	}
	return success;
}

//--- pre_process string --//
int Pre_Process_String(
	int master_count,vector <string> &input,          //-> input data
	int &sco1,int &sco2,double &eval,                 //-> output score
	int &iden,int &posi,int &lali,int &gap,           //-> output lali
	string &nam,string &ali1_str,string &ali2_str,    //-> output alignment
	int &ali_start)                                   //-> output range

{
	//init
	sco1=-1;
	sco2=-1;
	eval=-1;
	iden=-1;
	posi=-1;
	lali=-1;
	gap=-1;
	nam="";
	ali1_str="";
	ali2_str="";
	ali_start=-1;
	//process
	int i,k,l;
	int cur;
	int size=(int)input.size();
	string content,buf,temp;
	//process name
	for(i=0;i<size;i++)
	{
		if(input[i]=="")continue;
		if(input[i][0]=='>')
		{
			content=input[i].substr(1,input[i].length()-1);
			istringstream www(content);
			www>>nam;
			//judge name
			if(nam=="")
			{
				fprintf(stderr,"name process error here !! %d \n",master_count);
				return -1;
			}
			//break
			cur=i;
			break;
		}
	}
	//process score
	for(i=cur+1;i<size;i++)
	{
		if(input[i]=="")continue;
		content=input[i];
		istringstream www(content);
		www>>temp;
		if(temp=="Score")
		{
			//get sco1,sco2,eval
			www>>temp>>sco1>>temp>>buf;
			temp=buf.substr(1,buf.length()-3);
			sco2=atoi(temp.c_str());
			www>>temp>>temp>>buf;
			temp=buf.substr(0,buf.length()-1);
			eval=atof(temp.c_str());
			//get iden,posi,lali,gap
			int cur_=i;
			i++;
			int cur_success=0;
			string content_;
			for(k=cur_+1;k<size;k++,i++)
			{
				if(input[k]=="")continue;
				content_=input[k];
				istringstream sss(content_);
				sss>>temp;
				if(temp=="Identities")
				{
					// iden/lali
					sss>>temp>>buf;
					for(l=0;l<(int)buf.length();l++)
					{
						if(buf[l]=='/')break;
					}
					temp=buf.substr(0,l);
					iden=atoi(temp.c_str());
					temp=buf.substr(l+1,(int)buf.length()-l-1);
					lali=atoi(temp.c_str());
					// posi/lali
					sss>>temp>>temp>>temp>>buf;
					for(l=0;l<(int)buf.length();l++)
					{
						if(buf[l]=='/')break;
					}
					temp=buf.substr(0,l);
					posi=atoi(temp.c_str());
					// gap/lali
					sss>>temp>>temp>>temp>>buf;
					for(l=0;l<(int)buf.length();l++)
					{
						if(buf[l]=='/')break;
					}
					temp=buf.substr(0,l);
					gap=atoi(temp.c_str());
					//final
					cur=k;
					cur_success=1;
					break;
				}
			}
			//break judge
			if(cur_success==1)break;
			else
			{
				fprintf(stderr,"score process error here !! %d,%s \n",master_count,nam.c_str());
				return -1;
			}
		}
	}
	//process alignment
	int first=1;
	for(i=cur+1;i<size;i++)
	{
		if(input[i]=="")continue;
		content=input[i];
		istringstream www(content);
		www>>temp;
		//judge name
		if(temp=="Query:")
		{
			www>>buf;
			if(first==1)
			{
				ali_start=atoi(buf.c_str());
				if(ali_start>0) first=0;
			}
//			www>>temp;
			buf=content.substr(8,content.length()-8);
			istringstream zzz(buf);
			string tttmp;
			zzz>>tttmp;
			//check digit
			int check_digit=Check_Digit(tttmp);
			if(check_digit==0)ali1_str+=tttmp;
			else
			{
				zzz>>tttmp;
				ali1_str+=tttmp;
			}
		}
		if(temp=="Sbjct:")
		{
//			www>>temp>>temp;
			buf=content.substr(8,content.length()-8);
			istringstream zzz(buf);
			string tttmp;
			zzz>>tttmp;
			//check digit
			int check_digit=Check_Digit(tttmp);
			if(check_digit==0)ali2_str+=tttmp;
			else
			{
				zzz>>tttmp;
				ali2_str+=tttmp;
			}
		}
	}
	//judge alignment
	if(ali1_str.length()!=ali2_str.length())
	{
		fprintf(stderr,"alignment length [%d!=%d] not equal here !! %d,%s \n",
			 (int)ali1_str.length(),(int)ali2_str.length(),master_count,nam.c_str());
		return -1;
	}
	//success return
	return 0;  //success
}

//========= main process ==========//
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

//-------- kill gap ----------//
void Kill_Gap(string &input,string &output)
{
	int i;
	int len=(int)input.length();
	output.clear();
	for(i=0;i<len;i++)
	{
		if(input[i]!='-')output.push_back(input[i]);
	}
}

//-------- post process ---------//
void Post_Process(string &ali1_str,string &ali2_str,
	int start,int length,string &out_str)
{
	int i;
	out_str="";
	//calc len
	string pure_str;
	Kill_Gap(ali1_str,pure_str);
	int len=(int)pure_str.length();
	//add head and tail gap
	string head_gap="";
	for(i=0;i<start;i++)head_gap+='-';
	string tail_gap="";
	for(i=0;i<length-start-len;i++)tail_gap+='-';
	//turn inserts into smaller case
	int size=(int)ali1_str.length();
	for(i=0;i<size;i++)
	{
		if(ali1_str[i]=='-')ali2_str[i]=ali2_str[i]+32;  // to lower case
	}
	//final
	out_str=head_gap+ali2_str+tail_gap;
}


//[note]: we should consider the following "no hit found" cases
/*
***** No hits found *****
*/
int Main_Process(string &infile,string &query,
	vector <string> &nam_rec,vector <string> &seq_rec)
{
	//init
	nam_rec.clear();
	seq_rec.clear();
	//read
	ifstream fin;
	string buf,temp,bof,tomp;
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",infile.c_str());
		return -1;
	}
	//data structure
	int sco1,sco2;
	double eval;                    //-> output score
	int iden,posi,lali,gap;         //-> output lali
	string nam,ali1_str,ali2_str;   //-> output alignment
	int ali_start;                  //-> output range
	//skip to find the first QUERY
	vector <string> temp_rec;
	temp_rec.clear();
	int first=1;
	int sco_duplicate=0;
	int rel_next=0;
	int master_count=0;
	int query_len=(int)query.length();
	int retv;
	string name_temp_rec;
	for(;;)
	{
		if(!getline(fin,bof,'\n'))
		{
			fprintf(stderr,"file %s format bad!\n",infile.c_str());
			return -1;
		}
		if(bof=="")continue;
		if(bof==" ***** No hits found ******")
		{
//			fprintf(stderr,"file %s contains no hits!\n",infile.c_str());
			nam_rec.push_back("QUERY");
			seq_rec.push_back(query);
			return 1;
		}

		//real process
		if(bof[0]=='>') //start to record
		{
			buf=bof;
			name_temp_rec=buf;
			//push_back
			nam_rec.push_back("QUERY");
			seq_rec.push_back(query);
gobak:
			//process previous record
			if(first==1)first=0;
			else
			{
				//proc
				retv=Pre_Process_String(master_count,temp_rec,sco1,sco2,eval,
					iden,posi,lali,gap,nam,ali1_str,ali2_str,ali_start);
				if(retv==0)
				{
					string out_str;
					Post_Process(ali1_str,ali2_str,ali_start-1,query_len,out_str);
					//push_back
					nam_rec.push_back(nam);
					seq_rec.push_back(out_str);
				}
				//clear
				temp_rec.clear();
				if(rel_next==1)sco_duplicate=0;
				if(sco_duplicate>1)temp_rec.push_back(name_temp_rec);
				//count
				master_count++;
			}

			//record current
			temp_rec.push_back(buf);
			for(;;)
			{
				if(!getline(fin,buf,'\n'))
				{
					fprintf(stderr,"file %s format bad!\n",infile.c_str());
					return -1;
				}
				if(buf=="")continue;
				int lambda_len=(int)buf.length();
				if(lambda_len<6)continue;
				string lambda_str=buf.substr(0,6);
				istringstream lll(buf);
				string lll1,lll2;
				lll>>lll1>>lll2;
				if(lambda_str=="Lambda" && lll1=="Lambda" && lll2=="K") //final process
				{
					//proc
					retv=Pre_Process_String(master_count,temp_rec,sco1,sco2,eval,
						iden,posi,lali,gap,nam,ali1_str,ali2_str,ali_start);
					if(retv==0)
					{
						string out_str;
						Post_Process(ali1_str,ali2_str,ali_start-1,query_len,out_str);
						//push_back
						nam_rec.push_back(nam);
						seq_rec.push_back(out_str);
					}
					//clear
					temp_rec.clear();
					//count
					master_count++;
					//goto end
					goto end;
				}
				istringstream www(buf);
				www>>temp;
				if(temp=="Score")
				{
					sco_duplicate++;
					if(sco_duplicate>1)
					{
						rel_next=0;
						goto gobak;
					}
				}
				if(buf[0]=='>')
				{
					name_temp_rec=buf;
					rel_next=1;
					goto gobak;
				}
				temp_rec.push_back(buf);
			}
		}
	}
end:
	return (int)nam_rec.size();
}


//========= validate sequence ==========//
int Ori_AA_Map_WS[26]=
{ 0,20,2,3,4,5,6,7,8,20,10,11,12,13,20,15,16,17,18,19,20, 1, 9,20,14,20};
// A B C D E F G H I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
// 0 1 2 3 4 5 6 7 8  9 10 11 12 14 14 15 16 17 18 19 20 21 22 23 24 25

void Validate_Sequence(string &instr,string &outstr)
{
	int i;
	int len=(int)instr.length();
	outstr=instr;
	for(i=0;i<len;i++)
	{
		if(instr[i]=='-')continue;
		char a=instr[i];
		if(a>='a' && a<='z')continue;
		if(a<'A' || a>='Z')
		{
			outstr[i]='X';
			continue;
		}
		int retv=Ori_AA_Map_WS[a-'A'];
		if(retv==20)
		{
			outstr[i]='X';
			continue;
		}
	}
}


//-------- main -------//
int main(int argc,char **argv)
{
	//------ BLAST_To_A3M -------//
	{
		if(argc<4)
		{
			fprintf(stderr,"Version: 1.03 \n");
			fprintf(stderr,"BLAST_To_A3M <query_fasta> <blast_input> <a3m_output> \n");
			fprintf(stderr,"[note]: BLAST out is pariwise mode. \n");
			fprintf(stderr,"        Must use new version of blastpgp (say, 2.2.26 ) \n");
			exit(-1);
		}
		string query_fasta=argv[1];
		string blast_input=argv[2];
		string a3m_output=argv[3];
		//process
		int retv;
		//-> read query fasta
		string query;
		retv=Read_FASTA_SEQRES(query_fasta,query);
		if(retv<=0)exit(-1);
		//-> process blast input
		vector <string> nam_rec;
		vector <string> seq_rec;
		retv=Main_Process(blast_input,query,nam_rec,seq_rec);
		if(retv<=0)exit(-1);
		//-> output
		FILE *fp=fopen(a3m_output.c_str(),"wb");
		for(int i=0;i<(int)nam_rec.size();i++)
		{
			fprintf(fp,">%s\n",nam_rec[i].c_str());
			string outstr;
			Validate_Sequence(seq_rec[i],outstr);
			fprintf(fp,"%s\n",outstr.c_str());
		}
		fclose(fp);
		//exit
		exit(0);
	}
}


