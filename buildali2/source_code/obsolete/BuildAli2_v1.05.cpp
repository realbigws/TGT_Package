#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "getopt.h"
using namespace std;



//--- macro ---//
int NEW_or_OLD;  // new -> psiblast [1]; old-> blastpgp [0]

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


//======= get PSI/A3M number ========//
int Get_PSI_Number(string &multi_fasta)
{
	ifstream fin;
	string buf,temp;
	fin.open(multi_fasta.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",multi_fasta.c_str());
		return -1;
	}
	int firstlen=0;
	int first=1;
	string seq;
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		seq=buf.substr(33,buf.length()-33);
		count++;
		if(first==1)
		{
			firstlen=(int)seq.length();
			first=0;
		}
		else
		{
			int curlen=(int)seq.length();
			if(curlen!=firstlen)
			{
				fprintf(stderr,"length not equal at %s, [%d!=%d] \n",buf.c_str(),curlen,firstlen);
				return -1;
			}
		}
	}
	return count;
}
int Get_A3M_Number(string &multi_fasta)
{
	ifstream fin;
	string buf,temp;
	fin.open(multi_fasta.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",multi_fasta.c_str());
		return -1;
	}
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='>')count++;
	}
	return count;
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


//================ BuildAli2 ===============//
//[purpose]:
//1. Given input fasta sequence, output K iteration psiblast to construct A3M.
//2. We may keep the A3M file from each iteration.
//3. Given A3M file, we could generate .ss2, .psp, .pssm, and .chk files.

//[requirement]:
//1. MSA_To_PSSM
//2. AlignHits
//3. psiblast (from BLAST+ 2.2.29)
//4. formatted NR90, NR70 database

//[note]:
//1. We omit PSIPRED calculation here.
//2. For canonical A3M, (i.e., contains SSE and SSE_conf), try other package.
//3. To correctly invoke AlignHits, we should provide the .mtx file from the core A3M.


// main part for PSI-BLAST search  -- added by Sheng Wang at 2014.02.01
int BuildAli2_Main(int max_iter,int cpu_num,
	string &nr90_root,string &nr70_root,
	string &query_file,string &out_dir,string &util_dir,int break_or_not)
{
	int iter;
	int retv;
	char command[30000];
	string query_name;
	getBaseName(query_file,query_name,'/','.');

	//----- assign parameter -----//
	//-> blast related
	int nmax=10000;            //-> alignments to break BuildAli2
	int num_alignments=20000;  //-> alignments to output during BLAST
	double evalue=0.001;       //-> BLAST e-value threshold
	//-> MSA_To_PSSM related
	int cdhit_cutnum=20;       //-> if MSA number is less than cdhit_cutnum, then omit cd-hit
	//-> alignhits related
	//core MSA
	int cov_core=80;           //-> AlignHits coverage for core MSA
	double bl_core=0.33;       //-> AlignHits lenient HSP pruning for core MSA
	double bs_core=0.67;       //-> AlignHits strict HSP pruning for core MSA
	int bg_core=30;            //-> AlignHits gap pruning for core MSA
	//normal MSA
	int cov_normal=20;         //-> AlignHits coverage for normal MSA
	double bl_normal=0.0;      //-> AlignHits lenient HSP pruning for core MSA
	double bs_normal=0.167;    //-> AlignHits strict HSP pruning for core MSA
	int bg_normal=20;          //-> AlignHits gap pruning for core MSA


	//---- load sequence file -----//
	string query_sequence;
	retv=Read_FASTA_SEQRES(query_file,query_sequence);
	if(retv<=0)
	{
		fprintf(stderr,"seqfile %s contains no query_sequence !! \n",query_file.c_str());
		return -1;
	}


	//------- if max_iter==1 -----------//
	int nhits_prev=-1;
	int nhits;
	for(iter=1;iter<=max_iter;iter++)
	{

//wstest
printf("%s -> iter %d/%d \n",query_name.c_str(),iter,max_iter);

		//----------- run psiblast ----------//
		//-> 1. get nr and run blast
		string nr_root;
		if(iter==1)  //using query_file
		{
			nr_root=nr90_root;
			if(NEW_or_OLD==1)
			{
				sprintf(command,"%s/psiblast -show_gis -use_sw_tback -num_threads %d -num_alignments %d -num_descriptions 0 -evalue %lf -db %s -query %s -out %s/%s.%d.blast -comp_based_stats 1 1> ws1 2> ws2",
					util_dir.c_str(),cpu_num,num_alignments,evalue,nr_root.c_str(),query_file.c_str(),out_dir.c_str(),query_name.c_str(),iter);
			}
			else
			{
				sprintf(command,"%s/BLAST/bin/blastpgp -I T -s T -a %d -b %d -v 0 -e %lf -d %s -i %s -o %s/%s.%d.blast 1> ws1 2> ws2",
					util_dir.c_str(),cpu_num,num_alignments,evalue,nr_root.c_str(),query_file.c_str(),out_dir.c_str(),query_name.c_str(),iter);
			}
		}
		else         //using pssm_file from previous round 
		{
			nr_root=nr70_root;
			if(NEW_or_OLD==1)
			{
				sprintf(command,"%s/psiblast -show_gis -use_sw_tback -num_threads %d -num_alignments %d -num_descriptions 0 -evalue %lf -db %s -in_pssm %s/%s.chk -out %s/%s.%d.blast -comp_based_stats 1 1> ws1 2> ws2",
					util_dir.c_str(),cpu_num,num_alignments,evalue,nr_root.c_str(),out_dir.c_str(),query_name.c_str(),out_dir.c_str(),query_name.c_str(),iter);
			}
			else
			{
				sprintf(command,"%s/BLAST/bin/blastpgp -I T -s T -a %d -b %d -v 0 -e %lf -d %s -i %s -R %s/%s.chk -q 1 -o %s/%s.%d.blast 1> ws1 2> ws2",
					util_dir.c_str(),cpu_num,num_alignments,evalue,nr_root.c_str(),query_file.c_str(),out_dir.c_str(),query_name.c_str(),out_dir.c_str(),query_name.c_str(),iter);
			}
		}
		retv=system(command);
		if(retv!=0)
		{
			fprintf(stderr,"seqfile %s failed at psiblast on %d/%d iteration !! \n",query_file.c_str(),iter,max_iter);
			return -1;
		}
		sprintf(command,"rm -f ws1 ws2 error.log");
		retv=system(command);


		//-> 2. run BLAST_To_A3M
		if(NEW_or_OLD==1)
		{
			sprintf(command,"%s/BLAST_To_A3M_new %s %s/%s.%d.blast %s/%s.%d.a3m_ini",
				util_dir.c_str(),query_file.c_str(),out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str(),iter);
		}
		else
		{
			sprintf(command,"%s/BLAST_To_A3M_old %s %s/%s.%d.blast %s/%s.%d.a3m_ini",
				util_dir.c_str(),query_file.c_str(),out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str(),iter);
		}
		retv=system(command);
		if(retv!=0)
		{
			fprintf(stderr,"seqfile %s failed at BLAST_To_A3M on %d/%d iteration !! \n",query_file.c_str(),iter,max_iter);
			return -1;
		}

		//-> check for zero or single hit
		{
			sprintf(command,"%s/%s.%d.a3m_ini",out_dir.c_str(),query_name.c_str(),iter);
			string a3m_file=command;
			int number=Get_A3M_Number(a3m_file);
			if(number==1)
			{
				if(iter==1)
				{
					sprintf(command,"cat %s/%s.%d.a3m_ini > %s/%s.a3m",
						out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str());
					retv=system(command);
				}
				else
				{
					sprintf(command,"cat %s/%s.%d.a3m_ini >> %s/%s.a3m",
						out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str());
					retv=system(command);
				}
				sprintf(command,"cat %s/%s.a3m > %s/%s.a3m.%d",
					out_dir.c_str(),query_name.c_str(),out_dir.c_str(),query_name.c_str(),iter);
				retv=system(command);
				break;
			}
		}

//printf("BLAST_To_A3M done ! \n");

		//-> 3. if max_iter==1, then break
		if(max_iter==1)
		{
			//-> alignhits to filter
			if(NEW_or_OLD==1)
			{
				sprintf(command,"%s/AlignHits_new -b 0 -c %d -l %lf -s %lf -g %d -j %s -i %s/%s.%d.a3m_ini -o %s/%s.%d.a3m_fil",
					util_dir.c_str(),cov_normal,bl_core,bs_core,bg_normal,query_file.c_str(),out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str(),iter);
			}
			else
			{
				sprintf(command,"%s/AlignHits_old -b 0 -c %d -l %lf -s %lf -g %d -j %s -i %s/%s.%d.a3m_ini -o %s/%s.%d.a3m_fil",
					util_dir.c_str(),cov_normal,bl_core,bs_core,bg_normal,query_file.c_str(),out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str(),iter);
			}
			retv=system(command);
			if(retv!=0)
			{
				fprintf(stderr,"seqfile %s failed at AlignHits Max_Iter=1 on %d/%d iteration !! \n",query_file.c_str(),iter,max_iter);
				return -1;
			}
			//-> generate current iterate a3m file
			sprintf(command,"cat %s/%s.%d.a3m_fil > %s/%s.a3m",
				out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str());
			retv=system(command);
			sprintf(command,"cat %s/%s.a3m > %s/%s.a3m.%d",
				out_dir.c_str(),query_name.c_str(),out_dir.c_str(),query_name.c_str(),iter);
			retv=system(command);
			break;
		}

		//------- generate core region ------//
		if(iter==1) 
		{
			//-> 1. generate core_a3m
			if(NEW_or_OLD==1)
			{
				sprintf(command,"%s/AlignHits_new -b 1 -c %d -l %lf -s %lf -g %d -j %s -i %s/%s.%d.a3m_ini -o %s/%s.core_a3m",
					util_dir.c_str(),cov_core,bl_core,bs_core,bg_core,query_file.c_str(),out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str());
			}
			else
			{
				sprintf(command,"%s/AlignHits_old -b 1 -c %d -l %lf -s %lf -g %d -j %s -i %s/%s.%d.a3m_ini -o %s/%s.core_a3m",
					util_dir.c_str(),cov_core,bl_core,bs_core,bg_core,query_file.c_str(),out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str());
			}
			retv=system(command);
			if(retv!=0)
			{
				fprintf(stderr,"seqfile %s failed at AlignHits Core on %d/%d iteration !! \n",query_file.c_str(),iter,max_iter);
				return -1;
			}
			//-> 2. generate core_psi
			sprintf(command,"%s/A3M_To_PSI %s/%s.core_a3m %s/%s.core_psi",
				util_dir.c_str(),out_dir.c_str(),query_name.c_str(),out_dir.c_str(),query_name.c_str());
			retv=system(command);
			if(retv!=0)
			{
				fprintf(stderr,"seqfile %s failed at A3M_To_PSI Core on %d/%d iteration !! \n",query_file.c_str(),iter,max_iter);
				return -1;
			}
			//-> 3. generate core_mtx
			if(NEW_or_OLD==1)
			{
				sprintf(command,"%s/MSA_To_PSSM -i %s/%s.core_psi -m %s/%s.core_mtx -c %d",
					util_dir.c_str(),out_dir.c_str(),query_name.c_str(),out_dir.c_str(),query_name.c_str(),cdhit_cutnum);
			}
			else
			{
				sprintf(command,"%s/BLAST/bin/blastpgp -d %s/dummydb -i %s -B %s/%s.core_psi -C %s.chk 1> ws1 2> ws2",
					util_dir.c_str(),util_dir.c_str(),query_file.c_str(),out_dir.c_str(),query_name.c_str(),query_name.c_str());
			}
			retv=system(command);
			if(retv!=0)
			{
				fprintf(stderr,"seqfile %s failed at MSA_To_PSSM Core on %d/%d iteration !! \n",query_file.c_str(),iter,max_iter);
				return -1;
			}

			//->3.1 blasgpgp process only
			if(NEW_or_OLD!=1)
			{
				sprintf(command,"cp %s %s.seq_makemat; echo %s.chk > %s.pn; echo %s.seq_makemat > %s.sn; %s/BLAST/bin/makemat -P %s; mv %s.mtx %s/%s.core_mtx",
					query_file.c_str(),query_name.c_str(),query_name.c_str(),query_name.c_str(),query_name.c_str(),query_name.c_str(),
					util_dir.c_str(),query_name.c_str(),query_name.c_str(),out_dir.c_str(),query_name.c_str());
				retv=system(command);
				if(retv!=0)
				{
					fprintf(stderr,"seqfile %s failed at MSA_To_PSSM Makemat on %d/%d iteration !! \n",query_file.c_str(),iter,max_iter);
					return -1;
				}
				sprintf(command,"rm -f %s.seq_makemat %s.pn %s.sn %s.chk %s.mn %s.aux ws1 ws2 error.log",
					query_name.c_str(),query_name.c_str(),query_name.c_str(),query_name.c_str(),query_name.c_str(),query_name.c_str());
				retv=system(command);
			}
		}


		//------ filter a3m_ini by core_mtx ------//
		int best_cur;
		if(iter==1)best_cur=1;
		else best_cur=0;
		int cov_cur;
		if(iter==1)cov_cur=cov_core;
		else cov_cur=cov_normal;
		//-> 1. filter a3m_ini to a3m_fil
		if(NEW_or_OLD==1)
		{
			sprintf(command,"%s/AlignHits_new -b %d -c %d -l %lf -s %lf -g %d -j %s -z %s/%s.core_mtx -i %s/%s.%d.a3m_ini -o %s/%s.%d.a3m_fil",
				util_dir.c_str(),best_cur,cov_cur,bl_normal,bs_normal,bg_normal,query_file.c_str(),out_dir.c_str(),query_name.c_str(),
				out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str(),iter);
		}
		else
		{
			sprintf(command,"%s/AlignHits_old -b %d -c %d -l %lf -s %lf -g %d -j %s -z %s/%s.core_mtx -i %s/%s.%d.a3m_ini -o %s/%s.%d.a3m_fil",
				util_dir.c_str(),best_cur,cov_cur,bl_normal,bs_normal,bg_normal,query_file.c_str(),out_dir.c_str(),query_name.c_str(),
				out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str(),iter);
		}
		retv=system(command);
		if(retv!=0)
		{
			fprintf(stderr,"seqfile %s failed at AlignHits on %d/%d iteration !! \n",query_file.c_str(),iter,max_iter);
			return -1;
		}
		//-> 2. generate psi_fil
		sprintf(command,"%s/A3M_To_PSI %s/%s.%d.a3m_fil %s/%s.%d.psi_fil",
			util_dir.c_str(),out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str(),iter);
		retv=system(command);
		if(retv!=0)
		{
			fprintf(stderr,"seqfile %s failed at A3M_To_PSI on %d/%d iteration !! \n",query_file.c_str(),iter,max_iter);
			return -1;
		}

		//------ generate checkpoint file for the next iteration ------//
		//------- check for termination -------//
		sprintf(command,"%s/%s.%d.psi_fil",out_dir.c_str(),query_name.c_str(),iter);
		string multi_fasta=command;
		nhits=Get_PSI_Number(multi_fasta);

printf("%s -> nhits = %d, nhits_prev = %d \n",query_name.c_str(),nhits,nhits_prev);



//-> early break
if(NEW_or_OLD==0)  //old (blasgpgp), then early break
{
		if(nhits==1)
		{
			if(iter==1)
			{
				sprintf(command,"cat %s/%s.%d.a3m_fil > %s/%s.a3m",
					out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str());
				retv=system(command);
			}
			break;
		}
		if(iter>1)
		{
			if(nhits<=nhits_prev && break_or_not==0)break;
			if(nhits>=nmax)break;
		}
		nhits_prev=nhits;
}

		//-> 1. cat to %.psi file
		if(iter==1)
		{
			sprintf(command,"cat %s/%s.%d.psi_fil > %s/%s.psi",
				out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str());
		}
		else
		{
			sprintf(command,"cat %s/%s.%d.psi_fil >> %s/%s.psi",
				out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str());
		}
		retv=system(command);
		//-> 2. cat to %.a3m file
		if(iter==1)
		{
			sprintf(command,"cat %s/%s.%d.a3m_fil > %s/%s.a3m",
				out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str());
		}
		else
		{
			sprintf(command,"cat %s/%s.%d.a3m_fil >> %s/%s.a3m",
				out_dir.c_str(),query_name.c_str(),iter,out_dir.c_str(),query_name.c_str());
		}
		retv=system(command);
		//-> 3. generate %s.chk file from %s.psi file
		if(NEW_or_OLD==1)
		{
			sprintf(command,"%s/MSA_To_PSSM -i %s/%s.psi -t %s/%s.chk -c %d",
				util_dir.c_str(),out_dir.c_str(),query_name.c_str(),out_dir.c_str(),query_name.c_str(),cdhit_cutnum);
		}
		else
		{
			sprintf(command,"%s/BLAST/bin/blastpgp -d %s/dummydb -i %s -B %s/%s.psi -C %s/%s.chk -u 1 -J 1> ws1 2> ws2",
				util_dir.c_str(),util_dir.c_str(),query_file.c_str(),out_dir.c_str(),query_name.c_str(),out_dir.c_str(),query_name.c_str());
		}
		retv=system(command);
		if(retv!=0)
		{
			fprintf(stderr,"seqfile %s failed at MSA_To_PSSM on %d/%d iteration !! \n",query_file.c_str(),iter,max_iter);
			return -1;
		}
		sprintf(command,"rm -f ws1 ws2 error.log");
		retv=system(command);

		//-> 4. generate current iterate a3m file
		sprintf(command,"cat %s/%s.a3m > %s/%s.a3m.%d",
			out_dir.c_str(),query_name.c_str(),out_dir.c_str(),query_name.c_str(),iter);
		retv=system(command);



//-> late break
if(NEW_or_OLD==1)  //new (psiblast), then late break
{
		if(nhits==1)break;
		if(iter>1)
		{
			if(nhits<=nhits_prev && break_or_not==0)break;
			if(nhits>=nmax)break;
		}
		nhits_prev=nhits;
}

	}
	//return
	if(iter>max_iter)iter=max_iter;
	return iter;
}


//---------- usage ---------//
void Usage() 
{
	fprintf(stderr,"Version: 1.05 \n");
	fprintf(stderr,"BuildAli2 -i query_file [-o a3m_out] [-d nr90] [-D nr70] [-u util_root] \n");
	fprintf(stderr,"          [-m max_iter] [-c cpu_num] [-r tmp_root] [-k remove] [-b] [-n] \n\n");
	fprintf(stderr,"Usage : \n\n");
	fprintf(stderr,"-i query_file :        Input query file in FASTA format. \n\n");
	fprintf(stderr,"-o a3m_out :           Output A3M file. ( default = <query_name>.a3m ) \n\n");
	fprintf(stderr,"-d nr90 :              Location of formatted nr90 database by formatdb. \n");
	fprintf(stderr,"                       ( default = databases/NR_new/nr90 ) \n\n");
	fprintf(stderr,"-D nr70 :              Location of formatted nr70 database by formatdb. \n");
	fprintf(stderr,"                       ( default = databases/NR_new/nr70 ) \n\n");
	fprintf(stderr,"-u util_root :         Directory for Utility Functions. (default=util/ ) \n\n");
	fprintf(stderr,"-m max_iter :          Maximal iteration. (default=5) \n\n");
	fprintf(stderr,"-c cpu_num :           CPU number. (default=1) \n\n");
	fprintf(stderr,"-r tmp_root :          Directory to put temporary files. ( default=tmp/ ) \n\n");
	fprintf(stderr,"-k remove :            If specified, then remove all temporary files. \n");
	fprintf(stderr,"                       ( [0] for remove all, 1 for keep a3m, 2 for keep all ) \n\n");
	fprintf(stderr,"-b :                   If specified, then will NOT break during iteration. \n");
	fprintf(stderr,"                       (by default, will break during iteration) \n\n");
	fprintf(stderr,"-n :                   If specified, then apply psiblast. (default=blasgpgp) \n\n");
	fprintf(stderr,"Note: the following Utility Functions, say\n");
	fprintf(stderr,"      BLAST/psiblast, BLAST_To_A3M, AlignHits, A3M_To_PSI, MSA_To_PSSM, hh_filter \n");
	fprintf(stderr,"      should be located at the [util_root] directory. \n\n");
}

//------------ main -----------//
int main(int argc,char **argv)
{
	//------ BuildAli2 -------//
	{
		if(argc<3)
		{
			Usage();
			exit(-1);
		}
		string query_file="";
		string a3m_out="";
		string nr90="databases/NR_new/nr90";
		string nr70="databases/NR_new/nr70";
		int max_iter=5;
		int cpu_num=1;
		int kill_temp=0;
		int break_or_not=0;       //default: DON'T break
		NEW_or_OLD=0;             //default: use blasgpgp
		string tmp_root="tmp/";
		string util_root="util/";

		//command-line arguments process
		extern char* optarg;
		char c = 0;
		while ((c = getopt(argc, argv, "i:o:d:D:u:m:c:r:k:bn")) != EOF) {
			switch (c) {
			case 'i':
				query_file = optarg;
				break;
			case 'o':
				a3m_out = optarg;
				break;
			case 'd':
				nr90 = optarg;
				break;
			case 'D':
				nr70 = optarg;
				break;
			case 'u':
				util_root = optarg;
				break;
			case 'm':
				max_iter = atoi(optarg);
				break;
			case 'c':
				cpu_num = atoi(optarg);
				break;
			case 'r':
				tmp_root = optarg;
				break;
			case 'k':
				kill_temp = atoi(optarg);
				break;
			case 'b':
				break_or_not = 1;
				break;
			case 'n':
				NEW_or_OLD = 1;
				break;
			default:
				Usage();
				exit(-1);
			}
		}

		//check arguments
		if(query_file=="")
		{
			fprintf(stderr,"query_file should be specified. \n");
			exit(-1);
		}
		if(max_iter<=0 || max_iter>10)
		{
			fprintf(stderr,"max_iter %d is invalid. Should be from 1 to 10. \n",max_iter);
			exit(-1);
		}
		if(cpu_num<=0)
		{
			fprintf(stderr,"cpu_num %d is invalid. Should be positive. \n",cpu_num);
			exit(-1);
		}
		if(kill_temp<0 || kill_temp>2)
		{
			fprintf(stderr,"remove %d is invalid. Should be 0,1,2. \n",kill_temp);
			exit(-1);
		}

		//---- create out_dir in the begining ---//
		int retv;
		char command[30000];
		sprintf(command,"mkdir -p %s",tmp_root.c_str());
		retv=system(command);

		//check outname
		string query_name;
		getBaseName(query_file,query_name,'/','.');
		if(a3m_out=="")a3m_out=query_name+".a3m";

		//process
		int cur_iter;
		cur_iter=BuildAli2_Main(max_iter,cpu_num,nr90,nr70,query_file,tmp_root,util_root,break_or_not);
		if(cur_iter<=0)exit(-1);

		//final hh_filter
		if(NEW_or_OLD==1)
		{
			sprintf(command,"%s/hh_filter_old -i %s/%s.a3m -o %s",
				util_root.c_str(),tmp_root.c_str(),query_name.c_str(),a3m_out.c_str());
		}
		else
		{
			sprintf(command,"%s/hh_filter_old -i %s/%s.a3m -o %s",
				util_root.c_str(),tmp_root.c_str(),query_name.c_str(),a3m_out.c_str());
		}
		retv=system(command);
		if(retv!=0)
		{
			fprintf(stderr,"seqfile %s failed at hh_filter !! \n",query_file.c_str());
			exit(-1);
		}

		//remove temporary
		if(kill_temp==0) //-> kill all temp
		{
			sprintf(command,"rm -f %s/%s.*",tmp_root.c_str(),query_name.c_str());
			retv=system(command);
		}
		else if(kill_temp==1) //-> retain iteration a3m
		{
			//-> output cur_iter a3m
			for(int i=1;i<=cur_iter;i++)
			{
				if(NEW_or_OLD==1)
				{
					sprintf(command,"%s/hh_filter_old -i %s/%s.a3m.%d -o %s.%d",
						util_root.c_str(),tmp_root.c_str(),query_name.c_str(),i,a3m_out.c_str(),i);
				}
				else
				{
					sprintf(command,"%s/hh_filter_old -i %s/%s.a3m.%d -o %s.%d",
						util_root.c_str(),tmp_root.c_str(),query_name.c_str(),i,a3m_out.c_str(),i);
				}
				retv=system(command);
				if(retv!=0)
				{
					fprintf(stderr,"seqfile %s failed at hh_filter cur_iter %d !! \n",query_file.c_str(),i);
					exit(-1);
				}
			}
			//-> remove temporary
			sprintf(command,"rm -f %s/%s.*",tmp_root.c_str(),query_name.c_str());
			retv=system(command);
		}

		//exit
		exit(0);
	}
}

