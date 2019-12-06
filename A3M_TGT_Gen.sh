#!/bin/bash


# ----- usage ------ #
usage()
{
	echo "A3M_TGT_Gen v1.00 [Dec-06-2019] "
	echo "    Generate A3M and TGT file from a given sequence in FASTA format. "
	echo ""
	echo "USAGE:  ./A3M_TGT_Gen.sh <-i input_fasta> [-o out_root] [-c CPU_num] [-m memory] [-h package] [-d database] "
	echo "                         [-n iteration] [-e evalue] [-E neff] [-C coverage] [-M min_cut] [-N max_num] "
	echo "                         [-A addi_meff] [-V addi_eval] [-D addi_db] [-K remove_tmp] [-f force] [-H home] "
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i input_fasta  : Query protein sequence in FASTA format. "
	echo ""
	echo "***** optional arguments *****"
	echo "#--| misc parameters"
	echo "-o out_root     : Default output would the current directory. [default = './\${input_name}_A3MTGT'] "
	echo ""
	echo "-c CPU_num      : Number of processors. [default = 4] "
	echo ""
	echo "-m memory       : Maximal allowed memory (for hhsuite2 or hhsuite3 only). [default = 3.0 (G)] "
	echo ""
	echo "#--| search engine and database"
	echo "-h package      : The selected package to generate A3M file. [default = hhsuite2] "
	echo "                  users may use other packages: hhsuite3, jackhmm, or buildali2."
	echo ""
	echo "-d database     : The selected database for sequence search. [default = uniprot20_2016_02] "
	echo "                  users may use other uniprot20 databases to run hhsuite2 or hhsuite3,"
	echo "                  or use uniref90 for jackhmm, and nr_databases for buildali2."
	echo ""
	echo "#--| search strategy"
	echo "-n iteration    : Maximal iteration to run the seleced package. [default = 2] "
	echo ""
	echo "-e evalue       : E-value cutoff for the selected package. [default = 0.001] "
	echo ""
	echo "-E neff         : Neff cutoff for threading purpose (i.e., -C -2). [default = 7] (for hhsuite only) "
	echo ""
	echo "-C coverage     : Coverage for hhsuite only. [default = -2 (i.e., NOT use -cov in HHblits)] "
	echo "                  if set to -1, then automatically determine coverage value. "
	echo "                  if set to any other positive value, then use this -cov in HHblits. "
	echo ""
	echo "#--| filter strategy"
	echo "-M min_cut      : Minimal coverage of sequences in the generated MSA. [default = -1] "
	echo "                  -1 indicates that we DON'T perform any filtering. Please set from 50 to 70. "
	echo ""
	echo "-N max_num      : Maximal number of sequences in the generated MSA. [default = -1] "
	echo "                  -1 indicates that we DON'T perform any filtering. For example, set 20000 here. "
	echo ""
	echo "#--| additional A3M"
	echo "-A addi_meff    : run additional A3M only if the previous ln(meff) is lower than this. [default = -1] "
	echo "                  -1 indicates that we DON'T search for additional A3M. "
	echo ""
	echo "-V addi_eval    : run additional A3M with a given e-value. [default = 0.001] "
	echo ""
	echo "-D addi_db      : run additional A3M using a given database. [default = metaclust50] "
	echo ""
        echo "#--| other options"
	echo "-K remove_tmp   : Remove temporary folder or not. [default = 1 to remove] "
	echo ""
	echo "-f force        : If specificied, then FORCE overwrite existing files. [default = 0 NOT to] "
	echo ""
	echo "***** home relevant **********"
	echo "-H home         : home directory of TGT_Package."
	echo "                  [default = `dirname $0`] "
	echo ""
	exit 1
}

#------------------------------------------------------------#
##### ===== get pwd and check BlastSearchHome ====== #########
#------------------------------------------------------------#

#------ current directory ------#
curdir="$(pwd)"

#-------- check usage -------#
if [ $# -lt 1 ];
then
        usage
fi



#---------------------------------------------------------#
##### ===== All arguments are defined here ====== #########
#---------------------------------------------------------#


# ----- get arguments ----- #
#-> required arguments
query_seq=""
#-> optional arguments
#--| misc parameters
out_root=""     #-> output to current directory
cpu_num=4       #-> use 4 CPUs
memory="3.0"    #-> use 3.0G memory
#--| search engine and database
hhsuite=hhsuite2              #-> can be hhsuite2, hhsuite3, jackhmm, and buildali2
uniprot20=uniprot20           #-> can be uniprot20, uniclust30, uniref90, and NR_New
#--| search strategies
iteration=2     #-> default is 2 iterations, for threading purpose
e_value=0.001   #-> default is 0.001, for threading purpose
neffmax=7       #-> default is 7, for threading purpose
coverage=-2     #-> automatic determine the coverage on basis of input sequence length (i.e., for threading)
#--| filter strategies
min_cut=-1      #-> default is -1. If set, then run Cov_Filter
max_num=-1      #-> default is -1. If set, then run Meff_Filter
#--| additional a3m
addi_meff=-1    #-> -1 means that we DON'T search additional a3m
addi_eval=0.001 #-> default is 0.001
addi_db=metaclust50           #-> can be ANY sequence database in plain text format
#--| other options
kill_tmp=1      #-> default: kill temporary root
force=0         #-> default: NOT force overwrite
#--| home relevant
home=`dirname $0`  #-> home directory


#-> parse arguments
while getopts ":i:o:c:m:h:d:n:e:E:C:M:N:A:V:D:K:f:H:" opt;
do
	case $opt in
	#-> required arguments
	i)
		input_fasta=$OPTARG
		;;
	#-> optional arguments
	#--| misc parameters
	o)
		out_root=$OPTARG
		;;
	c)
		cpu_num=$OPTARG
		;;
	m)
		memory=$OPTARG
		;;
	#--| search engine and database
	h)
		hhsuite=$OPTARG
		;;
	d)
		uniprot20=$OPTARG
		;;
	#--| search strategies
	n)
		iteration=$OPTARG
		;;
	e)
		e_value=$OPTARG
		;;
	E)
		neffmax=$OPTARG
		;;
	C)
		coverage=$OPTARG
		;;
	#--| filter strategies
	M)
		min_cut=$OPTARG
		;;
	N)
		max_num=$OPTARG
		;;
	#--| additional a3m
	A)
		addi_meff=$OPTARG
		;;
	V)
		addi_eval=$OPTARG
		;;
	D)
		addi_db=$OPTARG
		;;
	#--| others options
	K)
		kill_tmp=$OPTARG
		;;
	f)
		force=$OPTARG
		;;
	#-> home relevant
	H)
		home=$OPTARG
		;;
	#-> others
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done


#---------------------------------------------------------#
##### ===== Part 0: initial argument check ====== #########
#---------------------------------------------------------#

# ------ check home directory ---------- #
if [ ! -d "$home" ]
then
	echo "home directory $home not exist " >&2
	exit 1
fi
home=`readlink -f $home`

# ------ check input fasta ------#
if [ ! -s "$input_fasta" ]
then
	echo "input_fasta $input_fasta not found !!" >&2
	exit 1
fi
input_fasta=`readlink -f $input_fasta`
fulnam=`basename $input_fasta`
relnam=${fulnam%.*}

# ------ check output directory ------#
if [ "$out_root" == "" ]
then
	out_root=${relnam}_A3MTGT
fi
mkdir -p $out_root
out_root=`readlink -f $out_root`



#-------------------------------------------------#
##### ===== Part 1: A3MTGT process ====== #########
#-------------------------------------------------#


# --- create temporary folder --#
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp_root="${out_root}/TMP_A3MTGT_${relnam}_${RANDOM}_${DATE}"
mkdir -p $tmp_root


# ---- verify FASTA file -------- #
seq_file=$relnam.seq
$home/util/Verify_FASTA $input_fasta $out_root/$seq_file
OUT=$?
if [ $OUT -ne 0 ]
then
	echo "failed in util/Verify_FASTA $input_fasta $seq_file"
	exit 1
fi


# ----- determine coverage ---- #
if [ $coverage -eq -1 ]
then
	a=60
	b=`tail -n1 $out_root/$seq_file | wc | awk '{print int(7000/($3-1))}'`
	if [ $a -gt $b ]
	then
		coverage=$b
	else
		coverage=$a
	fi
fi


# ---- generate A3M file -------- #
a3m_file=$relnam.a3m
if [ ! -f "$out_root/$a3m_file" ] || [ $force -eq 1 ] 
then
	#---- this is the default home ----#
	HHSUITE=$home/$hhsuite
	#---- run HHpred (v2 or v3) ----------#
	if [ "$hhsuite" == "hhsuite2" ] || [ "$hhsuite" == "hhsuite3" ]
	then
		HHLIB=$HHSUITE/lib/hh
		echo "hhblits start with database $uniprot20 with evalue $e_value and iteration $iteration with cpu $cpu_num"
		if [ $coverage -eq -2 ]
		then
			echo "run HHblits with default parameter without -cov "
			$HHSUITE/bin/hhblits -i $out_root/$seq_file -cpu $cpu_num -d $home/databases/$uniprot20/$uniprot20 -o $tmp_root/$relnam.hhr \
				-oa3m $out_root/$relnam.a3m -n $iteration -e $e_value -neffmax $neffmax -maxmem $memory
		else
			echo "run HHblits with -maxfilt 500000 -diff inf -id 99 -cov $coverage"
			$HHSUITE/bin/hhblits -i $out_root/$seq_file -cpu $cpu_num -d $home/databases/$uniprot20/$uniprot20 -o $tmp_root/$relnam.hhr \
				-oa3m $out_root/$relnam.a3m -n $iteration -e $e_value -maxfilt 500000 -diff inf -id 99 -cov $coverage -maxmem $memory
		fi
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "failed in $HHSUITE/bin/hhblits -i $seq_file -cpu $cpu_num -d databases/$uniprot20/$uniprot20 -n $iteration"
			exit 1
		fi
		echo "hhblits done"
	else
		#---- run JackHmmer ---------#
		if [ "$hhsuite" == "jackhmm" ]
		then
			echo "jackhmm start with database $uniprot20 with evalue $e_value and iteration $iteration with cpu $cpu_num"
			#-> run jackhmmer
			$HHSUITE/jackhmmer -A $tmp_root/$relnam.sto -N $iteration --cpu $cpu_num -E $e_value --domE $e_value --incE $e_value --incdomE $e_value \
				--noali --notextw --tblout $tmp_root/$relnam.tblout --domtblout $tmp_root/$relnam.domtblout -o $tmp_root/$relnam.output \
				$out_root/$seq_file $home/databases/$uniprot20.fasta
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "failed in $HHSUITE/jackhmmer -A $tmp_root/$relnam.sto $out_root/$seq_file $home/databases/$uniprot20.fasta"
				exit 1
			fi
			#-> transform from sto to a3m
			$HHSUITE/util/reformat.pl $tmp_root/$relnam.sto $tmp_root/$relnam.a3m -M first
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "failed in $HHSUITE/util/reformat.pl $tmp_root/$relnam.sto $tmp_root/$relnam.a3m -M first"
				exit 1
			fi
			#-> run hhfilter or not
			if [ $coverage -ne -2 ]
			then
				$HHSUITE/util/self_filter -i $tmp_root/$relnam.a3m -o $out_root/$relnam.a3m -s 0.99 -d 0.5
			else
				cp $tmp_root/$relnam.a3m $out_root
			fi
			echo "jackhammer done"
		fi
		#---- run BLAST wrapped buildali2 ---------#
		if [ "$hhsuite" == "buildali2" ]
		then
			echo "BuildAli2 with database $uniprot20 with evalue $e_value and iteration $iteration with cpu $cpu_num"
			#-> run Sheng's BuildAli2
			$HHSUITE/BuildAli2 -i $out_root/$seq_file -o $out_root/$relnam.a3m -d $home/databases/$uniprot20/nr90 -D $home/databases/$uniprot20/nr70 \
				-u $HHSUITE/util -m $iteration -e $e_value -c $cpu_num -r $tmp_root
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "failed in $HHSUITE/BuildAli2 -i $out_root/$seq_file -o $out_root/$relnam.a3m -d $home/databases/$uniprot20/nr90"
				exit 1
			fi
			#-> replace QUERY with relnam
			sed -i "1 s/^>QUERY/>$relnam/" $out_root/$relnam.a3m
			echo "BuildAli2 done"
		fi
	fi
fi



# ========================================= post process ====================================#

# ---- search for additional a3m ---- #
if [ $addi_meff -ne -1 ]
then
	#-> calculate original Meff
	meff=`$home/util/meff_cdhit -i $out_root/$relnam.a3m`
	echo "ln(meff) of the original A3M is $meff"
	#-> judge Meff
	if [ 1 -eq "$(echo "$meff < $addi_meff" | bc)" ]
	then
		#--| additional a3m
		$home/AddiA3M_Gen.sh -i $out_root/$seq_file -I $out_root/$relnam.a3m -o $out_root/${relnam}_AddiA3M
		cp $out_root/${relnam}_AddiA3M/$relnam.a3m $out_root
		#--| more Meff
		meff_addi=`$home/util/meff_cdhit -i $out_root/$relnam.a3m`
		echo "ln(meff) of the additional A3M is $meff_addi"
	fi
fi


# ---- filter sequences in MSA --- #
#-> coverage filter
if [ $min_cut -ne -1 ]
then
	#-> calculate the number of sequences in A3M
	numLines=`grep "^>" $out_root/$relnam.a3m | wc | awk '{print $1}'`
	echo "number of lines in the original A3M before CovFilt is $numLines"
	#-> coverage filter
	cp $out_root/$relnam.a3m $out_root/$relnam.a3m_prev
	$home/util/MSA_CovFilter $out_root/$relnam.a3m_prev $out_root/$relnam.a3m $min_cut
	#-> calculate the number of sequences in A3M
	numLines=`grep "^>" $out_root/$relnam.a3m | wc | awk '{print $1}'`
	echo "number of lines in the cov_filt A3M before CovFilt is $numLines"
fi
#-> maximal number filter
if [ $max_num -ne -1 ]
then
	#-> calculate the number of sequences in A3M
	numLines=`grep "^>" $out_root/$relnam.a3m | wc | awk '{print $1}'`
	echo "number of lines in the original A3M is $numLines"
	#-> judge numLines
	if [ $numLines -gt $max_num ]
	then
		#--| filter
		$home/util/meff_filter -i $out_root/$relnam.a3m -o $out_root/$relnam.a3m_filter -n $max_num
		a3m_file=$relnam.a3m_filter
		#--| filter result
		numLines=`grep "^>" $out_root/$a3m_file | wc | awk '{print $1}'`
		echo "number of lines in the filtered A3M is $numLines"
	fi
fi

# ---- generate TGT file ------ #
tgt_file=$relnam.tgt
if [ ! -f "$out_root/$tgt_file" ] || [ $force -eq 1 ]
then
	$home/A3M_To_TGT -i $out_root/$seq_file -I $out_root/$a3m_file -o $out_root/$tgt_file -t $tmp_root -H $home 
	OUT=$?
	if [ $OUT -ne 0 ]
	then
		echo "failed in ./A3M_To_TGT -i $seq_file -I $a3m_file -o $tgt_file -t $tmp_root"
		exit 1
	fi
fi

# ---- post process ----- #
if [ $kill_tmp -eq 1 ]
then
	rm -rf $tmp_root
else
	rm -rf $out_root/"TMP_A3MTGT_"${relnam}
	mv $tmp_root $out_root/"TMP_A3MTGT_"${relnam}
fi
cp $input_fasta $out_root/$relnam.fasta_raw


# ========= exit 0 =========== #
exit 0


