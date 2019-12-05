#!/bin/bash

# ----- usage ------ #
usage()
{
	echo "AddiA3M_Gen v1.01 [Oct-01-2019] "
	echo "    Generate additional A3M from a given A3M and a sequence_db. "
	echo ""
	echo "USAGE:  ./AddiA3M_Gen.sh <-i input_fasta> <-I input_a3m> [-d database] [-o out_root] [-c CPU_num] "
	echo "                         [-m merge] [-e evalue] [-s seqid] [-K remove_tmp] [-H home] "
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i input_fasta  : Query protein sequence in FASTA format. "
	echo ""
	echo "-I input_a3m    : Initial multiple sequence alignment in A3M format. "
	echo ""
	echo "***** optional arguments *****"
	echo "-d database     : The selected database for sequence search. [default = metaclust50] "
	echo "                  The sequence database shall be put in databases/ with '.fasta' as the suffix. "
	echo ""
	echo "-o out_root     : Default output would the current directory. [default = './\${input_name}_AddiA3M'] "
	echo ""
	echo "-c CPU_num      : Number of processors. [default = 4] "
	echo ""
	echo "-m merge        : Merge the additional A3M (1) into the input A3M or not (0). [default = 0] "
	echo "                  Set -1 to re-align the input A3M with the additional A3M. "
	echo ""
	echo "-e evalue       : E-value cutoff for the selected package. [default = 0.001] "
	echo ""
	echo "-s seqid        : Remove the redundancy in the additional A3M. [default = 0.90] "
	echo "                  Set 0 to disable this module. If set, this value MUST > 0.65 "
	echo ""
	echo "-K remove_tmp   : Remove temporary folder or not. [default = 1 to remove] "
	echo ""
	echo "-H home         : home directory of MetaGenome_Gen."
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
in_seq=""
in_a3m=""

#-> optional arguments
data_db=metaclust50  #-> can be uniref90, or other sequence databases, with suffix ".fasta"
out_root=""          #-> output to current directory
cpu=4                #-> use 4 CPUs
merge=0              #-> default: 0 for NOT merge
e_value=0.001        #-> default: 0.001
seqid=0.90           #-> default: 0.90
#-> others
kill_tmp=1           #-> default: kill temporary root
home=`dirname $0`    #-> home directory


#-> parse arguments
while getopts ":i:I:d:o:c:m:e:s:K:H:" opt;
do
	case $opt in
	#-> required arguments
	i)
		in_seq=$OPTARG
		;;
	I)
		in_a3m=$OPTARG
		;;
	#-> optional arguments
	d)
		data_db=$OPTARG
		;;
	o)
		out_root=$OPTARG
		;;
	c)
		cpu=$OPTARG
		;;
	m)
		merge=$OPTARG
		;;
	e)
		e_value=$OPTARG
		;;
	s)
		seqid=$OPTARG
		;;
	K)
		kill_tmp=$OPTARG
		;;
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

# ------ check input sequence ------#
if [ ! -s "$in_seq" ]
then
	echo "in_seq $in_seq not found !!" >&2
	exit 1
fi
in_seq=`readlink -f $in_seq`
if [ ! -s "$in_a3m" ]
then
	echo "in_a3m $in_a3m not found !!" >&2
	exit 1
fi
in_a3m=`readlink -f $in_a3m`
#-> get input name
fulnam=`basename $in_seq`
relnam=${fulnam%.*}

# ------ check output directory ------#
if [ "$out_root" == "" ]
then
	out_root=${relnam}_AddiA3M
fi
mkdir -p $out_root
out_root=`readlink -f $out_root`

# ------ check SeqID ---------#
thres=0.65
if [ "$seqid" != "0" ]
then
	retv=`echo $seqid | awk '{if($0<a){print 0}else{print 1}}' a=$thres`
	if [ $retv -eq 0 ]
	then
		echo "seqid $seqid MUST > 0.65 !!">&2
		exit 1
	fi
fi


#--------------------------------------------------#
##### ===== Part 1: AddiA3M process ====== #########
#--------------------------------------------------#


#---- step 0: create temporary folder ---#
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp_root="${out_root}/TMP_AddiA3M_${relnam}_${RANDOM}_${DATE}"
mkdir -p $tmp_root

#---- step 1: convert A3M to A2M -----#
$home/util/A3M_To_A2M $in_a3m $tmp_root/$relnam.a2m 0

#---- step 2: build HMM ----#
$home/jackhmm/hmmbuild --symfrac 0 $tmp_root/$relnam.hmm $tmp_root/$relnam.a2m

#---- step 3: search against sequence database  ---#
database=$home/databases/${data_db}.fasta
$home/jackhmm/hmmsearch -A $tmp_root/$relnam.sto -o $tmp_root/$relnam.out --cpu $cpu --noali --notextw \
	-E $e_value --domE $e_value --incE $e_value --incdomE $e_value \
	--tblout $tmp_root/$relnam.tblout --domtblout $tmp_root/$relnam.domtblout $tmp_root/$relnam.hmm $database

#---- step 4: process alinged sequences ---#
grep -v "^#" $tmp_root/$relnam.sto | awk '{if(NF==2){print $0}}' > $tmp_root/$relnam.sto2
$home/util/reformat.pl sto a2m $tmp_root/$relnam.sto2 $tmp_root/$relnam.a2m
$home/util/MSA_To_SEQ $tmp_root/$relnam.a2m $tmp_root/$relnam.msa_addi
if [ $merge -eq -1 ]
then
	$home/util/MSA_To_SEQ $in_a3m $tmp_root/$relnam.msa_orig
	cat $tmp_root/$relnam.msa_orig $tmp_root/$relnam.msa_addi > $tmp_root/$relnam.msa
else
	cp $tmp_root/$relnam.msa_addi $tmp_root/$relnam.msa
fi
cat $in_seq $tmp_root/$relnam.msa > $tmp_root/$relnam.msa_seq

#---- step 4.5: remove redundancy ----#
if [ "$seqid" != "0" ]
then
	$home/jackhmm/cd-hit -i $tmp_root/$relnam.msa_seq -o $tmp_root/$relnam.msa_nored -c $seqid -s 0.9
else
	cp $tmp_root/$relnam.msa_seq $tmp_root/$relnam.msa_nored
fi

#---- step 5: re-run hmmsearch ----#
$home/jackhmm/hmmsearch -A $tmp_root/$relnam.sto_ii -o $tmp_root/$relnam.out_ii --cpu $cpu --noali --notextw \
	-E $e_value --domE $e_value --incE $e_value --incdomE $e_value \
	--tblout $tmp_root/$relnam.tblout_ii --domtblout $tmp_root/$relnam.domtblout_ii $tmp_root/$relnam.hmm $tmp_root/$relnam.msa_nored

#---- step 6: put self to first ---#
$home/util/proc_sto.sh $tmp_root/$relnam.sto_ii $relnam $tmp_root/$relnam.sto_fin

#---- step 7: generate final A3M ---#
$home/util/reformat.pl sto a3m $tmp_root/$relnam.sto_fin $tmp_root/$relnam.a3m -M first
$home/util/A3M_Seq_Refine $tmp_root/$relnam.a3m $in_seq $tmp_root/$relnam.a3m_fin

#---- step 8: merge or not ----#
if [ $merge -eq 1 ]
then
	sed '1,2d' $tmp_root/$relnam.a3m_fin > $tmp_root/$relnam.a3m_fin_noquery
	cat $in_a3m $tmp_root/$relnam.a3m_fin_noquery > $out_root/$relnam.a3m
else
	cp $tmp_root/$relnam.a3m_fin $out_root/$relnam.a3m
fi


#------- remove and exit -----#
if [ $kill_tmp -eq 1 ]
then
	rm -rf $tmp_root
else
	rm -rf $out_root/"TMP_AddiA3M_"${relnam}
	mv $tmp_root $out_root/"TMP_AddiA3M_"${relnam}
fi
cp $in_seq $out_root/$relnam.fasta_raw
cp $in_a3m $out_root/$relnam.a3m_orig

# ========= exit 0 =========== #
exit 0

