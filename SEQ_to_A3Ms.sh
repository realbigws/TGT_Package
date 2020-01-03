#!/bin/bash

#---------------------------------------------------------#
##### ===== All functions are defined here ====== #########
#---------------------------------------------------------#


# ----- usage ------ #
usage()
{
	echo "SEQ_to_A3Ms.sh v1.00 [Dec-30-2019] "
	echo "    Generate a large set of A3Ms from a given FASTA sequence by certain strategies. "
	echo "    [note]: the strategy string shall be 'Evalue:Iteration:AddiMeff|...|...', set 'null' to omit. "
	echo ""
	echo "USAGE:  ./SEQ_to_A3Ms.sh <-i input_fasta> [-x UC_strategy] [-y UR_strategy] [-z NR_strategy] "
	echo "               [-o out_root] [-c CPU_num] [-X UC_database] [-Y UR_database] [-z NR_database] "
	echo "               [-M min_cut] [-N max_num] [-V addi_eval] [-D addi_db] [-K remove_tmp] [-f force] [-H home] "
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i input_fasta    : The input sequence in FASTA format. "
	echo ""
	echo "***** strategy arguments *****"
	echo "#--| strategy string"
	echo "-x UC_strategy    : The strategy string for HHblits against UniClust30.  [default = '1e-3:3:6|1:3:6'] "
	echo ""
	echo "-y UR_strategy    : The strategy string for JackHMMER against UniRef90.  [default = '1e-3:3:6|1e-5:3:6'] "
	echo ""
	echo "-z NR_strategy    : The strategy string for BuildAli2 against NR90+NR70. [default = '1e-3:5:6|1:5:6'] "
	echo ""
	echo "#--| relevant database"
	echo "-X UC_database    : The UniClust30 database for HHblits.  [default = uniclust30 ] "
	echo ""
	echo "-Y UR_database    : The UniRef90 database for JackHMMER.  [default = uniref90 ] "
	echo ""
	echo "-Z NR_database    : The NR90+NR70 database for BuildAli2. [default = nr_databases ] "
	echo ""
	echo "***** optional arguments *****"
	echo "#--| general options"
	echo "-o out_root       : Default output would the current directory. [default = './\${input_name}_MultA3Ms'] "
	echo ""
	echo "-c CPU_num        : CPU number. [default = 4] "
	echo ""
	echo "#--| filter strategy"
	echo "-M min_cut        : Minimal coverage of sequences in the generated MSA. [default = -1] "
	echo "                    -1 indicates that we DON'T perform any filtering. Please set from 50 to 70. "
	echo ""
	echo "-N max_num        : Maximal number of sequences in the generated MSA. [default = 20000] "
	echo "                    -1 indicates that we DON'T perform any filtering. By default, we set 20000 here. "
	echo ""
	echo "#--| additional A3M"
	echo "-V addi_eval      : run additional A3M with a given e-value. [default = 0.001] "
	echo ""
	echo "-D addi_db        : run additional A3M using a given database. [default = metaclust50] "
	echo ""
	echo "#--| other options"
	echo "-K remove_tmp     : Remove temporary folder or not. [default = 1 to remove] "
	echo ""
	echo "-f force          : If specificied, then FORCE overwrite existing files. [default = 0 NOT to] "
	echo ""
	echo "***** home relevant roots ******"
	echo "-H home           : home directory of SEQ_to_A3Ms.sh "
	echo "                    [default = `dirname $0`]"
	echo ""
	exit 1
}



#-------------------------------------------------------------#
##### ===== get pwd and check HomeDirectory ====== #########
#-------------------------------------------------------------#

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
input_fasta=""

#-> search strategy
#--| strategy string
UC_strategy="1e-3:3:6|1:3:6"
UR_strategy="1e-3:3:6|1e-5:3:6"
NR_strategy="1e-3:5:6|1:5:6"
#--| relevant database
UC_database="uniclust30"
UR_database="uniref90"
NR_database="nr_databases"

#-> optional arguments
#--| general options
out_root=""            #-> output root 
CPU_num=4              #-> default CPU number is 4
#--| filter strategies
min_cut=-1             #-> default is -1. If set, then run Cov_Filter
max_num=20000          #-> default is 20000. If set, then run Meff_Filter
#--| additional a3m
addi_eval="1e-3"       #-> default is 0.001
addi_db="metaclust50"  #-> can be ANY sequence database in plain text format
#--| others
kill_tmp=1             #-> default: kill temporary root
force=0                #-> default: NOT force overwrite

#-> home relevant
home=`dirname $0`      #-> home directory


#-> parse arguments
while getopts "i:x:y:z:X:Y:Z:o:c:M:N:V:D:K:f:H:" opt;
do
	case $opt in
	#-> required arguments
	i)
		input_fasta=$OPTARG
		;;
	#-> strategy arguments
	#--| strategy string
	x)
		UC_strategy=$OPTARG
		;;
	y)
		UR_strategy=$OPTARG
		;;
	z)
		NR_strategy=$OPTARG
		;;
	#--| relevant database
	X)
		UC_database=$OPTARG
		;;
	Y)
		UR_database=$OPTARG
		;;
	Z)
		NR_database=$OPTARG
		;;
	#-> optional arguments
	#--| general options
	o)
		out_root=$OPTARG
		;;
	c)
		CPU_num=$OPTARG
		;;
	#--| filter strategies
	M)
		min_cut=$OPTARG
		;;
	N)
		max_num=$OPTARG
		;;
	#--| additional a3m
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
	#-> help
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

#------- check input_fasta  -----------#
if [ -z "$input_fasta" ]
then
	echo "input input_fasta is null !!" >&2
	exit 1
fi
if [ ! -s "$input_fasta" ]
then
	echo "input_fasta $input_fasta not found !!" >&2
	exit 1
fi
input_fasta=`readlink -f $input_fasta`
#-> get query_name
fulnam=`basename $input_fasta`
relnam=${fulnam%.*}

# ------ check output directory -------- #
if [ "$out_root" == "" ]
then
	out_root=$curdir/${relnam}_MultA3Ms
fi
mkdir -p $out_root
out_root=`readlink -f $out_root`




#------------------------------------------------------#
##### ===== Part 1: Generate multi A3Ms ====== #########
#------------------------------------------------------#

#--------------init step: run strategies --------------------#
echo "#================ init step: run strategies to generate A3Ms =====================#"
for str in UC_strategy UR_strategy NR_strategy
do
	#----- initialization -----#
	st_code=${str:0:2}
	if [ "$st_code" == "UC" ]
	then
		method=hhsuite2
		strategies=$UC_strategy
		database=$UC_database
		echo "#------------------- step 1: run $method on $database with strategy $UC_strategy -------------------- #"
	elif [ "$st_code" == "UR" ]
	then
		method=jackhmm
		strategies=$UR_strategy
		database=$UR_database
		echo "#------------------- step 2: run $method on $database with strategy $UR_strategy -------------------- #"
	else
		method=buildali2
		strategies=$NR_strategy
		database=$NR_database
		echo "#------------------- step 3: run $method on $database with strategy $NR_strategy -------------------- #"
	fi
	#------- check null --------#
	if [ "$strategies" != "null" ] && [ "$strategies" != "" ]
	then
		st_num=`echo $strategies | awk -F"|" '{print NF}'`
		#-> for each strategy
		for ((i=1;i<=$st_num;i++))
		do
			#--| get strategy block
			strategy=`echo $strategies | cut -d '|' -f $i`
			#--| get each item
			e_value=`echo $strategy   | cut -d ':' -f 1`
			iteration=`echo $strategy | cut -d ':' -f 2`
			addi_a3m=`echo $strategy  | cut -d ':' -f 3`
			#--| check final output
			outn=${relnam}_${st_code}_${e_value}_${iteration}_${addi_a3m}
			if [ -s "$out_root/$outn.a3m" ]
			then
				continue
			fi
			#--| create output
			outr=$out_root/$outn
			mkdir -p $outr
			#--| run job
			if [ ! -s $outr/$relnam.a3m ] || [ $force -eq 1 ]
			then
				#---- submit jobs in parallel
				$home/A3M_TGT_Gen.sh -i $input_fasta -h $method -d $database -c $CPU_num -e $e_value -n $iteration -A $addi_a3m -C -1 \
					-o $outr -M $min_cut -N $max_num -V $addi_eval -D $addi_db -K $kill_tmp -f $force &
				#---- submit jobs in consideration of their strategy types
				#if [ "$st_code" == "NR" ]   #-> for NR, as the memory issue, we HAVE TO run jobs sequentially
				#then
				#	$home/A3M_TGT_Gen.sh -i $input_fasta -h $method -d $database -c $CPU_num -e $e_value -n $iteration -A $addi_a3m -C -1 \
				#		-o $outr -M $min_cut -N $max_num -V $addi_eval -D $addi_db -K $kill_tmp -f $force &
				#else                        #-> for UC and UR, we MAY run jobs in parallel
				#	$home/A3M_TGT_Gen.sh -i $input_fasta -h $method -d $database -c $CPU_num -e $e_value -n $iteration -A $addi_a3m -C -1 \
				#		-o $outr -M $min_cut -N $max_num -V $addi_eval -D $addi_db -K $kill_tmp -f $force &
				#fi
			fi
		done
		wait
	fi
done

#-------------- final collections -----------------#
echo "#============= final step: collect ALL successfully generated A3Ms ============#"
for str in UC_strategy UR_strategy NR_strategy
do
	#-> initialization
	st_code=${str:0:2}
	strategies=`eval echo \\$$str`
	#-> check null
	if [ "$strategies" != "null" ] && [ "$strategies" != "" ]
	then
		st_num=`echo $strategies | awk -F"|" '{print NF}'`
		#-> for each strategy
		for ((i=1;i<=$st_num;i++))
		do
			#--| get strategy block
			strategy=`echo $strategies | cut -d '|' -f $i`
			#--| get each item
			e_value=`echo $strategy   | cut -d ':' -f 1`
			iteration=`echo $strategy | cut -d ':' -f 2`
			addi_a3m=`echo $strategy  | cut -d ':' -f 3`
			#--| get output file
			outn=${relnam}_${st_code}_${e_value}_${iteration}_${addi_a3m}
			outr=$out_root/$outn
			if [ -s $outr/$relnam.a3m ]
			then
				mv $outr/$relnam.a3m $out_root/$outn.a3m
			fi
			#--| remove temporary folder 
			if [ $kill_tmp -eq 1 ]
			then
				rm -rf $outr
			fi
		done
	fi
done

#======================= exit ====================#
exit 0

