#!/bin/bash
if [ $# -ne 5 ]
then
        echo "Usage: ./genMAIN.sh <input_a3m> <out_dir> <home_root> <cpu_num> <cut_num=20> "
        exit
fi

#----- input files ------#
INPUTA3M=$1
DESTDIR=$2
home_root=$3
cpu_num=$4
cut_num=$5
fulnam=`basename $1`
bname=${fulnam%.*}

#------ directory ------#
A3M_To_PSI=$home_root/util/A3M_To_PSI
MSA_To_PSSM=$home_root/util/MSA_To_PSSM
HHMAKE=$home_root/util/hhmake

#----- generate PSP and MTX -----#
$A3M_To_PSI $INPUTA3M $DESTDIR/$bname.psi_tmp
grep -v "ss_pred\|ss_conf\|ss_dssp" $DESTDIR/$bname.psi_tmp > $DESTDIR/$bname.psi
$MSA_To_PSSM -i $DESTDIR/$bname.psi -o $DESTDIR/$bname.psp -m $DESTDIR/$bname.mtx -c $cut_num -C $cpu_num  #-> by default, cut_num shold be set to 20

#======== we must consider LONG FILE NAME here !!! =======# start
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmpnam=${bname}_${RANDOM}_${DATE}
cp $INPUTA3M /tmp/$tmpnam.a3m
$HHMAKE -i /tmp/$tmpnam.a3m -o /tmp/$tmpnam.hhm
cp /tmp/$tmpnam.hhm $DESTDIR/$bname.hhm
rm -f /tmp/$tmpnam.a3m /tmp/$tmpnam.hhm
#======== we must consider LONG FILE NAME here !!! =======# end

rm -f $DESTDIR/$bname.psi_tmp

