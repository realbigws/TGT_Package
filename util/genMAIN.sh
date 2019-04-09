#!/bin/bash
if [ $# -ne 4 ]
then
        echo "Usage: ./genMAIN.sh <input_a3m> <out_dir> <home_root> <cut_num=20> "
        exit
fi

#----- input files ------#
INPUTA3M=$1
DESTDIR=$2
home_root=$3
cut_num=$4
fulnam=`basename $1`
bname=${fulnam%.*}

#------ directory ------#
A3M_To_PSI=$home_root/util/A3M_To_PSI
MSA_To_PSSM=$home_root/util/MSA_To_PSSM
HHMAKE=$home_root/util/hhmake

#----- generate PSP and MTX -----#
$A3M_To_PSI $INPUTA3M $DESTDIR/$bname.psi_tmp
grep -v "ss_pred\|ss_conf\|ss_dssp" $DESTDIR/$bname.psi_tmp > $DESTDIR/$bname.psi
$MSA_To_PSSM -i $DESTDIR/$bname.psi -o $DESTDIR/$bname.psp -m $DESTDIR/$bname.mtx -c $cut_num  #-> by default, cut_num shold be set to 20
$HHMAKE -i $INPUTA3M -o $DESTDIR/$bname.hhm
rm -f $DESTDIR/$bname.psi_tmp

