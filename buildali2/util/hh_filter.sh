#!/bin/bash
if [ $# -ne 4 ]
then
	echo "Usage: ./hhfilter_wrapper <input_a3m> <output_a3m> <home_root> <old_or_new>"
	echo "[note]: if old_or_new is set to 0, then use hh_filter_old; otherwise use hh_filter_new"
	exit
fi 

#----- input files ------#
IN_A3M=$1
OUT_A3M=$2
home_root=$3
old_or_new=$4
#-> get bname
fulnam=`basename $1`
bname=${fulnam%.*}
#-> get hhfilter
if [ $old_or_new -eq 0 ]
then
	hh_filter=$home_root/hh_filter_old
else
	hh_filter=$home_root/hh_filter_new
fi
#======== we must consider LONG FILE NAME here !!! =======# start
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmpnam=${bname}_${RANDOM}_${DATE}
cp $IN_A3M /tmp/$tmpnam.a3m
$hh_filter -i /tmp/$tmpnam.a3m -o /tmp/$tmpnam.a3m_out
cp /tmp/$tmpnam.a3m_out $OUT_A3M
rm -f /tmp/$tmpnam.a3m /tmp/$tmpnam.a3m_out
#======== we must consider LONG FILE NAME here !!! =======# end
