#!/bin/bash
if [ $# -ne 3 ]
then
        echo "Usage: ./genSS8 <tgt_name> <tmp_root> <home_root>"
        exit
fi

# ---- get arguments ----#
tgt_name=$1
tmp_root=$2
home_root=$3

# ---- process -----#
SS8_Pred=$home_root/util/SS8_Predict/bin/run_raptorx-ss8.pl
$SS8_Pred $tmp_root/$tgt_name.seq -pssm $tmp_root/$tgt_name.psp -outdir $tmp_root/ -home $home_root/

