#!/bin/bash
if [ $# -ne 3 ]
then
        echo "Usage: ./genACC <tgt_name> <tmp_root> <home_root>"
        exit
fi

# ---- get arguments ----#
tgt_name=$1
tmp_root=$2
home_root=$3

# ---- process -----#
ACCPred=$home_root/util/ACC_Predict/acc_pred
$ACCPred $tmp_root/$tgt_name.hhm $tmp_root/$tgt_name.ss2 $tmp_root/$tgt_name.ss8 $home_root/util/ACC_Predict/model.accpred $tmp_root/$tgt_name.acc

