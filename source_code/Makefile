

# make the subdir first
SHELL :=/bin/bash
CXX = g++
FLAG = -O3

all: meff tgt util


meff:
	#---- meff_app ------ #
	cd Meff_App&&make
	mv meff_cdhit meff_filter self_filter ../util

tgt:
	# --- A3M_To_TGT ---- #
	cd MSA_To_PSSM_cdhit_src&&make
	mv MSA_To_PSSM_cdhit_src/MSA_To_PSSM ../util
	$(CXX) A3M_To_TGT.cpp $(FLAG) -o A3M_To_TGT
	mv A3M_To_TGT ../

util:
	# ---- util ---- #
	cd util_src&&make
	#-> for A3M_To_TGT
	mv Verify_FASTA A3M_ReFormat A3M_To_A2M A3M_To_PSI SS8_To_SS3 ACC_ReFormat ../util
	mv MSA_To_SEQ A3M_Seq_Refine MSA_CovFilter ../util 

# clean all the temp files
clean:

