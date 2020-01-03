

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

all:
	# ------ A3M to TGT -------#
	$(CXX) Verify_FASTA.cpp $(FLAG) -o Verify_FASTA
	$(CXX) A3M_ReFormat.cpp $(FLAG) -o A3M_ReFormat
	$(CXX) A3M_To_A2M.cpp $(FLAG) -o A3M_To_A2M
	$(CXX) A3M_To_PSI.cpp $(FLAG) -o A3M_To_PSI
	$(CXX) SS8_To_SS3.cpp $(FLAG) -o SS8_To_SS3
	$(CXX) ACC_ReFormat.cpp $(FLAG) -o ACC_ReFormat
	# ------ A3M addi ---------#
	$(CXX) MSA_To_SEQ.cpp $(FLAG) -o MSA_To_SEQ
	$(CXX) A3M_Seq_Refine.cpp $(FLAG) -o A3M_Seq_Refine
	$(CXX) MSA_CovFilter.cpp $(FLAG) -o MSA_CovFilter
	# move to ../ directory
	mv Verify_FASTA A3M_ReFormat A3M_To_A2M A3M_To_PSI SS8_To_SS3 ACC_ReFormat ../
	mv MSA_To_SEQ A3M_Seq_Refine MSA_CovFilter ../

# clean all the temp files
clean:
