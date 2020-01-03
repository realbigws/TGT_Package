
# make the subdir first
CXX = g++
FLAG = -O3 
SHELL :=/bin/bash

all:
	#-> util
	$(CXX) util_src/AlignHits_new.cpp $(FLAG) -o AlignHits_new
	$(CXX) util_src/AlignHits_old.cpp $(FLAG) -o AlignHits_old
	$(CXX) util_src/BLAST_To_A3M_new.cpp $(FLAG) -o BLAST_To_A3M_new
	$(CXX) util_src/BLAST_To_A3M_old.cpp $(FLAG) -o BLAST_To_A3M_old
	mv AlignHits_new AlignHits_old BLAST_To_A3M_new BLAST_To_A3M_old ../util
	#-> BuildAli2
	$(CXX) BuildAli2.cpp $(FLAG) -lrt -o BuildAli2
	mv BuildAli2 ../


hhprofile:
	cd hhprofile_src&&make
	mv hhprofile_src/hh_filter hhprofile_src/hh_make ../util

BuildAli2:
	$(CXX) BuildAli2.cpp $(FLAG) -lrt -o BuildAli2
	mv BuildAli2 ../

# clean all the temp files
clean:

