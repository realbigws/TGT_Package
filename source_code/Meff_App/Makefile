
# make the subdir first
CXX = g++
FLAG = -O3 
SHELL :=/bin/bash

all:
	#-> meff_cdhit
	cd meff_cdhit_src&&make
	#-> meff_filter
	cd meff_filter_src&&make
	#-> self_filter
	cd self_filter_src&&make
	#-> move
	mv meff_cdhit_src/meff_cdhit meff_filter_src/meff_filter self_filter_src/self_filter ../

# clean all the temp files
clean:

