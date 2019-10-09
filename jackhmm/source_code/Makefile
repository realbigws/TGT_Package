

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

all:
	#-> self_filter
	cd self_filter_src&&make
	mv self_filter_src/self_filter .
	# move to ../util directory
	mv self_filter ../util

# clean all the temp files
clean:
