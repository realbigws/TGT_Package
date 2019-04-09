#!/bin/bash

#------------ install jackhmmer ----------------#
if [ ! -s jackhmmer ]
then
	echo "2.1 download hmmer"
	wget -q http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2.tar.gz
	tar xzf hmmer-3.1b2.tar.gz
	echo "2.2 compile hmmer"
	cd hmmer-3.1b2
		./configure 1> ws1 2> ws2
		make 1> ws1 2> ws2
		rm -f ws1 ws2
		cp src/jackhmmer ../
	cd ../
	#-> remove
	rm -f hmmer-3.1b2.tar.gz
	rm -rf hmmer-3.1b2
else
	echo "jackhmmer has already been installed in the local directory"
fi

