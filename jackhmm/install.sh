#!/bin/bash

#------------ 1) install jackhmmer ----------------#
if [ ! -s "jackhmmer" ] || 
   [ ! -s "hmmbuild" ]  ||
   [ ! -s "hmmsearch" ]
then
	echo "1.1 download hmmer"
	wget -q http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2.tar.gz
	tar xzf hmmer-3.1b2.tar.gz
	echo "1.2 compile hmmer"
	cd hmmer-3.1b2
		./configure 2> /dev/null
		make 2> /dev/null
		cp src/jackhmmer src/hmmbuild src/hmmsearch  ../
	cd ../
	#-> remove
	rm -f hmmer-3.1b2.tar.gz
	rm -rf hmmer-3.1b2
else
	echo "jackhmmer has already been installed in the local directory"
fi

#------------- 2) install cd-hit ------------------------#
if [ ! -s "cd-hit" ]
then
	echo "2.1 download cd-hit"
	git clone https://github.com/weizhongli/cdhit
	echo "2.2 compile cd-hit"
	cd cdhit
		make 2> /dev/null
		cp cd-hit ../
	cd ../
	#-> remove
	rm -rf cdhit
else
	echo "cd-hit has already been installed in the local directory"
fi
