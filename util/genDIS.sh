#!/bin/bash
if [ $# -ne 3 ]
then
        echo "Usage: ./DISOPRED <input_mtx> <out_dir> <home_root>"
        exit
fi

# $1 the sequence name with suffix .seq
# $2 out directory

INPUTMTX=$1
DESTDIR=$2
home_root=$3
PSIPREDDIR=$home_root/util/DISOPRED


###### BY TINA, May 6, 2003
# first make PSP directory if it doesn't exist
if [ ! -d $DESTDIR ] ; then
    mkdir $DESTDIR
fi
###### BY TINA, May 6, 2003


#----- we must create a unique name here !!!! ------- ## 2018.10.22 by Sheng Wang 
fulnam=`basename $1`
bname=${fulnam%.*}
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
rootname="TMP_DISOPRED_${RANDOM}_${DATE}_${bname}"
#----- we must create a unique name here !!!! ------- ## over


if [ ! -f $DESTDIR/$rootname.mtx ] ; then
    cp $INPUTMTX $DESTDIR/$rootname.mtx
fi

echo Pass1 ...
echo Pass2 ...
$PSIPREDDIR/bin/disopred $rootname $DESTDIR/$rootname.mtx $PSIPREDDIR/data/


echo "Final output files:" $rootname.diso $rootname.horiz_d
mv $rootname.diso $DESTDIR/$bname.diso
mv $rootname.horiz_d $DESTDIR/$bname.horiz_d

#remove temporary files
echo Cleaning up ....
rm -f $DESTDIR/$rootname.mtx

