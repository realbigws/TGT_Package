#!/bin/bash
if [ $# -ne 3 ]
then
        echo "Usage: ./PSIPRED <input_mtx> <out_dir> <home_root> "
        exit
fi

# $1 the sequence name with suffix .seq
# $2 out directory

INPUTMTX=$1
DESTDIR=$2
home_root=$3
PSIPREDDIR=$home_root/util/PSIPRED



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
rootname="TMP_PSIPRED_${RANDOM}_${DATE}_${bname}"
#----- we must create a unique name here !!!! ------- ## over


if [ ! -f $DESTDIR/$rootname.mtx ] ; then
    cp $INPUTMTX $DESTDIR/$rootname.mtx
fi

echo Pass1 ....
$PSIPREDDIR/bin/psipred $DESTDIR/$rootname.mtx $PSIPREDDIR/data/weights.dat $PSIPREDDIR/data/weights.dat2 $PSIPREDDIR/data/weights.dat3 > $DESTDIR/$rootname.ss

echo Pass2
$PSIPREDDIR/bin/psipass2 $PSIPREDDIR/data/weights_p2.dat 1 0.98 1.09 $DESTDIR/$rootname.ss2 $DESTDIR/$rootname.ss > $DESTDIR/$rootname.horiz

echo "Final output files:" $rootname.ss2 $rootname.horiz $rootname.ss
mv $DESTDIR/$rootname.ss $DESTDIR/$bname.ss
mv $DESTDIR/$rootname.ss2 $DESTDIR/$bname.ss2
mv $DESTDIR/$rootname.horiz $DESTDIR/$bname.horiz


#remove temporary files
echo Cleaning up ....
rm -f $DESTDIR/$rootname.mtx

