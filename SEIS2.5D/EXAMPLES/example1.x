#!/bin/sh

export SEISDIR=/home/fpollitz/fred4/CIDER2015/SEIS2.5D/
export EXDIR=$SEISDIR/MAINPROG


echo "Running forward simulation on Stephenson et al. 2017 structure"
\cp seisstruc.paramSTE-sharp seisstruc.param
\cp seis2pt5d-sphericalEX1.info seis2pt5d-spherical.info

\cp $EXDIR/seis2pt5d .
\rm seis2pt5dsource-spherical.param
ln -sf seis2pt5dsource-spherical.param17Jul2019 seis2pt5dsource-spherical.param
ln -sf topo_EX1 topofile

sed s/AAA/"Monroe_eq_aftershock"/ laststep-25proc.x > laststep.x

\cp receivers-latlondep17Jul2019.txt receivers-latlondep.txt
\rm ../Work*/rec-j-m* ../Work*/displ-j-m*

date > ../Work1/seis2pt5d-progress1.txt
date > ../Work2/seis2pt5d-progress2.txt
date > ../Work3/seis2pt5d-progress3.txt
date > ../Work4/seis2pt5d-progress4.txt
date > ../Work5/seis2pt5d-progress5.txt
date > ../Work6/seis2pt5d-progress6.txt
date > ../Work7/seis2pt5d-progress7.txt
date > ../Work8/seis2pt5d-progress8.txt
date > ../Work9/seis2pt5d-progress9.txt
date > ../Work10/seis2pt5d-progress10.txt
date > ../Work11/seis2pt5d-progress11.txt
date > ../Work12/seis2pt5d-progress12.txt
date > ../Work13/seis2pt5d-progress13.txt
date > ../Work14/seis2pt5d-progress14.txt
date > ../Work15/seis2pt5d-progress15.txt
date > ../Work16/seis2pt5d-progress16.txt
date > ../Work17/seis2pt5d-progress17.txt
date > ../Work18/seis2pt5d-progress18.txt
date > ../Work19/seis2pt5d-progress19.txt
date > ../Work20/seis2pt5d-progress20.txt
date > ../Work21/seis2pt5d-progress21.txt
date > ../Work22/seis2pt5d-progress22.txt
date > ../Work23/seis2pt5d-progress23.txt
date > ../Work24/seis2pt5d-progress24.txt
date > ../Work25/seis2pt5d-progress25.txt

cp seisstruc.param ../Work1/.
cp seisstruc.param ../Work2/.
cp seisstruc.param ../Work3/.
cp seisstruc.param ../Work4/.
cp seisstruc.param ../Work5/.
cp seisstruc.param ../Work6/.
cp seisstruc.param ../Work7/.
cp seisstruc.param ../Work8/.
cp seisstruc.param ../Work9/.
cp seisstruc.param ../Work10/.
cp seisstruc.param ../Work11/.
cp seisstruc.param ../Work12/.
cp seisstruc.param ../Work13/.
cp seisstruc.param ../Work14/.
cp seisstruc.param ../Work15/.
cp seisstruc.param ../Work16/.
cp seisstruc.param ../Work17/.
cp seisstruc.param ../Work18/.
cp seisstruc.param ../Work19/.
cp seisstruc.param ../Work20/.
cp seisstruc.param ../Work21/.
cp seisstruc.param ../Work22/.
cp seisstruc.param ../Work23/.
cp seisstruc.param ../Work24/.
cp seisstruc.param ../Work25/.

cp seis2pt5d-spherical.info ../Work1/.
cp seis2pt5d-spherical.info ../Work2/.
cp seis2pt5d-spherical.info ../Work3/.
cp seis2pt5d-spherical.info ../Work4/.
cp seis2pt5d-spherical.info ../Work5/.
cp seis2pt5d-spherical.info ../Work6/.
cp seis2pt5d-spherical.info ../Work7/.
cp seis2pt5d-spherical.info ../Work8/.
cp seis2pt5d-spherical.info ../Work9/.
cp seis2pt5d-spherical.info ../Work10/.
cp seis2pt5d-spherical.info ../Work11/.
cp seis2pt5d-spherical.info ../Work12/.
cp seis2pt5d-spherical.info ../Work13/.
cp seis2pt5d-spherical.info ../Work14/.
cp seis2pt5d-spherical.info ../Work15/.
cp seis2pt5d-spherical.info ../Work16/.
cp seis2pt5d-spherical.info ../Work17/.
cp seis2pt5d-spherical.info ../Work18/.
cp seis2pt5d-spherical.info ../Work19/.
cp seis2pt5d-spherical.info ../Work20/.
cp seis2pt5d-spherical.info ../Work21/.
cp seis2pt5d-spherical.info ../Work22/.
cp seis2pt5d-spherical.info ../Work23/.
cp seis2pt5d-spherical.info ../Work24/.
cp seis2pt5d-spherical.info ../Work25/.

cp seis2pt5dsource-spherical.param ../Work1/.
cp seis2pt5dsource-spherical.param ../Work2/.
cp seis2pt5dsource-spherical.param ../Work3/.
cp seis2pt5dsource-spherical.param ../Work4/.
cp seis2pt5dsource-spherical.param ../Work5/.
cp seis2pt5dsource-spherical.param ../Work6/.
cp seis2pt5dsource-spherical.param ../Work7/.
cp seis2pt5dsource-spherical.param ../Work8/.
cp seis2pt5dsource-spherical.param ../Work9/.
cp seis2pt5dsource-spherical.param ../Work10/.
cp seis2pt5dsource-spherical.param ../Work11/.
cp seis2pt5dsource-spherical.param ../Work12/.
cp seis2pt5dsource-spherical.param ../Work13/.
cp seis2pt5dsource-spherical.param ../Work14/.
cp seis2pt5dsource-spherical.param ../Work15/.
cp seis2pt5dsource-spherical.param ../Work16/.
cp seis2pt5dsource-spherical.param ../Work17/.
cp seis2pt5dsource-spherical.param ../Work18/.
cp seis2pt5dsource-spherical.param ../Work19/.
cp seis2pt5dsource-spherical.param ../Work20/.
cp seis2pt5dsource-spherical.param ../Work21/.
cp seis2pt5dsource-spherical.param ../Work22/.
cp seis2pt5dsource-spherical.param ../Work23/.
cp seis2pt5dsource-spherical.param ../Work24/.
cp seis2pt5dsource-spherical.param ../Work25/.

cp receivers-latlondep.txt ../Work1/.
cp receivers-latlondep.txt ../Work2/.
cp receivers-latlondep.txt ../Work3/.
cp receivers-latlondep.txt ../Work4/.
cp receivers-latlondep.txt ../Work5/.
cp receivers-latlondep.txt ../Work6/.
cp receivers-latlondep.txt ../Work7/.
cp receivers-latlondep.txt ../Work8/.
cp receivers-latlondep.txt ../Work9/.
cp receivers-latlondep.txt ../Work10/.
cp receivers-latlondep.txt ../Work11/.
cp receivers-latlondep.txt ../Work12/.
cp receivers-latlondep.txt ../Work13/.
cp receivers-latlondep.txt ../Work14/.
cp receivers-latlondep.txt ../Work15/.
cp receivers-latlondep.txt ../Work16/.
cp receivers-latlondep.txt ../Work17/.
cp receivers-latlondep.txt ../Work18/.
cp receivers-latlondep.txt ../Work19/.
cp receivers-latlondep.txt ../Work20/.
cp receivers-latlondep.txt ../Work21/.
cp receivers-latlondep.txt ../Work22/.
cp receivers-latlondep.txt ../Work23/.
cp receivers-latlondep.txt ../Work24/.
cp receivers-latlondep.txt ../Work25/.

cp topofile ../Work1/.
cp topofile ../Work2/.
cp topofile ../Work3/.
cp topofile ../Work4/.
cp topofile ../Work5/.
cp topofile ../Work6/.
cp topofile ../Work7/.
cp topofile ../Work8/.
cp topofile ../Work9/.
cp topofile ../Work10/.
cp topofile ../Work11/.
cp topofile ../Work12/.
cp topofile ../Work13/.
cp topofile ../Work14/.
cp topofile ../Work15/.
cp topofile ../Work16/.
cp topofile ../Work17/.
cp topofile ../Work18/.
cp topofile ../Work19/.
cp topofile ../Work20/.
cp topofile ../Work21/.
cp topofile ../Work22/.
cp topofile ../Work23/.
cp topofile ../Work24/.
cp topofile ../Work25/.

cp seis2pt5d ../Work1/.
cp seis2pt5d ../Work2/.
cp seis2pt5d ../Work3/.
cp seis2pt5d ../Work4/.
cp seis2pt5d ../Work5/.
cp seis2pt5d ../Work6/.
cp seis2pt5d ../Work7/.
cp seis2pt5d ../Work8/.
cp seis2pt5d ../Work9/.
cp seis2pt5d ../Work10/.
cp seis2pt5d ../Work11/.
cp seis2pt5d ../Work12/.
cp seis2pt5d ../Work13/.
cp seis2pt5d ../Work14/.
cp seis2pt5d ../Work15/.
cp seis2pt5d ../Work16/.
cp seis2pt5d ../Work17/.
cp seis2pt5d ../Work18/.
cp seis2pt5d ../Work19/.
cp seis2pt5d ../Work20/.
cp seis2pt5d ../Work21/.
cp seis2pt5d ../Work22/.
cp seis2pt5d ../Work23/.
cp seis2pt5d ../Work24/.
cp seis2pt5d ../Work25/.

ZEROETH=`qsub   -wd $SEISDIR/EXAMPLES -l nodes=1 step0.x`
echo $ZEROETH

FIRST=`qsub -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work1 -l nodes=1 step1.x`
echo $FIRST
SECOND=`qsub -W depend=afterok:$ZEROETH   -wd $SEISDIR/Work2 -l nodes=1 step2.x`
echo $SECOND
THIRD=`qsub -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work3 -l nodes=1 step3.x`
echo $THIRD
FOURTH=`qsub -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work4 -l nodes=1 step4.x`
echo $FOURTH
FIFTH=`qsub -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work5 -l nodes=1 step5.x`
echo $FIFTH
SIXTH=`qsub -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work6 -l nodes=1 step6.x`
echo $SIXTH
SEVENTH=`qsub -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work7 -l nodes=1 step7.x`
echo $SEVENTH 
EIGHTH=`qsub -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work8 -l nodes=1 step8.x`
echo $EIGHTH
NINTH=`qsub -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work9 -l nodes=1 step9.x`
echo $NINTH
TENTH=`qsub -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work10 -l nodes=1 step10.x`
echo $TENTH
ELEVENTH=`qsub -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work11 -l nodes=1 step11.x`
echo $ELEVENTH
TWELVTH=`qsub -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work12 -l nodes=1 step12.x`
echo $TWELVTH
THIRTEENTH=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work13 -l nodes=1 step13.x`
echo $THIRTEENTH
FOURTEENTH=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work14 -l nodes=1 step14.x`
echo $FOURTEENTH
FIFTEENTH=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work15 -l nodes=1 step15.x`
echo $FIFTEENTH
SIXTEENTH=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work16 -l nodes=1 step16.x`
echo $SIXTEENTH
SEVENTEENTH=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work17 -l nodes=1 step17.x`
echo $SEVENTEENTH
EIGHTEENTH=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work18 -l nodes=1 step18.x`
echo $EIGHTEENTH
NINETEENTH=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work19 -l nodes=1 step19.x`
echo $NINETEENTH
TWENTIETH=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work20 -l nodes=1 step20.x`
echo $TWENTIETH
TWENTYFIRST=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work21 -l nodes=1 step21.x`
echo $TWENTYFIRST
TWENTYSECOND=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work22 -l nodes=1 step22.x`
echo $TWENTYSECOND
TWENTYTHIRD=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work23 -l nodes=1 step23.x`
echo $TWENTYTHIRD
TWENTYFOURTH=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work24 -l nodes=1 step24.x`
echo $TWENTYFOURTH
TWENTYFIFTH=`qsub  -W depend=afterok:$ZEROETH  -wd $SEISDIR/Work25 -l nodes=1 step25.x`
echo $TWENTYFIFTH

TWENTYSIXTH=`qsub  -W depend=afterok:$FIRST:$SECOND:$THIRD:$FOURTH:$FIFTH:$SIXTH:$SEVENTH:$EIGHTH:$NINTH:$TENTH:$ELEVENTH:$TWELVTH:$THIRTEENTH:$FOURTEENTH:$FIFTEENTH:$SIXTEENTH:$SEVENTEENTH:$EIGHTEENTH:$NINETEENTH:$TWENTIETH:$TWENTYFIRST:$TWENTYSECOND:$TWENTYTHIRD:$TWENTYFOURTH:$TWENTYFIFTH -wd $SEISDIR/EXAMPLES -l nodes=1 laststep.x`
echo $TWENTYSIXTH

