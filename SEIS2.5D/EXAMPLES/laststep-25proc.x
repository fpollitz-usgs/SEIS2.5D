#!/bin/sh
cat ../Work1/vertp-j ../Work2/vertp-j ../Work3/vertp-j ../Work4/vertp-j ../Work5/vertp-j ../Work6/vertp-j ../Work7/vertp-j \
../Work8/vertp-j ../Work9/vertp-j ../Work10/vertp-j ../Work11/vertp-j ../Work12/vertp-j ../Work13/vertp-j ../Work14/vertp-j ../Work15/vertp-j ../Work16/vertp-j ../Work17/vertp-j \
../Work18/vertp-j ../Work19/vertp-j ../Work20/vertp-j ../Work21/vertp-j ../Work22/vertp-j ../Work23/vertp-j \
../Work24/vertp-j ../Work25/vertp-j > /nas/users/fpollitz/scratch2/vertp-j
cat ../Work1/rec-j-m1 ../Work2/rec-j-m2 ../Work3/rec-j-m3 ../Work4/rec-j-m4 ../Work5/rec-j-m5 ../Work6/rec-j-m6 ../Work7/rec-j-m7 \
../Work8/rec-j-m8 ../Work9/rec-j-m9 ../Work10/rec-j-m10 ../Work11/rec-j-m11 ../Work12/rec-j-m12 ../Work13/rec-j-m13 ../Work14/rec-j-m14 ../Work15/rec-j-m15 ../Work16/rec-j-m16 ../Work17/rec-j-m17 \
../Work18/rec-j-m18 ../Work19/rec-j-m19 ../Work20/rec-j-m20 ../Work21/rec-j-m21 ../Work22/rec-j-m22 ../Work23/rec-j-m23 \
../Work24/rec-j-m24 ../Work25/rec-j-m25 > /nas/users/fpollitz/scratch2/rec-j-m
cat ../Work1/displ-j-m1 ../Work2/displ-j-m2 ../Work3/displ-j-m3 ../Work4/displ-j-m4 ../Work5/displ-j-m5 ../Work6/displ-j-m6 ../Work7/displ-j-m7 \
../Work8/displ-j-m8 ../Work9/displ-j-m9 ../Work10/displ-j-m10 ../Work11/displ-j-m11 ../Work12/displ-j-m12 ../Work13/displ-j-m13 ../Work14/displ-j-m14 ../Work15/displ-j-m15 ../Work16/displ-j-m16 ../Work17/displ-j-m17 \
../Work18/displ-j-m18 ../Work19/displ-j-m19 ../Work20/displ-j-m20 ../Work21/displ-j-m21 ../Work22/displ-j-m22 ../Work23/displ-j-m23 \
../Work24/displ-j-m24 ../Work25/displ-j-m25 > /nas/users/fpollitz/scratch2/displ-j-m

\rm vertp-j rec-j-m displ-j-m
ln -sf /nas/users/fpollitz/scratch2/vertp-j vertp-j
ln -sf /nas/users/fpollitz/scratch2/rec-j-m rec-j-m
ln -sf /nas/users/fpollitz/scratch2/displ-j-m displ-j-m

./seis2pt5d << ! > /dev/null
2
not-needed
not-needed
not-needed
0. 1.0000
!

\mv seis2pt5d.outxyz_vertp seis2pt5d.outxyz_vertp-AAA
\mv seis2pt5d.outstrains_vertp seis2pt5d.outstrains_vertp-AAA
\mv seis2pt5d.outCURL+DIVxyz_vertp seis2pt5d.outCURL+DIVxyz_vertp-AAA
\mv seis2pt5d.outxyz_rec seis2pt5d.outxyz_rec-AAA
\rm seis2pt5d.outxyz
