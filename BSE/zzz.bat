#!/bin/sh

###mass z tphysf	in the shell loop !!!	150.0 0.001 10000.0

for mi in 5 10 16 25 40 60 80 100 120 150 180 200 ### in units of 1 Msol
###for mi in  0.64 0.8 0.91 2 3 4 5 6 7 8 9 10 11
do

for zi in 5 10 50 100 200 	                  ### in units of 0.0001
do

for ti in 1.2			                  ### in units of 10000 Myr
do

mmm=`echo $mi*1.0 | bc -l`

zzz=`echo $zi*0.0001 | bc -l`

ttt=`echo $ti*10000.0 | bc -l`

echo M = ${mmm} [Msol] Z = ${zzz} [Zsol = 0.02] T_end = ${ttt} [Myr]

echo ${mmm} ${zzz} ${ttt} > 111

cat single.zero >> 111
mv 111 single.in

./sse.exe > sse.out 

mv single.in single.in_${mmm}_${zzz}_${ttt}
mv sse.out sse.out_${mmm}_${zzz}_${ttt}
mv single.dat single.dat_${mmm}_${zzz}_${ttt}

done
done
done
