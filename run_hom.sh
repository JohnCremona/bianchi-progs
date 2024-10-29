#!/bin/bash
d=$1

export PARI_SIZE=20000000000

echo "Starting field" $d " at " `date`
T1=`date +%s`
echo ${d} | ./swan_hom_test > out/hom.out.new.${d}
T2=`date +%s`
echo "Finished field:" $d " at " `date`
let 'T = T2 - T1'
echo "Elapsed time for field $d : $T seconds"

#parallel -j 4 --progress ./run_hom.sh {1} ::: 1103 1108 1111 1112 1115 1119 1123 1124 1128 1131 1135 1139 1140 1144 1147 1151 1155 1159 1160 1163 1167 1171 1172 1187 1191 1192 1195 1199

#parallel -j 5 --progress ./run_hom.sh {1} ::: 1 2 3 7 11
