#!/bin/bash
d=$1

echo "Starting field:" $d
echo ${d} | ./swan_hom_test > hom.out.${d}
echo "Finished field:" $d
grep homology hom.out.${d}
echo "============================================"

#parallel -j 4 --progress ./run_hom.sh {1} ::: 1103 1108 1111 1112 1115 1119 1123 1124 1128 1131 1135 1139 1140 1144 1147 1151 1155 1159 1160 1163 1167 1171 1172 1187 1191 1192 1195 1199

#parallel -j 5 --progress ./run_hom.sh {1} ::: 1 2 3 7 11
