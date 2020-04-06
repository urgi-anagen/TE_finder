!# /bin/bash

~/src/git/github/duster2.29 -w 15 -k 4 -d 5 -f 100 -S 7 -n 1 DmelChr4.fa DmelChr4_denovoLibTEs.fa
mv DmelChr4.fa.1.duster.bed expDmelChr4.fa.1.duster.bed
rm DmelChr4.fa.*
rm *.kidx