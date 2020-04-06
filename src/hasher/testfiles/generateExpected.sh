#! /bin/bash

#EXE="~/src/git/github/hasher2.29"
EXE="../../../cmake-build-debug/src/hasher/hasher2.30"
$EXE -w 15 -k 4 -d 5 -S 7 -s 50 -c 100 -v 1  DmelChr4.fa DmelChr4_denovoLibTEs.fa

mv DmelChr4.fa.hasher.align expDmelChr4.fa.hasher.align
rm DmelChr4.fa.*
rm *.kidx