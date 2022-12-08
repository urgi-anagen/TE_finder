# TE Finder 2.31

## Overview

 The TE Finder suite gathers C++ programs developed for transposable element search and their annotation in 
 large eukaryotic genome sequence. But we think that the tools are generic enough to be used in other contexts.
 
 The suite is mainly composed by three programs BLASTER, MATCHER, GROUPER, and DUSTER. 
 The results of BLASTER or HASHER can be treated further by the MATCHER and GROUPER program, 
 but these two last programs can be also used in conjunction with other programs such as RepeatMasker, 
 CENSOR, BLAST, BLAT, ... if their outputs are correctly formated as a BLASTER output file (i.e. *.align).
 See the [documentation](Documentation.md)

The binaries are part of the [REPET package](http://urgi.versailles.inrae.fr/Tools/REPET).

## Containers
Singularity images are available [here](https://cloud.sylabs.io/library/hquesneville)

You can run by simply typing:

-if the image is localy downloaded, type for example :

`$ singularity run te_finder2.31.sif blaster2.31 -h`

-if the image is remote, type for example :

`$ singularity run library://hquesneville/default/te_finder:2.31 blaster2.31 -h`
