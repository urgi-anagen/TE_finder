# TE Finder 2.30

 The TE Finder suite gathers C++ programs developed for transposable element search and their annotation in 
 large eukaryotic genome sequence. But we think that the tools are generic enough to be used in other contexts.
 
 The suite is mainly composed by four programs BLASTER, MATCHER, GROUPER, and DUSTER. 
 The results of BLASTER can be treated further by the MATCHER and GROUPER program, 
 but these two last programs can be also used in conjunction with other programs such as RepeatMasker, 
 CENSOR, BLAST, BLAT, ... if their outputs are correctly formated as a BLASTER output file (i.e. *.align).
 See the [documentation](Documentation.md)

The binaries are part of the [REPET package](http://urgi.versailles.inra.fr/Tools/REPET).

## Containers
Singularity images are available [here](https://cloud.sylabs.io/library/hquesneville)

You can run by simply typing:

-if the image is localy downloaded, type for example :

`$ singularity run te_finder2.30.sif blaster2.30 -h`

-if the image is remote, type for example :

`$ singularity run library://hquesneville/default/te_finder:2.30 blaster2.30 -h`
