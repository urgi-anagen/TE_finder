# The TE_finder suite Documentation
**Hadi Quesneville**


*Last modification: 12/03/19*

The TE_finder programs suite are C++ programs developed for transposable element search and their annotation in large 
eukaryotic genome sequence. But we think that the tools are generic 
enough to be used in other contexts. 

The suite is composed by four main programs: BLASTER, MATCHER, GROUPER and DUSTER. The results of 
BLASTER can be treated further by the MATCHER and GROUPER program, but these two last programs 
can be also used in conjunction with other programs such as RepeatMasker, CENSOR, BLAST, BLAT, ... 
if their outputs are correctly formated as a BLASTER output file (i.e. *.align).

Besides, several utility tools are given as-is in the TOOLS section. 

## BLASTER

BLASTER compares two sets of sequences: a query databank against a subject databank. For each sequence in the query
databank, BLASTER launches one of the BLAST programs (blastn, tblastn, blastx, tblastx, blastp, megablast) (Altschul et
al., 1990, 1997) to search the subject databank. BLASTER can use either WUBLAST or NCBIBLAST programs. It cuts long
sequences before launching BLAST and reassembles the results afterwards. Hence, it can work on whole genomes, in
particular to compare a genome with itself to detect repeats. The results of BLASTER can then be treated by the MATCHER
and GROUPER programs described below.

### Commands 
``` 
usage: blaster [option] 
Options are: 
-q: 
query bank name, default=<not set> 
-s: 
subject bank name, default=<not set> 
-a: 
all_by_all, default=FALSE 
-W: 
use WU-BLAST, default=NCBI-BLAST 
-n: 
blast name, default=blastn 
-p: 
blast options, default= 
-l: 
fragment cutter length, default=50000 
-o: 
fragment cutter overlap, default=100 
-w: 
minimum cutter word, default=11 
-e: 
extention of cutted file, default=_cut 
-S: 
sensitivity, default=0 (3 levels for blastn, 3 being the most 
sensitive) 
-E: 
e-value filter, default=1e-10 
-I: 
identity filter, default=0 
-L: 
length filter, default=20 
-b: 
output filename prefix, default=<not set> 
-B: 
output filename prefix (no time stamp), default=<not set> 
-r: 
force re-run all, default=Off 
-P:
only prepare the data default=Off 
```
### Outputs: 
Several files are produced by BLASTER. 

* *_cut : preprocessed file for use with blast (cut in chunk and remove_self_hits N stretches) 
* *_cut.nhr, *_cut.nin, *_cut.nsq, ... : bank preprocessed by the NCBI blast formatdb 
program 
* *_cut.xnd, *_cut.xns, *_cut.xnt, ... : bank preprocessed by the WU blast xdformat 
program 
* *.Nstretch.map : coordinates of the removed N stretches 
* *.param : dump of the used parameters 
* *.raw : raw BLAST results 
* *.seq_treated: sequences treated by BLASTER. Used to restart in case of failure. 
* last_time_stamp.log: a file containing the last_time stamp number used. 
* *.align : results of BLASTER in a tabular format. 

Columns from '*.align' files are separated by tabulations 

* col1: query sequence name 
* col2: match start coordinate on the query sequence 
* col3: match end coordinate on the query sequence 
* col4: subject sequence name 
* col5: match start coordinate on the subject sequence 
* col6: match end coordinate on the subject sequence 
* col7: BLAST Evalue 
* col8: BLAST score 
* col9: identity percentage

## MATCHER
MATCHER map the matches (HSPs) of the subject sequences on the queries. Cross hits are filtered as overlapping hits.
Here, overlapping means also included. So when two HSPs overlap on the genomic sequence, the one with the best alignment
score is kept, the other is truncated such as nonoverlapping region remains on the HSP. As a result of this procedure an
HSPs is totally removed only if it is included in a longer one with a best score. By this way, nested elements are kept.
Long insertions (or deletions) in one of two homologous sequences result in two HSPs, instead of one with a long gap.
Thus the remaining HSP are linked by dynamic programming as in FASTA algorithm. A score is calculated by summing HSP
scores and substracting a gap penalty (x times the gap length) and a mismatch penalty (x times the mismatch length
region). The algorithm is inspired by the first part of the algorithm of sim2 described in Chao et al. (1995) CABIOS,
Vol. 11, pages 147153. But it is modified to make local alignment as follows.

1) An HSP is chained with a chain of others only if its score is less than the score of 
the resulting chain that links it. Thus, the chaining is stopped if the score of the 
resulting chain of HSP is less than the HSP not chained. The best scoring chain is 
kept. 
2) Then to identify other HSP chains, the chain previously found is removed, and we 
search again for the next best HSP chain. This is done iteratively until no chain is 
found. 
3) This algorithm is repeated independently for HSP matching in strand +/+, -/-, and -/+ 
A maximum of x% of overlap between the HSP is authorized. 
### Commands 
```
usage: matcher [option] 
Options are: 
-m: 
matches text file (in BLASTER *.align format) 
-q: 
query bank name, default=<not set> 
-s: 
subject bank name, default=<not set> 
-j: 
join matches, default=false 
-i: 
identity tolerance to join matches, default=2 
-g: 
gap penalty to join matches, default=0.05 
-d: 
distance penalty to join matches, default=0.2 
-c: 
authorized overlap to join matches, default=20 
-E: 
e-value filter, default=1e-10 
-I: 
identity filter, default=0 
-L: 
length filter, default=20 
-b: 
output filename prefix, default=<not set> 
-B: 
output filename prefix (no time stamp), default=<not set> 
-a, 
all conflicting subject default=FALSE
```
### Outputs 
*.clean_match.fa : sequences of the matching region. This is a concatenation of the joined 
fragments 

*.clean_match.map : coordinates of the whole matching region (boundaries of the joined 
fragments) on the query in the following tabulated format (the separator is a tabulation): 

* col1: the subject sequence name 
* col2: the query sequence name 
* col3: start coordinate 
* col4: end coordinate

*.clean_match.param : dump of the used parameters for MATCHER 

*.clean_match.path : results of MATCHER in a tabular format. Columns are separated by 
tabulations. 
* col1: path number. Joined fragments have the same path number 
* col2: query sequence name 
* col3: match start coordinate on the query sequence 
* col4: match end coordinate on the query sequence 
* col5: subject sequence name 
* col6: match start coordinate on the subject sequence 
* col7: match end coordinate on the subject sequence 
* col8: BLAST Evalue 
* col9: BLAST score 
* col10: identity percentage 

*.clean_match.tab : summary results of MATCHER in a tabular format. One line per joined 
fragment. Columns are separated by tabulations. 
* col1: query sequence name 
* col2: whole match start coordinate on the query sequence 
* col3: whole match end coordinate on the query sequence 
* col4: length on the query sequence 
* col5: length in percentage of the query sequence 
* col6: length on the query relative to the subject length in percentage 
* col7: subject sequence name 
* col8: whole match start coordinate on the subject sequence 
* col9: whole match end coordinate on the subject sequence 
* col10: BLAST Evalue 
* col11: BLAST score 
* col12: identity percentage 
* col13: path number 

## GROUPER
GROUPER uses HSPs to gather similar sequences into groups by simple link clustering. An alignment belongs to a group if
one of the two aligned sequences already belongs to this group over a given coverage percentage (*-C* parameter). Groups
that share sequence locations are regrouped into what we called a cluster. As a result of these procedures, each group
contains sequences that are homogeneous in length. A given region may belong to several groups, but all of these groups
belong to the same cluster.
### Commands 
```
usage: grouper [option] 
Options are: 
-m: 
match text file(in BLASTER *.align format) 
-q:
query bank name, default=<not set> 
-s: 
subject bank name, default=<not set> 
-j: 
join matches, default=false 
-i: 
identity tolerance to join matches, default=2 
-g: 
gap penalty to join matches, default=0.05 
-d: 
distance penalty to join matches, default=0.2 
-c: 
authorized overlap to join matches, default=20 
-E: 
e-value filter, default=1e-10 
-I: 
identity filter, default=0 
-L: 
length filter, default=20 
-b: 
output filename prefix, default=<not set> 
-B: 
output filename prefix (no time stamp), default=<not set> 
-G: 
min coverage length for connecting groups in a cluster,default=100 
-S: 
connecting groups in a cluster if all members of one group overlap 
members of another group, default=FALSE 
-C: 
coverage for group construction default=0.95 
```
### Outputs
Output files have the following prefix *.align.group.c0.[0-9]+. The number after 
*.align.group.c indicates the coverage threshold used. 

*.align.group.c0.95.fa: repeated sequences found in fasta format. Each repeated sequence is 
named according to the following nomenclature: Mb[SQ][09]+
Gr[09]+
Cl[09]+ 
* "Mb" followed by "S"or "Q" according to the fact that this sequence was found in a 
query (Q) or a subject (S) sequence and a number identifying the sequence 
* "Gr" followed by the group number 
* "Cl" followed by the cluster number 
After this name, follows the "fasta header" of the source sequence and the coordinates on it. 

*.align.group.c0.95.map: coordinates of the sequences found. Columns are separated by 
tabulations. 
* col1: name of the repeated sequence fragment (as for fasta) 
* col2: query sequence name 
* col3: start coordinate on the query sequence 
* col4: end coordinate on the sequence 

*.align.group.c0.95.set: coordinates in set format, indicating connected fragments. Columns 
are separated by tabulations.
* col1: path number. Joined fragments have the same path number 
* col2: name of the repeated sequence fragment (as for fasta) 
* col3: query sequence name 
* col4: start coordinate on the query sequence 
* col5: end coordinate on the sequence 

*.align.group.c0.95.param : parameters used for grouper 

*.align.group.c0.95.txt: detailed description of the groups and clusters 

*.align.group.c0.95.cluster.dot: Graph representation of groups links within clusters in 
dot format (use the dot program of graphviz package to visualize http://www.graphviz.org/) 

*.align.group.c0.95.overlap.dot: Graph representation of sequences links in dot format 
(useful only for small example) 

## DUSTER

DUSTER is able to compare a genome sequence, here considered as a query sequence, to a large amount of TE sequences,
i.e. a sequence library. Its algorithm used k-mers to search for similar sequences without performing nucleotid
aligments. Sensitivity is obtained allowing one mismatches in k-mers every *n* consecutive nucleotids. In summary the
algorithm compares k-mers between the genome and each sequence from the library and reports matches when at least two
k-mers are found on the same alignment diagonal (i.e. the difference between the coordinates on the query and the
sequence library are identical) with a maximal distance of *d* k-mers. The region bounded by the two-extreme k-mer
position are reported as matching. Two matching regions on the genome separated by less than *x* k-mers are merged. At
the end of this first pass, the region identified on the genome can be used as a new sequence library for a new search
(the *-n* parameter). This procedure is repeated until genome coverage increased by less than 1% if *-n* is set to 0.
### Commands 
```
usage: duster [<options>] <fasta query sequence> [<fasta sequence model>].
 options:
   -h, --help:
         this help
   -w, --kmer:
         kmer length (<16), default: 15
   -S, --step_q:
         step on query sequence, default: 7
   -k, --mask_hole_period:
         period of k-mer mask, default: 4
   -d, --kmer_dist:
         max number of kmer between two matching kmer to connect, default: 5
   -f, --frag_connect_dist:
         max distance between two fragments to connect, default: 100
   -s, --min_size:
         min size range to report, default: 20
   -C, --filter_cutoff:
         filter kmer with counts over a percentile (Value [0-1]), default: 1
   -D, --diversity_cutoff:
         filter kmer with diversity measure of kmer size used for background probability (Value [0-1]), default: 0
   -m, --min_count:
         filter kmer with counts less than this value, default: 0
   -b, --background_kmer_size:
         kmer size to compute kmer background probability, default: 1
   -o, --file_out:
         filename for output,
   -c, --chunk_size:
         sequence chunk size in kb, default: None
   -n, --nb_iter:
         number of iteration. A value of 0 make iteration stop if coverage variation is less than 1%. default:0
   -a, --analysis:
         compute kmer statistics only
   -v, --verbosity:
         verbosity level, default:0
   -x, --filter_ssr_no:
         suppress SSR filter: No
```
### Output
Outputs have the genome filename as prefix, followed by a dot and a number indicating the iteration of the result, or
'final' for the last one. The following suffixes correspond to:

*.duster.bed : The coordinates of the regions on the genome in bed format.

*.duster.bed.fa : The sequences of the found regions.

A file with the TE sequence filename as prefix followed by '.kidx' is an index file of the kmer positions on the TE sequences. 
It is used when the analysis is re-run to avoid recomputing them to save time.


## TOOLS
Various basic tools.
 
## References
Main publications: 

BLASTER, MATCHER, GROUPER

* Quesneville H, Nouaud D, Anxolabehere D. Detection of new transposable element families in 
Drosophila melanogaster and Anopheles gambiae genomes. J Mol Evol. 2003;57 Suppl 1:S509. 
PMID: 15008403
* Quesneville H, Bergman CM, Andrieu O, Autard D, Nouaud D, Ashburner M, Anxolabehere D. 
Combined evidence annotation of transposable elements in genome sequences. PLoS Comput 
Biol. 2005 Jul;1(2):16675. 
Epub 2005 Jul 29. PMID: 16110336

DUSTER:

* Baud, Agnès; Wan, Mariène; Nouaud, Danielle; Francillonne, Nicolas; Anxolabéhère, Dominique; Quesneville, Hadi. Traces of transposable elements in genome dark matter co-opted by flowering gene regulation networks. Peer Community Journal, Volume 2 (2022), article no. e14. doi : 10.24072/pcjournal.68. https://peercommunityjournal.org/articles/10.24072/pcjournal.68/


