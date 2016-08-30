/* Modified by Hadi Quesneville 

from

   A MULTIPLE ALIGNMENT PROGRAM (MAP):

Copyright (c) 1992 Xiaoqiu Huang
All Rights Reserved.

Permission to use, copy, modify, and distribute this software and its
documentation for educational, research and non-profit purposes, without
fee, and without a written agreement is hereby granted, provided that the
above copyright notice, this paragraph and the following three paragraphs
appear in all copies.

Permission to incorporate this software into commercial products may be
obtained from the Intellectual Property Office, 1400 Townsend Drive, 301
Administration Building, Houghton, MI  49931, phone (906) 487-3429.

IN NO EVENT SHALL THE AUTHOR OR MICHIGAN TECHNOLOGICAL UNIVERSITY BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
DOCUMENTATION, EVEN IF MICHIGAN TECHNOLOGICAL UNIVERSITY HAS BEEN ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

MICHIGAN TECHNOLOGICAL UNIVERSITY SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN
"AS IS" BASIS, AND MICHIGAN TECHNOLOGICAL UNIVERSITY HAS NO OBLIGATIONS TO
PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

    Proper attribution of the author as the source of the software would
    be appreciated: "On global sequence alignment" (to appear in CABIOS).
	      Xiaoqiu Huang
	      Department of Computer Science
	      Michigan Technological University
	      Houghton, MI 49931
              E-mail: huang@cs.mtu.edu

    The MAP program computes a multiple global alignment of sequences using
    iterative pairwise method. The underlying algorithm for aligning
    two sequences computes a best overlapping alignment bewteen
    two sequences without penalizing terminal gaps. In addition,
    long internal gaps in short sequences are not heavily penalized.
    So MAP is good at producing an alignment where there are long
    terminal or internal gaps in some sequences. The MAP program is
    designed in a space-efficient manner, so long sequences can be aligned. 

    Users supply scoring parameters. In the simplest form, users
    provide 3 integers: ms, q and r, where ms is the score of a mismatch
    and the score of an i-symbol indel is -(q + r * i). Each match
    automatically receives score 10. In addition, an integer gs is
    provided so that any gap of length > gs in a short sequence is 
    given a penalty of -(q + r * gs), the linear penalty for a gap of
    length gs. In other words, long gaps in the short sequence are
    given a constant penalty. This simple scoring scheme may be used
    for DNA sequences.  NOTE: all scores are integers.

    In general, users can define an alphabet of characters to appear
    in the sequences and a matrix that gives the substitution score
    for each pair of symbols in the alphabet. The 127 ASCII characters
    are eligible. The alphabet and matrix are given in a file, where
    the first line lists the characters in the alphabet and the lower
    triangle of the matrix comes next. An example file looks as follows:

    ARNDC	       
     13
    -15  19
    -10 -22  11
    -20 -10 -20  18
    -10 -20 -10 -20  12

    Here the -22 at position (3,2) is the score of replacing N by R.
    This general scoring scheme is useful for protein sequences where the
    set of protein characters and Dayhoff matrix are specified in the file.
    Note that the characters in the alphabet must be exactly the same
    (including lower or upper cases) as ones appearing in sequences.

    The MAP program is written in C and runs under Unix systems on
    Sun workstations and under DOS systems on PCs.
    We think that the program is portable to many machines.

    Sequences to be aligned are stored in one file.
    A sample file of sequences looks like:
>Human-beta
VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV
KAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK
EFTPPVQAAYQKVVAGVANALAHKYH
>Horse-beta
VQLSGEEKAAVLALWDKVNEEEVGGEALGRLLVVYPWTQRFFDSFGDLSNPGAVMGNPKV
KAHGKKVLHSFGEGVHHLDNLKGTFAALSELHCDKLHVDPENFRLLGNVLVVVLARHFGK
DFTPELQASYQKVVAGVANALAHKYH
>Human-alpha
VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGK
KVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPA
VHASLDKFLASVSTVLTSKYR
>Horse-alpha
VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHFDLSHGSAQVKAHGK
KVGDALTLAVGHLDDLPGALSNLSDLHAHKLRVDPVNFKLLSHCLLSTLAVHLPNDFTPA
VHASLDKFLSSVSTVLTSKYR
>Sea-lamprey
PIVDTGSVAPLSAAEKTKIRSAWAPVYSDYETSGVDILVKFFTSTPAAEEFFPKFKGLTT
ADELKKSADVRWHAERIIDAVDDAVASMDDTEKMSSMKDLSGKHAKSFEVDPEYFKVLAA
VIADTVAAGDAGFEKLLRMICIL
LRSAY
>Sperm-whale
VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASED
LKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHP
GDFGADAQGAMNKALELFRKDIAAKYKELG
YQG
>Yellow-lupin
GALTESQAALVKSSWEEFNANIPKHTHRFFILVLEIAPAAKDLFSSFLKGGTSEVPQNNPE
LQAHAGKVFKLVYEAAIQLEVTGVVASDATLKNLGSVHVSKGVVADAHFPVVKEAILKTIK
EVVGAKWSEELNSAWTIAYDELAIVIKKEMDDAA
    The string after ">" is the name of the following sequence.

    To find the best alignment of sequences in file A,
    use a command of form

	   map  A  gs  ms  q  r > result

    where map is the name of the object code, gs is the minimum length
    of any gap in a short sequence charged with a constant gap penalty,
    ms is a negative integer specifying mismatch weight, q and r are
    non-negative integers specifying gap-open and gap-extend penalties,
    respectively. Output alignments are saved in the file "result".

    For using a scoring matrix defined in file S, use a command of form

	   map  A  gs  S  q  r > result

    Note that ms is replaced by the file S.

    Acknowledgments
    The function diff() evolved from that by Gene Myers.
    The author thanks Chunwei Wang for pointing out the problem
    with existing multiple alignment software.
    The author also thanks Dave Gordon and John Hunt for suggesting
    that the alignment be produced in flat and interleaved formats
    so that it can be read by some phylogenetic analysis programs.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define	NAMELEN  40	/* can't exceed 40 */

static int match, mismh;	/* max and min substitution weights */
static int SCORE;	/* score of a best alignment */
static int STARI;	/* start row number */
static int STARJ;	/* start column number */
static int ENDI;	/* end row number */
static int ENDJ;	/* end column number */
static int v[128][128];	/* substitution scores */
static int  q, r;       /* gap penalties */
static int  qr;         /* qr = q + r */
static int  gaplen;     /* minimum length for constant-cost insertion */
static int  pay;        /* constant-cost for long insertion */

typedef struct OVERLAP  /* of 5' and 3' segments */
{  int    id1; 	        /* id of sequence 1   */
   int    id2;          /* id of sequence 2   */
   int    score;        /* score of overlap alignment */
   struct  OVERLAP  *next; } over, *overptr;
struct SEG
{  char    *name;       /* name of sequence */
   int     len;         /* length of ssequence name */
   char    *seq;        /* sequence */
   int     length;      /* length of sequence */
   overptr list;        /* list of overlapping edges */
        } *segment;     /* array of sequence records */
int   seg_num;          /* The number of segments   */
int   maxlen;		/* maximum of sequence lengths */
overptr   *edge;        /* set of overlapping edges */
int   edge_num;         /* The number of overlaps */
struct ALG 
{  char    *row[2];	/* one row of multiple alignment */
   int     len[2];      /* maximum length of row */
   int     flag;        /* row[flag] is being used; flag = 0 or 1 */
   int     next;        /* id of next sequence in alignment */
   int     head;        /* id of the first sequence in alignment */
   int     num;    	/* number of sequences in this group if it is head */
   int     length;	/* length of alignment */
   int     pos;		/* current position in layout */
   char	   **rect;	/* pointer to 2-D alignment matrix */
           } *group;    /* sequence groups */

int main(argc, argv) int argc; char *argv[];
{ int   M;        			/* Sequence length */
  int  ms;				/* User-supplied weights      */
  FILE *Ap, *Sp, *ckopen();
  char *ckalloc();			/* space-allocating function  */
  char  alph[129], *s;			/* alphabet */
  int  size;				/* size of alphabet */
  int   total;			        /* Total of sequence lengths */
  int   number;                         /* The number of sequences   */
  char  *sequence;			/* Storing all sequences     */
  int   symbol, prev;			/* The next character         */
  int  i, j;				/* index variables	      */
  short  heading;			/* 1: heading; 0: sequence    */

	if ( argc != 6 )
	   fatalf("The proper form of command is: \n%s file gap-size(>0) mismatch(<0 /match=10) gap-open(>=0) gap-extend(>=0)", argv[0]);
	/* determine number of sequences and total lengths */
	j = maxlen = 0;
	Ap = ckopen(argv[1], "r");
	prev = '\n';
	for (total = 3, number = 0; ( symbol = getc(Ap)) != EOF ; total++ )
	  { if ( symbol == '>' && prev == '\n' )
	      number++;
	    prev = symbol;
	  }
	if ( number == 0 )
	   fatal("There are no sequences or sequences are in wrong format");
	total  += number * 20;
	/* allocate space for sequences */
	sequence = ( char * ) ckalloc( total * sizeof(char));
	segment = ( struct SEG * ) ckalloc( number * sizeof(struct SEG));
	/* read the sequences into sequence */
	M = 0;
	Ap = ckopen(argv[1], "r");
	number = -1;
        heading = 0;
	prev = '\n';
	for ( i = 0; ( symbol = getc(Ap)) != EOF ; )
	  { if ( symbol != '\n' )
	       sequence[++i] = symbol;
	    if ( symbol == '>' && prev == '\n' )
	      { heading = 1;
		if ( number >= 0 )
		  { segment[number].length = i - j - 1;
		    if ( maxlen < i - j - 1 )
		       maxlen = i - j - 1;
		    if ( i - j - 1 > M ) M = i - j -1;
		  }
		number++;
		j = i;
		segment[number].name = &(sequence[i]);
		segment[number].list = NULL;
	      }
	    if ( heading && symbol == '\n' )
	      { heading = 0;
		segment[number].len = i - j;
		segment[number].seq = &(sequence[i]);
		j = i;
	      }
	    prev = symbol;
	  }
	segment[number].length = i - j;
	if ( maxlen < i - j )
	    maxlen = i - j;
        if ( i - j > M ) M = i - j;
	seg_num = ++number;
 (void)	fclose(Ap);
        edge_num = 0;

	(void) sscanf(argv[2],"%d", &gaplen);
	if ( gaplen < 1 )
	  fatal("The minimum length for constant-cost insertion is a positive integer");

	(void) sscanf(argv[argc-2],"%d", &q);
	if ( q < 0 )
	   fatal("The gap-open penalty is a nonnegative integer");

	(void) sscanf(argv[argc-1],"%d", &r);
	if ( r < 0 )
	   fatal("The gap-extend penalty is a nonnegative integer");

	pay = q + r * gaplen;
	qr = q + r;
	/* check if the argument represents a negative integer */
	s = argv[argc-3];
	if ( *s == '-' ) s++;
	for ( ; *s >= '0' && *s <= '9' ; s++ );
	if ( *s == '\0' )
	  { (void) sscanf(argv[argc-3],"%d", &ms);
	    if ( ms >= 0 )
	       fatal("The mismatch weight is a negative integer");
	    match = 10;
	    mismh = ms;
	    /* set match and mismatch weights */
	    for ( i = 0; i < 128 ; i++ )
	      for ( j = 0; j < 128 ; j++ )
	         if (i == j )
	            v[i][j] = 10;
	         else
	            v[i][j] = mismh;
	  }
	else
	  { /* read a file containing alphabet and substitution weights */
	    Sp = ckopen(argv[argc-3], "r");
	    (void) fscanf(Sp, "%s", alph);
	    size = strlen(alph);
	    match = mismh = 0;
	    for ( i = 0; i < 128 ; i++ )
	      for ( j = 0; j < 128 ; j++ )
                  v[i][j] = 0;
	    for ( i = 0; i < size ; i++ )
	      for ( j = 0; j <= i ; j++ )
		{ (void) fscanf(Sp, "%d", &ms);
		  v[alph[i]][alph[j]] = v[alph[j]][alph[i]] = ms;
		  if ( ms > match ) match = ms;
		  if ( ms < mismh ) mismh = ms;
		}
	  }
	for ( i = 0; i < 128 ; i++ )
	    v['-'][i] = v[i]['-'] = - r;
	v['-']['-'] = 0;
        Pairwise(total);
        Multiple();
        /* Show(); */
	FlatFormat();
	/* InterFormat(); */
	return 0;
}

static int *CC, *DD;			/* saving matrix scores */
static int *RR, *SS;		 	/* saving start-points */
static int  *S;				/* saving operations for diff */

/* The following definitions are for function diff() */

int  diff();
static int  zero = 0;				/* int type zero        */

#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel score */

#define gap2(k)  ((k) <= 0 ? 0 : ((k) <= gaplen ? q+r*(k) : pay))
/* k-symbol insertion score */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */

						/* Append "Delete k" op */
#define DEL(k)				\
{ if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
}
						/* Append "Insert k" op */
#define INS(k)				\
{ if (last > 0)				\
    last = sapp[-1] += (k);             \
  else					\
    last = *sapp++ = (k);		\
}

						/* Append "Replace" op */
#define REP 				\
{ last = *sapp++ = 0; 			\
}

/* Perform pair-wise comparisons of sequences. */
Pairwise(total)
int  total;		/* total sequence length */
{ int   i, j;		/* row and column indices  */
  char *A, *B;          /* pointers to sequences */
  int  M, N;		/* sequence lengths */
  overptr node1;	/* pointer to overlap */
  char *ckalloc();		/* space-allocating function */

	/* allocate space for all vectors */
	j = (total + 1) * sizeof(int);
	CC = ( int * ) ckalloc(j);
	DD = ( int * ) ckalloc(j);
	RR = ( int * ) ckalloc(j);
	SS = ( int * ) ckalloc(j);
	S = ( int * ) ckalloc(2 * j);
	for ( i = 0; i < seg_num - 1 ; i++ )
	 { A = segment[i].seq;
	   M = segment[i].length;
	   for ( j = i+1; j < seg_num ; j++ )
	    { B = segment[j].seq;
	      N = segment[j].length;
	      node1 = ( overptr ) ckalloc( (int ) sizeof(over));
	      SCORE = - ( 2 * q + (M + N) * r + 1000);
	      if ( M <= N )
	       { big_pass(A,B,M,N);
	         node1->id1 = i;
	         node1->id2 = j;
	         node1->score = SCORE;
	         node1->next = segment[i].list;
                 segment[i].list = node1;
	       }
	      else
	       { big_pass(B,A,N,M);
	         node1->id1 = j;
	         node1->id2 = i;
	         node1->score = SCORE;
	         node1->next = segment[j].list;
                 segment[j].list = node1;
	       }
	      edge_num++;
	    }
         }
}

/* find best overlap score between two sequences */
big_pass(A,B,M,N) char A[],B[]; int M,N;
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci;			/* end-point associated with c */ 
  register  int  di;			/* end-point associated with d */
  register  int  fi;			/* end-point associated with f */
  register  int  pi;			/* end-point associated with p */
  register  int  g;     /* best score ending with constant-cost insertion */
  register  int  gi;                    /* end-point associated with g */
  int  t;                               /* temporary variable */
  int  *va;				/* pointer to v(A[i], B[j]) */

	/* Compute the matrix.
	   CC : the scores of the current row
	   RR : the starting point that leads to score CC
	   DD : the scores of the current row, ending with deletion
	   SS : the starting point that leads to score DD        */
 	/* Initialize the 0 th row */
	for ( j = 1; j <= N ; j++ )
	  {  CC[j] = 0;
	     DD[j] = - (q);
	     RR[j] = SS[j] = -j;
	  }
	CC[0] = 0;
	for ( i = 1; i <= M; i++) 
	  {  c = 0;				/* Initialize column 0 */
	     f = - (q);
	     g = - pay;
	     ci = fi = gi = i;
	     p = 0;
	     pi = i - 1;
	     va = v[A[i]];
	     for ( j = 1 ; j <= N ; j++ )  
	       { if ( ( f = f - r ) < ( c = c - qr ) )
		    { f = c; fi = ci; }
		  di = SS[j];
		  if ( ( d = DD[j] - r ) < ( c = CC[j] - qr ) )
		    { d = c; di = RR[j]; } 
		  c = p+va[B[j]];		/* diagonal */
		  ci = pi;
		  if ( c < d )
		    { c = d; ci = di; }
		  if ( c < f )
		    { c = f; ci = fi; }
		  if ( j - gaplen > 0 )
		   { if ( g < ( t = CC[j-gaplen-1] - pay ) )
		       { g = t; gi = RR[j-gaplen-1]; }
		     if ( c < g )
	       	       { c = g; ci = gi; }
		   }
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  RR[j] = ci;
		  DD[j] = d;
		  SS[j] = di;
	          if ( ( j == N || i == M ) &&  c > SCORE )
		   { SCORE = c;
		     ENDI = i;
		     ENDJ = j;
		     STARI = ci;
		   }
	        }
	  }
	if ( STARI < 0 )
	 { STARJ = - STARI;
	   STARI = 0;
	 }
	else
	   STARJ = 0;
}

/* Construct mutiple alignments */
Multiple()
{ char *ckalloc();	/* space-allocating function */
  int   i, j, k, t;	/* index variables */
  overptr  node1;	/* temporary pointer */
  int   sorted;		/* boolean variable */
  char  *a, *b;		/* temporary pointers */
  int   head1, head2;	/* ids of first sequences in alignments */
  struct ALG *pa, *pb;	/* pointers to group elements */

	group = ( struct ALG * ) ckalloc( seg_num * sizeof(struct ALG));
	for ( i = 0; i < seg_num; i++ )
	 { group[i].row[0] = ( char * ) ckalloc( 4 * maxlen * sizeof(char));
	   group[i].row[1] = ( char * ) ckalloc( 4 * maxlen * sizeof(char));
	   group[i].len[0] = group[i].len[1] = 4 * maxlen; 
	   group[i].flag = 1;
	   group[i].next = -1;
	   group[i].head = i;
	   group[i].num = 1;
	   group[i].length = k = segment[i].length;
	   group[i].pos = 0;
	   group[i].rect = ( char ** ) ckalloc( sizeof(char *));
	   group[i].rect[0] = a = group[i].row[1];
	   a[0] = ' ';
	   for ( b = segment[i].seq, j = 1; j <= k; j++ )
	      a[j] = b[j];
	   a[j] = ' ';
	 }
	edge = ( overptr * ) ckalloc( edge_num * sizeof(overptr) );
	for ( j = 0, i = 0; i < seg_num; i++ )
	  for ( node1 = segment[i].list; node1 != NULL; node1 = node1->next )
	      edge[j++] = node1;
	edge_num = j;
	for ( i = edge_num - 1; i > 0; i-- )
	 { sorted = 1;
	   for ( j = 0; j < i; j++ )
	     if ( edge[j]->score < edge[j+1]->score )
	      { node1 = edge[j];
	        edge[j] = edge[j+1];
	        edge[j+1] = node1;
		sorted = 0;
	      }
	   if ( sorted )
	     break;
	 }
	for ( k = 0; k < edge_num; k++ )
	 { head1 = group[edge[k]->id1].head;
	   head2 = group[edge[k]->id2].head;
	   if ( head1 != head2 )
	    {  if ( group[head1].length > group[head2].length )
		{ t = head1;
		  head1 = head2;
		  head2 = t;
		}
	       pa = &group[head1];
	       pb = &group[head2];
	       sapp = S;
	       last = 0;
	       diff(pa->rect,pb->rect,pa->length,pb->length,q,q,
	                             pa->num,pb->num,0,0,0,0,0,0);
	       Merge(head1, head2, S);
	    }
	 }
}

/* Merge two sequence alignment according to script S */
Merge(head1, head2, S)
int	head1, head2;	/* ids of first sequences in two alignments */
int	S[];		/* script */
{  char	 **rect1, **rect2;	/* pointers to two input alignments */
   char  **rect;	/* pointer to resulting alignment */
   int   size1, size2;	/* number of sequences in two alignments */
   int   size;		/* number of sequences in resulting alignment */
   int   limit;		/* maximum length of resulting alignment */
   int   tail;		/* id of last sequence in alignment 1 */
   int   i, j, h, k;	/* index variables */
   int	 M, N;		/* lengths of alignments */
   int   op;		/* current script operation */
   char  *ckalloc();	/* space-allocating function */
   int   flag;		/* index of row */ 
	
	size1 = group[head1].num;
	size2 = group[head2].num;
	size = size1 + size2;
	rect = ( char ** ) ckalloc( size * sizeof(char *));
	limit = 2 + group[head1].length + group[head2].length;
	for ( h = 0, k = head1; k != -1; tail = k, k = group[k].next )
	 { group[k].flag = flag = 1 - group[k].flag;
	   if ( group[k].len[flag] < limit )
	    { group[k].row[flag] = (char *) ckalloc(2 * limit * sizeof(char));
	      group[k].len[flag] = 2 * limit;
	    }
	   rect[h++] = group[k].row[flag];
	 }
	group[tail].next = head2;
	for ( k = head2; k != -1; k = group[k].next )
	 { group[k].flag = flag = 1 - group[k].flag;
	   group[k].head = head1;
	   if ( group[k].len[flag] < limit )
	    { group[k].row[flag] = (char *) ckalloc(2 * limit * sizeof(char));
	      group[k].len[flag] = 2 * limit;
	    }
	   rect[h++] = group[k].row[flag];
	 }
	for ( h = 0; h < size; h++ )
	   rect[h][0] = ' ';
	rect1 = group[head1].rect;
	rect2 = group[head2].rect;
	k = 1;
	M = group[head1].length;
	N = group[head2].length;
	op = 0;
	i = 0;
	j = 0;
        while (i < M || j < N)
         { if (op == 0 && *S == 0)
            { op = *S++;
	      for ( ++i, h = 0; h < size1; h++ )
	         rect[h][k] = rect1[h][i];
	      for ( ++j, h = size1; h < size; h++ )
	         rect[h][k] = rect2[h-size1][j];
            }
           else
            { if (op == 0)
                   op = *S++;
              if (op > 0)
                   { for ( h = 0; h < size1; h++ )
		       if (rect1[h][i] == ' ' || rect1[h][i+1] == ' ')
	                  rect[h][k] = ' ';
		       else
	                  rect[h][k] = '-';
	             for ( ++j, h = size1; h < size; h++ )
	                rect[h][k] = rect2[h-size1][j];
                     op--;
                   }
              else
                   { for ( ++i, h = 0; h < size1; h++ )
	                rect[h][k] = rect1[h][i];
	             for ( h = size1; h < size; h++ )
		       if (rect2[h-size1][j] == ' ' || rect2[h-size1][j+1] == ' ')
	                  rect[h][k] = ' ';
		       else
	                  rect[h][k] = '-';
                     op++;
                   }
            }
	   k++;
         }
	for ( h = 0; h < size; h++ )
	   rect[h][k] = ' ';
	group[head1].num = size;
	group[head1].length = k - 1;
	group[head1].rect = rect;
}

/* diff(A,B,M,N,tb,te,U,W,mm,nn,sc,sr,ec,er) returns the score of an optimum conversion
   between A[0..U-1][mm+1..mm+M] and B[0..W-1][nn+1..nn+N] that begins(ends) with
   a delete if tb(te) is zero and appends such a conversion to the current script.
   If sc = 0, then the beginning deletion is not penalized; if sr = 0, the beginning
   insertion is not penalized; if ec = 0, the ending deletion is not charged;
   if er = 0; then the ending insertion is not charged.  Any insertion of length
   at least gaplen is given a constant cost */

int diff(A,B,M,N,tb,te,U,W,mm,nn,sc,sr,ec,er) char *A[], *B[]; int M, N;
int tb, te, U, W, mm, nn, sc, sr, ec, er;
{ int   midi, midj, type;	/* Midpoint, type, and cost */
  int midc;
  int  ss,cc;

{ register int   i, j;
  register int c, e, d, s;
           int t, *va;
	   int  g, temp;
  	   char  *ckalloc();
  int  h,k;			/* index variables */
  int  tt;			/* temporary variable */
  int  x,y;			/* number of non-blank symbols in columns */

/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) DEL(M)
      if ( !sc || !ec )
	 return 0;
      else
         return - gap(M);
    }
  if (M <= 1)
    { if (M <= 0)
        { INS(N);
	  if ( !sr || !er )
            return 0;
          else
            return - gap2(N);
        }
      midc = - (sc * (tb + r) + er * gap2(N) );
      midj = -1;
      if ( midc < ( c =  - (ec * (te + r) + sr * gap2(N) ) ) )
        { midc = c;
          midj = 0;
        }
      for ( x = h = 0; h < U; h++ )
        if ( A[h][mm+1] != ' ' )
           x += 1;
      for (j = 1; j <= N; j++)
	{ c = 0;
	  for ( h = 0; h < U; h++ )
	    if ( ( tt = A[h][mm+1] ) != ' ' )
	     { va = v[tt];
	       for ( k = 0; k < W; k++ )
		 if ( ( tt = B[k][nn+j] ) != ' ' )
		    c += va[tt];
             }
	  for ( y = k = 0; k < W; k++ )
	     if ( B[k][nn+j] != ' ' )
		y += 1;
	  c = c / (x * y) - ( sr * gap2(j-1) + er * gap2(N-j) );
          if (c > midc)
           { midc = c;
             midj = j;
           }
	}
      if (midj == -1)
        { DEL(1) INS(N) }
      else
      if (midj == 0)
        { INS(N) DEL(1) }
      else
        { if (midj > 1) INS(midj-1)
          REP
          if (midj < N) INS(N-midj)
        }
      return midc;
    }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;			/* Forward phase:                          */
  CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
  t = -q * sr;
  if ( N <= gaplen )
    for (j = 1; j <= N; j++)
      { CC[j] = t = (t-r) * sr;
        DD[j] = t-q;
      }
  else
   { for (j = 1; j <= gaplen; j++)
       { CC[j] = t = (t-r) * sr;
         DD[j] = t-q;
       }
     for (j = gaplen+1; j <= N; j++)
       { CC[j] = t = -pay * sr;
         DD[j] = t - q;
       }
   }
  if ( !ec ) DD[N] += q;
  t = -tb * sc;
  for (i = 1; i <= midi; i++)
    { s = CC[0];
      CC[0] = c = t = (t-r) * sc;
      e = t-q;
      g = t - pay;
      for ( x = h = 0; h < U; h++ )
        if ( A[h][mm+i] != ' ' )
           x += 1;
      for (j = 1; j <= N; j++)
        { if ((c = c - qr) > (e = e - r)) e = c;
	  if ( j == N && !ec )
	    { if ((c = CC[j] ) > (d = DD[j] )) d = c;}
	  else
            if ((c = CC[j] - qr) > (d = DD[j] - r)) d = c;
	  c = 0;
	  for ( h = 0; h < U; h++ )
	    if ( ( tt = A[h][mm+i] ) != ' ' )
	     { va = v[tt];
	       for ( k = 0; k < W; k++ )
		 if ( ( tt = B[k][nn+j] ) != ' ' )
		    c += va[tt];
             }
	  for ( y = k = 0; k < W; k++ )
	     if ( B[k][nn+j] != ' ' )
		y += 1;
	  c = c / (x * y) + s;
          if (c < d) c = d;
          if (c < e) c = e;
	  if ( j - gaplen > 0 )
	    { if ( g < ( temp = CC[j-gaplen-1] - pay ) )
	        g = temp;
	      if ( c < g ) c = g;
	    }
          s = CC[j];
          CC[j] = c;
          DD[j] = d;
        }
    }
  DD[0] = CC[0];

  RR[N] = 0;			/* Reverse phase:                          */
  t = -q * er;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
  if ( N <= gaplen )
    for (j = N-1; j >= 0; j--)
      { RR[j] = t = (t-r) * er;
        SS[j] = t-q;
      }
  else
    { temp = N - gaplen;
      for (j = N-1; j >= temp; j--)
        { RR[j] = t = (t-r) * er;
          SS[j] = t-q;
        }
      for (j = temp-1; j >= 0; j--)
        { RR[j] = t = -pay * er;
          SS[j] = t - q;
        }
    }
  if ( !sc ) SS[0] += q;
  t = -te * ec;
  for (i = M-1; i >= midi; i--)
    { s = RR[N];
      RR[N] = c = t = (t-r) * ec;
      g = t - pay;
      e = t-q;
      for ( x = h = 0; h < U; h++ )
        if ( A[h][mm+i+1] != ' ' )
           x += 1;
      for (j = N-1; j >= 0; j--)
        { if ((c = c - qr) > (e = e - r)) e = c;
          if ( !j && !sc )
            { if ((c = RR[j] ) > (d = SS[j] )) d = c;}
          else
            if ((c = RR[j] - qr) > (d = SS[j] - r)) d = c;
	  c = 0;
	  for ( h = 0; h < U; h++ )
	    if ( ( tt = A[h][mm+i+1] ) != ' ' )
	     { va = v[tt];
	       for ( k = 0; k < W; k++ )
		 if ( ( tt = B[k][nn+j+1] ) != ' ' )
		    c += va[tt];
             }
	  for ( y = k = 0; k < W; k++ )
	     if ( B[k][nn+j+1] != ' ' )
		y += 1;
	  c = c / (x * y) + s;
          if (c < d) c = d;
          if (c < e) c = e;
	  if ( j + gaplen < N )
	    { if ( g < ( temp = RR[j+gaplen+1] - pay ) )
	         g = temp;
 	      if ( c < g ) c = g;
	    }
          s = RR[j];
          RR[j] = c;
          SS[j] = d;
        }
    }
  SS[N] = RR[N];

  midc = CC[0]+RR[0];		/* Find optimal midpoint */
  midj = 0;
  type = 1;
  for (j = 0; j <= N; j++)
    if ((c = CC[j] + RR[j]) >= midc)
      if (c > midc || CC[j] != DD[j] && RR[j] == SS[j])
        { midc = c;
          midj = j;
        }
  for (j = N; j >= 0; j--)
    { if ( j == N )
	 d = q * ec;
      else
        if ( j == 0 )
          d = q * sc;
        else
          d = q;
      if ((c = DD[j] + SS[j] + d) > midc)
       { midc = c;
         midj = j;
         type = 2;
       }
    }
}

/* Conquer: recursively around midpoint */

  cc = midj == N ? ec : 1;
  ss = midj == 0 ? sc : 1;
  if (type == 1)
    { diff(A,B,midi,midj,tb,q,U,W,mm,nn,sc,sr,cc,1);
      diff(A,B,M-midi,N-midj,q,te,U,W,mm+midi,nn+midj,ss,1,ec,er);
    }
  else
    { diff(A,B,midi-1,midj,tb,zero,U,W,mm,nn,sc,sr,cc,1);
      DEL(2);
      diff(A,B,M-midi-1,N-midj,zero,te,U,W,mm+midi+1,nn+midj,ss,1,ec,er);
    }
  return midc;
}

/* Alignment display routine */

Show()
{ int   i, j, k, h, n;
  int  head;
  int  length;
  int  len;
  char name[NAMELEN+3];
  char *t;
  char **rect;
  char a;
  int  line;
  char stg[200];
  int  blank;

	head = group[0].head;
	rect = group[head].rect;
	length = group[head].length;
	line = 0;
	for ( i = 1; i <= length; )
	 { (void) printf("               ");
           for ( j = 0; j < 60; j += 10 )
	       (void) printf("    .    :");
	   line += 60 ;
	   (void) printf("%7d\n", line);
	   for ( j = head, h = 0; j != -1; j = group[j].next, h++ )
	    { len = segment[j].len;
              t = segment[j].name + 1;
              for ( k = 0; k < len && k < NAMELEN; k++ )
                 name[k] = *t++;
              name[k] = '\0';
	      t = stg;
	      blank = 1;
              (void) sprintf(t,"%-15s", name);
	      t += 15;
	      for ( k = i; k <= i+59 && k <= length; k++ )
	       { (void) sprintf(t++,"%c", (a = rect[h][k]) );
		 if ( a != ' ' && a != '-' )
		    group[j].pos++;
		 if ( a != ' ')
		   blank = 0;
	       }
	      for ( n = k; n <= i+59; n++ )
	        (void) sprintf(t++,"%c", ' ');
              //(void) sprintf(t,"%7d\n\0", group[j].pos);
	      (void) sprintf(t,"%7d\n", group[j].pos);
	      if ( ! blank )
                 (void) printf("%s", stg);
	    }
	   i = k; 
           (void) printf("\n");
	 }
}

/* Display output in a flat format */

FlatFormat()
{ int  i, j, k, h;
  int  head;
  int  length;
  int  len;
  char name[NAMELEN+3];
  char *t;
  char **rect;
  char a;
  char stg[200];

	head = group[0].head;
	rect = group[head].rect;
	length = group[head].length;
   	for ( j = head, h = 0; j != -1; j = group[j].next, h++ )
	  { len = segment[j].len;
            t = segment[j].name + 1;
            for ( k = 0; k < len && k < NAMELEN; k++ )
              name[k] = *t++;
            name[k] = '\0';
            (void) printf(">%s\n", name);
	    for ( i = 1; i <= length; )
	     { t = stg;
	       for ( k = i; k <= i+59 && k <= length; k++ )
	        { a = rect[h][k];
		  if ( a == ' ' )
		     a = '-';
		  *t++ = a;
	        }
               *t++ = '\0';
               (void) printf("%s\n", stg);
	       i = k; 
	     }
	  }
}

/* Display alignment in an interleaved format */

InterFormat()
{ int   i, j, k, h;
  int  head;
  int  length;
  int  num;
  int  len;
  char name[NAMELEN+3];
  char *t, *s;
  char **rect;
  char a;
  char stg[200];

	head = group[0].head;
	rect = group[head].rect;
	length = group[head].length;
	num = group[head].num;
        (void) printf("\nAlignment in an interleaved format\n\n");
        (void) printf(" %d %d\n", num, length);
	for ( i = 1; i <= length; )
	 { for ( j = head, h = 0; j != -1; j = group[j].next, h++ )
	    { t = stg;
	      if ( i == 1 )
	       { len = segment[j].len;
                 s = segment[j].name + 1;
                 for ( k = 0; k < len && k < NAMELEN && k < 10; k++ )
                    name[k] = *s++;
                 name[k] = '\0';
                 (void) sprintf(t,"%-10s   ", name);
		}
	      else
                 (void) sprintf(t,"             ");
	      t += 13;
	      for ( k = i; k <= i+59 && k <= length; k++ )
	       { a = rect[h][k];
		 if ( a == ' ')
		   a = '-';
		 *t++ = a;
	       }
              *t++ = '\0';
              (void) printf("%s\n", stg);
	    }
	   i = k; 
	 }
}

/* lib.c - library of C procedures. */

/* fatal - print message and die */
fatal(msg)
char *msg;
{
	fprintf(stderr, "%s\n", msg);
	exit(1);
}

/* fatalf - format message, print it, and die */
fatalf(msg, val)
char *msg, *val;
{
	fprintf(stderr, msg, val);
	putc('\n', stderr);
	exit(1);
}
	
/* ckopen - open file; check for success */
FILE *ckopen(name, mode)
char *name, *mode;
{
	FILE *fopen(), *fp;

	if ((fp = fopen(name, mode)) == NULL)
		fatalf("Cannot open %s.", name);
	return(fp);
}

/* ckalloc - allocate space; check for success */
char *ckalloc(amount)
int amount;
{
	char *p;

	if ((p = malloc( (unsigned) amount)) == NULL)
		fatal("Ran out of memory.");
	return(p);
}
