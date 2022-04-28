/* A PROGRAM FOR LOCAL SIMILARITIES WITH AFFINE WEIGHTS:

    copyright (c) 1990,1991 Xiaoqiu Huang and Webb Miller
    The distribution of the program is granted provided no charge is made
    and the copyright notice is included.
    E-mail: huang@cs.mtu.edu

    Please cite the paper
    "A Time-Efficient, Linear-Space Local Similarity Algorithm"
	(Advances in Applied Mathematics, 12: 337-357, 1991)

	      Xiaoqiu Huang
	      Department of Computer Science
	      Michigan Technological University
	      Houghton, MI 49931

	      Webb Miller
	      Department of Computer Science
	      The Pennsylvania State University
	      University Park, PA 16802

     An early version of this program is presented in the paper
     "A Space-Efficient Algorithm for Local Similarities"
	(Computer Applications in the Biosciences, 6(4): 373-381)

	      Xiaoqiu Huang
	      Department of Computer Science
	      Michigan Technological University
	      Houghton, MI 49931

	      Ross C. Hardison and Webb Miller
	      Department of Molecular and Cell Biology
	      Department of Computer Science
	      The Pennsylvania State University
	      University Park, PA 16802

    The SIM program finds k best non-intersecting alignments between
    two sequences or within one sequence. Using dynamic programming
    techniques, SIM is guaranteed to find optimal alignments. The
    alignments are reported in order of similarity score, with the
    highest scoring alignment first. The k best alignments share no
    aligned pairs. SIM requires space proportional to the sum of the
    input sequence lengths and the output alignment lengths. Thus
    SIM can handle sequences of tens of thousands, even a hundred of
    thousands, of base pairs on a workstation. For example,
    on 73,360-bp and 44,594-bp sequences, SIM took 15 hours to find
    100 best local alignments on a Sun4 workstation.

    Users supply scoring parameters. In the simplest form, users just
    provide 3 integers: ms, q and r, where ms is the score of a mismatch
    and the score of an i-symbol indel is -(q + r * i). Each match
    automatically receives score 10. This simple scoring scheme may be
    used for DNA sequences. NOTE: all scores are integers.

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

    Here the -22 at position (3,2) is the score of replastd::cing N by R.
    This general scoring scheme is useful for protein sequences where the
    set of protein characters and Dayhoff matrix are specified in the file.

    The SIM program is written in C and runs under Unix systems on
    Sun workstations and under DOS systems on PCs.
    We think that the program is portable to many machines.

    Sequences to be analyzed are stored in separate files.
    An input file contains all characters of a sequence, separated by
    newline characters, in linear order. No other characters are allowed.
    Since upper case and lower case characters are different, use the same
    case consistently. A sample sequence file of 4 lines is shown below.

GAATTCTAATCTCCCTCTCAACCCTACAGTCACCCATTTGGTATATTAAA
GATGTGTTGTCTACTGTCTAGTATCCCTCAAGTAGTGTCAGGAATTAGTC
ATTTAAATAGTCTGCAAGCCAGGAGTGGTGGCTCATGTCTGTAATTCCAG
CACTGGAGAGGTAGAAGTG

    The SIM program generates alignments to the standard output.
    Redirection is required to put alignments in a specified file. 
    A sample output file is shown below, where alignment 1 aligns
    the 3240-3313 region of A with the 26300-26374 region of B.

Match   Mismatch   Gap-Open Penalty   Gap-Extension Penalty
 10      -15            30                 5

                 Upper Sequence : A
                         Length : 36741
                 Lower Sequence : B
                         Length : 73360

*********************************************************
      Number 1 Local Alignment
      Similarity Score : 155
      Match Percentage : 70%
      Number of Matches : 54
      Number of Mismatches : 18
      Total Length of Gaps : 5
      Begins at (3240, 26300) and Ends at (3313, 26374)

    0     .    :    .    :    .    :    .    :    .    :
 3240 GAGGGGCATTTGAGGGTGTTTCCAATGTTCCTGTTATTCGGAATAGCGCT
      || || ||||||-|| || ||||||-||   || ||||  |||||| || 
26300 GATGGCCATTTG GGTTGGTTCCAA GTCTTTGGTATTGTGAATAGTGCC

   50     .    :    .    :    .
 3290 GGTGTGAACATTC   TGCACAGGTCT
      |   | ||||| |---|||||| ||||
26348 GCAATAAACATACGTGTGCACATGTCT

*********************************************************
      Number 2 Local Alignment
      Similarity Score : 145
      Match Percentage : 86%
      Number of Matches : 19
      Number of Mismatches : 3
      Total Length of Gaps : 0
      Begins at (4566, 47235) and Ends at (4587, 47256)

    0     .    :    .    :
 4566 GGGGAGGGAGGGGAGCAGAGCA
      |||||||  |||||| ||||||
47235 GGGGAGGTTGGGGAGAAGAGCA


    To find k best non-intersecting alignments of segments from two
    sequences in files A and B, respectively, use a command of form

	   sim  k  A  B  ms  q  r > result

    where sim is the name of the object code, k is a positive integer,
    ms is a negative integer specifying mismatch weight, q and r are
    non-negative integers specifying gap-open and gap-extend penalties,
    respectively. Output alignments are saved in the file "result".

    For using a scoring matrix defined in file S, use a command of form

	   sim  k  A  B  S  q  r > result

    Note that ms is replaced by the file S.

    To find repeats in one sequence in file A, use a command of form

	   sim  k  A  ms  q  r > result
    or
	   sim  k  A  S  q  r > result

    Acknowledgments
    The functions diff() and display() are from Gene Myers. We made
    the following modifications: similarity weights (integer), instead of
    distance weights (float), are used, the aligned pairs already output
    are not used in the subsequent computation, and the positions of
    sequence characters in output alignments are shown by increment of 50.
    T. Mark Reboul found truncation errors in the previous version
    of SIM and also made a few helpful suggestions.
*/

#include   <stdlib.h>
#include   <string.h>
#include <iostream>
#include <fstream>
#include   <SDGString.h>
#include   <SDGMemBioSeq.h>

class Lalign
{
  bool no_prev_align;
  SDGString name1,name2;           /* names of sequence files    */
  char *A,  *B;			       /* Storing two sequences      */
  int  M, N;			       /* Sequence lengths           */
  int  nseq;			       /* Number of sequences        */

  int v[128][128];		       /* substitution scores */
  int q,r,qr;			       /* qr = q + r */

  int  *S;			       /* saving operations for diff */
  int *CC, *DD;			       /* saving matrix scores */
  int *RR, *SS, *EE, *FF;	       /* saving start-points */
  int *HH, *WW;		 	       /* saving matrix scores */
  int *II, *JJ, *XX, *YY;	       /* saving start-points */
  int  m1, mm, n1, nn;		       /* boundaries of recomputed area */
  int  rl, cl;			       /* left and top boundaries */
  int  min;			       /* minimum score in LIST */
  bool flag;			       /* indicate if recomputation necessary*/

  int *sapp;				/* Current script append ptr */
  int  last;				/* Last script op appended */

  int max_align;                        //maximum number of alignment

  /* for saving used aligned pairs */
  struct ONE 
  {
    int COL ;
    struct ONE*  NEXT ;
  };
  typedef ONE pair;
  typedef ONE *pairptr;
  pairptr *row, z; 			
  short tt;

  // Structures for saving best scores
  typedef struct NODE
  { int  SCORE;
    int  STARI;
    int  STARJ;
    int  ENDI;
    int  ENDJ;
    int  TOP;
    int  BOT;
    int  LEFT;
    int  RIGHT; }  vertex, *vertexptr;
		
  vertexptr  *LIST;			/* an array for saving k best scores */
  vertexptr  low;			/* lowest score node in LIST */
  vertexptr  most;			/* latestly accessed node in LIST */
  int numnode;			        /* the number of nodes in LIST */

  //current alignement data
  int count;			        /* number of alignment found */	
  int endi, endj, stari, starj;	        /* endpoint and startpoint */ 
  int score;   			        /* the max score in LIST */
  int I, J;				/* current positions of A ,B */
  int no_mat; 				/* number of matches */ 
  int no_mis; 				/* number of mismatches */ 
  int al_len; 				/* length of alignment */


/* k-symbol indel score */
  int gap(int k) 
    {
      if(k <= 0)
	return 0;
      else
	return q+r*(k);
    };	

/* DIAG() assigns value to x if (ii,jj) is never used before */
  void DIAG(int ii, int jj, int& x, int value)				
    { 
      for ( tt = 1, z = row[(ii)]; z != 0; z = z->NEXT )	
	if ( z->COL == (jj) )				
	  { tt = 0; break; }				
      if ( tt )						
	x = ( value );					
    };

/* replace (ss1, xx1, yy1) by (ss2, xx2, yy2) if the latter is large */
  void ORDER(int& ss1, int& xx1, int& yy1, int& ss2, int& xx2, int& yy2)
    { 
      if ( ss1 < ss2 )					
	{ ss1 = ss2; xx1 = xx2; yy1 = yy2; }		
      else							
	if ( ss1 == ss2 )					
	  { 
	    if ( xx1 < xx2 )				
	      { xx1 = xx2; yy1 = yy2; }			
	    else						
	      if ( xx1 == xx2 && yy1 < yy2 )		
		yy1 = yy2;					
	  }							
    };

/* Append "Delete k" op */
  void DEL(int k)				
    { 
      I += k;				
      al_len += k;				
      if (last < 0)				
	last = sapp[-1] -= (k);		
      else					
	last = *sapp++ = -(k);		
    };

  /* Append "Insert k" op */
  void INS(int k)				
    {
      J += k;				
      al_len += k;				
      if (last < 0)				
	{ sapp[-1] = (k); *sapp++ = last; }	
      else					
	last = *sapp++ = (k);		
    };
  
  /* Append "Replace" op */
  void REP() 				
    { 
      last = *sapp++ = 0; 			
      al_len += 1;				
    };

/* The following definitions are for function diff() */

/* A big pass to compute K best classes */
void big_pass();

/* Determine the left and top boundaries of the recomputed area */
void locate();

/* recompute the area on forward pass */
void small_pass(int count);

/* Add a new node into list.  */
int addnode(int c, int ci, int cj, int i, int j, int K, int cost);

/* Find and remove_self_hits the largest score in list */
 vertexptr findmax();

/* return 1 if no node in LIST share vertices with the area */
 int no_cross();

/* diff(A,B,M,N,tb,te) returns the score of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */
 int diff(char A[], char B[], int M, int N, int tb, int te);

  /* Alignment display routine */
 void display(char A[], char B[], int M, int N, int S[],int AP,int BP);

  /* Alignment return routine */
 void getAlignedStr(char A[],char B[], int M, int N, int S[],SDGString& s1, SDGString& s2);


 public:

 Lalign(int nb_align=10): max_align(nb_align+1)
   {
     low = 0;
     most = 0;	

     nseq=0;

     setMismatch(1,3);
     setGap(5,2);
     
     A=NULL;
     B=NULL;
     N=0;
     M=0;
     no_prev_align=true;
   };

 Lalign(const Lalign& a)
   {
     q=a.q;
     r=a.r;
     qr=a.qr;

     nseq=a.nseq;

     max_align=a.max_align;
     score=a.score;

     no_mat=a.no_mat;
     no_mis=a.no_mis;
     al_len=a.al_len;
     
     sapp=S;
     last=a.last;

     
     for ( int i = 0; i < 128 ; i++ )
       for ( int j = 0; j < 128 ; j++ )
	 v[i][j] = a.v[i][j];
     
     
     name1=a.name1;
     name2=a.name2;

     M=a.M;
     A=new char[M+2];
     memcpy(A,a.A,(M+2)*sizeof(char));
     N=a.N;
     B=new char[N+2];
     memcpy(B,a.B,(N+2)*sizeof(char));

     int j = (N + 1) * sizeof(int);
     CC = new int[j];
     memcpy(CC,a.CC,j);
     DD = new int[j];
     memcpy(DD,a.DD,j);
     RR = new int[j];
     memcpy(RR,a.RR,j);
     SS = new int[j];
     memcpy(SS,a.SS,j);
     EE = new int[j];
     memcpy(EE,a.EE,j);
     FF = new int[j];
     memcpy(FF,a.FF,j);

     int i = (M + 1) * sizeof(int);
     S = new int[i + j];
     memcpy(S,a.S,i+j);


     HH = new int[i];
     memcpy(HH,a.HH,i);
     WW = new int[i];
     memcpy(WW,a.WW,i);
     II = new int[i];
     memcpy(II,a.II,i);
     JJ = new int[i];
     memcpy(JJ,a.JJ,i);
     XX = new int[i];
     memcpy(XX,a.XX,i);
     YY = new int[i];
     memcpy(YY,a.YY,i);

     row = new pairptr[(M + 1) * sizeof(pairptr)];
     memcpy(row,a.row,(M + 1) * sizeof(pairptr));
     for ( i = 1; i <= M; i++ )
       if ( nseq == 2 )
	 row[i] = a.row[i];
       else
	 { 
	   row[i] = z = new pair[(int)sizeof(pair)];
	   z->COL = a.row[i]->COL;			
	   z->NEXT = a.row[i]->NEXT;
	 }


     LIST = new vertexptr[max_align * sizeof(vertexptr)];
     for ( i = 0; i < max_align ; i++ )
       {
	  LIST[i] = new vertex[(int) sizeof(vertex)];
	  memcpy(LIST[i],a.LIST[i],(int) sizeof(vertex));
       }  

  m1=a.m1;
  mm=a.mm;
  n1=a.n1;
  nn=a.nn;	
  rl=a.rl;
  cl=a.cl;	
  min=a.min;	
  flag=a.flag;	

  max_align=a.max_align;

  z=a.z; 			
  tt=a.tt;

  low=a.low;
  most=a.most;
  numnode=a.numnode;

  count=a.count;
  endi=a.endi;
  endj=a.endj;
  stari=a.stari;
  starj=a.starj;
  I=a.I;
  J=a.J;

   };

 virtual ~Lalign()
   {
     if(M!=0)
       {
	 delete[] A;
	 M=0;
       }	 
     if(N!=0)
       {
	 delete[] B;
	 N=0;
       }	 
     reset_align();
   };

 void reset_align(void)
   {
     if(no_prev_align) return;

     low = 0;
     most = 0;	
     score = 0;

     delete[] CC;
     delete[] DD;
     delete[] RR;
     delete[] SS;
     delete[] EE;
     delete[] FF;
     delete[] HH;
     delete[] WW;
     delete[] II;
     delete[] JJ;
     delete[] XX;
     delete[] YY;
     delete[] S;
	
     for (int i = 1; i <= M; i++ )
       if ( nseq != 2 )
	 { 
	   delete row[i];
	 }
     delete[] row;
     
     for (int i = 0; i < max_align ; i++ )
       delete LIST[i];
     delete[] LIST;
   };

 void setMismatch(int match, int mismatch)
   {
     /* set match and mismatch weights */
     for ( int i = 0; i < 128 ; i++ )
       for ( int j = 0; j < 128 ; j++ )
	 if (((i == j) || (i=='N' && j!='N') || (i!='N' && j=='N'))
	     && (i!='X' || j!='X'))
	     v[i][j] = match;
	   else
	     v[i][j] = -mismatch;
   };

 void setMismatch(SDGString filename)
   {
     std::ifstream fin(filename);

     char  alph[129];		    

     fin.getline(alph,256);
     int size = strlen(alph);
     int match = 0;
     int mismh = 0;
     for ( int i = 0; i < size ; i++ )
       for ( int j = 0; j <= i ; j++ )
	 { 
	   int ms;
	   fin>>ms;
	   v[(unsigned)alph[i]][(unsigned)alph[j]] = v[(unsigned)alph[j]][(unsigned)alph[i]] = ms;
	   if ( ms > match ) match = ms;
	   if ( ms < mismh ) mismh = ms;
	 }
   };

 void setGap(int gapopen, int gapextend)
   {
     q=gapopen;
     r=gapextend;
     qr=q+r;
   };

 void setSeq(SDGBioSeq s,int num=1)
   {
     if(num==1)
       {
	 if(M!=0)
	   {
	     delete[] A;
	     M=0;
	   }	 
	 name1=s.getDE();
	 M=s.length();
	 A=new char[M+2];
	 strcpy(&A[1],s.toString().start());
       }
     else
       {
	 if(N!=0)
	   {
	     delete[] B;
	     N=0;
	   }	 
	 name2=s.getDE();
	 N=s.length();
	 B=new char[N+2];     
	 strcpy(&B[1],s.toString().start());
       }
     if(nseq<2) nseq++;
     if(nseq==1)
       {
	 if(N!=0)
	   {
	     delete[] B;
	     N=0;
	   }	 
	 name2=name1;
	 N=M;
	 B=new char[N+2];;     
	 strcpy(B,A);
       }
   };

 void setSeq(SDGBioSeq s1,SDGBioSeq s2)
   {
     setSeq(s1,1);
     setSeq(s2,2);
   };

 void setNbSeq(int n)
   {
     nseq=n;
   };

/* align() calculate n best non-intersecting alignments of
   the segments of A and B in order of similarity scores, where
   v[a][b] is the score of aligning a and b, and -(Q+R*i) is the score
   of an i-symbol indel.  						*/
  void align(void);
  void findNext(void);
  void view(void);
  void write_map(std::ostream& out=std::cout);
  //  void getAlignedSeq(SDGAlignedBioSeq_hdl& seq1,
  //			      SDGAlignedBioSeq_hdl& seq2 );
  void getAlignedSeq(SDGString& seq1,
		     SDGString& seq2 );

  SDGString getNameSeq1(){return name1;};
  SDGString getNameSeq2(){return name2;};
  
  int getStartSeq1(){return stari;};
  int getStartSeq2(){return starj;};
  int getEndSeq1(){return endi;};
  int getEndSeq2(){return endj;};
  long getScore(){return score;};
  //void setScore(int score){this.score = score;};
  //  double z_score(int nbrep);
  long getAlignmentLength(){return al_len;};
  long getNumberMismatch(){return no_mis;};
  long getNumberMatch(){return no_mat;};
  double getIdentity(){return (double)no_mat/(no_mat+no_mis);};
};

















