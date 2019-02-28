/*  A GLOBAL ALIGNMENT PROGRAM (GAP):

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
     be appreciated: "On global sequence alignment" (CABIOS).
	      Xiaoqiu Huang
	      Department of Computer Science
	      Michigan Technological University
	      Houghton, MI 49931
              E-mail: huang@cs.mtu.edu

    The GAP program computes a global alignment of two sequences
    without penalizing terminal gaps. It delivers the alignment in
    linear space, so long sequences can be aligned. 

    Users supply scoring parameters. In the simplest form, users just
    provide 3 integers: ms, q and r, where ms is the score of a mismatch
    and the score of an i-symbol indel is -(q + r * i). Each match
    automatically receives score 10. This simple scoring scheme may be
    used for DNA sequences. NOTE: all scores are integers.

    In general, users can define an alphabet of characters appearing
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

    The GAP program is written in C and runs under Unix systems on
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

    To find the best alignment of two sequences in files A and B,
    use a command of form

	   gap  A  B  gs  ms  q  r > result

    where gap is the name of the object code, gs is the minimum length
    of any gap in the short sequence receiving a constant gap penalty,
    ms is a negative integer specifying mismatch weight, q and r are
    non-negative integers specifying gap-open and gap-extend penalties,
    respectively. Output alignment is saved in the file "result".

    For using a scoring matrix defined in file S, use a command of form

	   gap  A  B  gs  S  q  r > result

    Note that ms is replaced by the file S.

    Acknowledgments
    The functions diff2() and display() evolved from those written by Gene Myers.
    We made the following modifications: similarity weights (integer), instead of
    distance weights (float), are used, terminal gaps are not penalized, and
    any gap of length at least gs in the short sequence is given a constant
    penalty.
*/

#include   <string.h>
#include   <stdlib.h>
#include <iostream>
#include <fstream>
#include   <SDGString.h>
#include   <SDGMemBioSeq.h>

class Galign
{

 bool no_prev_align;
 SDGString name1,name2;    /* names of sequence files    */
 char  *A,  *B;			/* Storing two sequences      */
 unsigned long  M, N;			/* Sequence lengths     */
 SDGString alignedStr1;
 SDGString alignedStr2;

 int match, mismh;	       /* max and min substitution weights */
 int v[128][128];	       /* substitution scores */
 int  q, r;                    /* gap penalties */
 int  qr;                      /* qr = q + r */
 int  gaplen;                  /* minimum length for constant-cost insertion */
 int  pay;	               /* constant-cost for long insertion */

 bool  change;	        

 int *CC, *DD;			/* saving matrix scores */
 int *RR, *SS;		 	/* saving start-points */
 int  *S;			/* saving operations for diff */


 int *sapp;				/* Current script append ptr */
 int  last;				/* Last script op appended */

 unsigned long no_mat; 				/* number of matches */ 
 unsigned long no_mis; 				/* number of mismatches */ 
 unsigned long al_len; 				/* length of alignment */

 unsigned short nseq;

 long score;   		         /* the max score */


/* The following definitions are for function diff() */

int  diff2(char *A,char *B,int M,int N,int tb,int te,int sc,int sr,int ec,int er);

/* k-symbol indel score */
 int gap(int k) 
   {
     if(k <= 0)
       return 0;
     else
       return q+r*(k);
   };	

  int gap2(int k) 
    {
      if(k <= 0)
	return 0;
      else
	return ((k) <= gaplen ? q+r*(k) : pay);
    };	

/* Append "Delete k" op */
  void DEL(int k)				
    { 
      al_len += k;				
      if (last < 0)				
	last = sapp[-1] -= (k);		
      else					
	last = *sapp++ = -(k);		
    };

/* Append "Insert k" op */
  void INS(int k)				
    {
      al_len += k;				
      if (last > 0)				
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

 void display(char A[], char B[], int M, int N, int S[], int AP,int BP);
 void getAlignedStr(char A[], char B[], int M, int N, int S[]);

public:

 Galign()
   {
     change=false;
     setMismatch(1,3);
     setGap(5,2,50);
     no_prev_align=true;
   };

 Galign(const Galign& a)
   {
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
     int i = (M + 1) * sizeof(int);
     S = new int[i + j];
     memcpy(S,a.S,i+j);
  
     match=a.match;
     mismh=a.mismh;
     q=a.q;
     r=a.r;
     qr=a.qr;
     gaplen=a.gaplen;
     pay=a.pay;
     change=a.change;
     nseq=a.nseq;
     score=a.score;
     last=a.last;
     no_mat=a.no_mat;
     no_mis=a.no_mis;
     al_len=a.al_len;

     sapp=S;

     for ( int i = 0; i < 128 ; i++ )
       for ( int j = 0; j < 128 ; j++ )
	     v[i][j] = a.v[i][j];

   };

 virtual ~Galign()
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
     delete[] CC;
     delete[] DD;
     delete[] RR;
     delete[] SS;
     delete[] S;
     alignedStr1="";
     alignedStr2="";
   }

 void setMismatch(int mtch, int msmtch)
   {
     /* set match and mismatch weights */
     for ( int i = 0; i < 128 ; i++ )
       for ( int j = 0; j < 128 ; j++ )
	 if ((i == j || (i=='N' && j!='N') || (i!='N' && j=='N')) 
	     && (i!='X' || j!='X'))
	     v[i][j] = mtch;
	   else
	     v[i][j] = -msmtch;

     match = mtch;
     mismh = -msmtch;
   };

 void setMismatch(SDGString filename)
   {
     std::ifstream fin(filename);
     char  alph[129];
     fin.getline(alph,256);
     int size = strlen(alph);
     match = 0;
     mismh = 0;
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

 void setGap(int gapopen, int gapextend, int glength)
   {
     q=gapopen;
     r=gapextend;
     gaplen=glength;
     pay = q + r * gaplen;
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
   };

 void setSeq(SDGBioSeq s1,SDGBioSeq s2)
   {
     setSeq(s1,1);
     setSeq(s2,2);
   };

 void setNbSeq(unsigned short n)
   {
     nseq=n;
   };

/* align() calculate n best non-intersecting alignments of
   the segments of A and B in order of similarity scores, where
   v[a][b] is the score of aligning a and b, and -(Q+R*i) is the score
   of an i-symbol indel.  						*/
  void align(void);
  void view(void);
  //  double z_score(int nbrep);

  SDGString getNameSeq1(){return name1;};
  SDGString getNameSeq2(){return name2;};
  
  //  void getAlignedSeq(SDGAlignedBioSeq_hdl& seq1,
  //			      SDGAlignedBioSeq_hdl& seq2 );
  void getAlignedSeq(SDGString& seq1,
		     SDGString& seq2 );

  int getScore(){return score;};
  unsigned long getAlignmentLength(){return al_len;};
  unsigned long getNumberMismatch(){return no_mis;};
  unsigned long getNumberMatch(){return no_mat;};


};





