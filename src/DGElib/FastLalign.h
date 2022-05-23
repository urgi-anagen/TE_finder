#ifndef FASTLALIGN_H
#define FASTLALIGN_H

#include   <limits>
#include   <stdio.h>
#include   <string.h>
#include   <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include   <vector>
#include   <deque>
#include   "CountStat.h"
#include   <SDGString.h>
#include   <SDGMemBioSeq.h>
#include   <BioSeq.h>

class FastLalign
{

  SDGString name1,name2;    /* names of sequence files    */
  char  *A,  *B;			/* Storing two sequences      */
  long  M, N;			/* Sequence lengths     */
  int *path;

  int match, mismh;	       /* max and min substitution weights */
  int v[128][128];	       /* substitution scores */
  int  q, r;                   /* gap penalties */
  int  gaplen;                 /* minimum length for constant-cost insertion */
  int  long_gap_pen;

  int endi,endj,starti,startj;
  int score;   		         /* the max score */
  int len_ali;

  double identity;

  std::vector<int> alignment;
 /* k-symbol indel score */

  int whosmax;
  int search_max(int a,int b)
    {
      if(a>=b)
	{
	  whosmax=1;
	  return a;
	}
      whosmax=2;
      return b;
    };

  int search_max(int a,int b,int c)
    {
      if(a>=b && a>=c)
	{
	  whosmax=1;
	  return a;
	}
      if(b>a && b>=c)
	{
	  whosmax=2;
	  return b;
	}
      whosmax=3;
      return c;
    };

  int search_max(int a,int b,int c,int d)
    {
      if(a>=b && a>=c && a>=d)
	{
	  whosmax=1;
	  return a;
	}
      if(b>a && b>=c && b>=d)
	{
	  whosmax=2;
	  return b;
	}
      if(c>a && c>b && c>=d)
	{
	  whosmax=3;
	  return c;
	}
      whosmax=4;
      return d;
    };

public:

 FastLalign()
   {
     setMismatch(1,3);
     setGap(5,2,50);
     A=B=nullptr;
     M=N=0;
     path=nullptr;
   };

 FastLalign(const FastLalign& a)
   {
     setMismatch(a.match,a.mismh);
     setGap(a.q,a.r,a.gaplen);
     A=B=nullptr;
     M=N=0;
     path=nullptr;

     name1=a.name1;
     name2=a.name2;
     
     if(M!=0)
       {     
	 M=a.M;
	 A=new char[M+2];
	 strcpy(A,a.A);
       }

     if(N!=0)
       {     
	 N=a.N; 
	 B=new char[N+2];
	 strcpy(B,a.B);
       }

     if(a.path!=nullptr)
       {
	 path=new int[(N+2)*(M+2)];
	 memcpy(path,a.path,(N+2)*(M+2)*sizeof(int));
       }
     match=a.match;
     mismh=a.mismh;
     q=a.q;
     r=a.r;
     gaplen=a.gaplen;
     long_gap_pen=a.long_gap_pen;
     score=a.score;
     
     alignment=a.alignment;
     for ( int i = 0; i < 128 ; i++ )
       for ( int j = 0; j < 128 ; j++ )
	     v[i][j] = a.v[i][j];

   };

 virtual ~FastLalign()
   {
     reset_align();
   };

 void reset_align(void)
   {
     alignment.clear();
     if(A!=nullptr)
       delete[] A;
     if(B!=nullptr)
       delete[] B;
     if(path!=nullptr)
       delete[] path;
     len_ali=0;
     score=0;
     identity=0;
     starti=startj=endi=endj=0;
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

 void setGap(int gapopen, int gapextend, int glength=0)
   {
     q=gapopen;
     r=gapextend;
     gaplen=glength;
     long_gap_pen=q+gaplen*r;
   };

 void setSeq(SDGBioSeq s1, SDGBioSeq s2)
   {
     reset_align();
     name1=s1.getDE();
     M=s1.length();
     A=new char[M+2];
     strcpy(&A[1],s1.toString().start());

     name2=s2.getDE();
     N=s2.length();
     B=new char[N+2];     
     strcpy(&B[1],s2.toString().start());

     path=new int[(N+2)*(M+2)];
   };

    void setSeq(BioSeq s1, BioSeq s2)
    {
        reset_align();
        name1=s1.header;
        M=s1.size();
        A=new char[M+2];
        strcpy(&A[1],s1.c_str());

        name2=s2.header;
        N=s2.size();
        B=new char[N+2];
        strcpy(&B[1],s2.c_str());

        path=new int[(N+2)*(M+2)];
    };

 void align(void)
   {
     if(gaplen==0)
	 align_pass_no_long_gap();
     else
	 align_pass_long_gap();
     traceback();
     len_ali=alignment.size();
   };

 int score_pass(void)
   {
     if(gaplen==0)
	 return score_pass_no_long_gap();
     else
	 return score_pass_long_gap();
   };

 void align_pass_no_long_gap(void);
 int score_pass_no_long_gap(void);
 void align_pass_long_gap(void);
 int score_pass_long_gap(void);
 void traceback(void);

 void view(void);

 int getStartSeq1(){return starti;};
 int getStartSeq2(){return startj;};
 int getEndSeq1(){return endi;};
 int getEndSeq2(){return endj;};

 int getGapOpen(){return q;};
 int getGapExtend(){return r;};

 SDGString getNameSeq1(){return name1;};
 SDGString getNameSeq2(){return name2;};

 void getAlignedSeq(SDGString& seq1,
		     SDGString& seq2 );

 int getScore(){return score;};
 double getIdentity(void);
 int getLength(void){return len_ali;};
};

#endif








