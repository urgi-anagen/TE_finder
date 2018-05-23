#ifndef FASTEXTALIGN_H
#define FASTEXTALIGN_H

#include   <stdio.h>
#include   <stdlib.h>
#include   <string.h>
#include   <string>
#include   <iostream>
#include   <iomanip>
#include   <sstream>
#include   <deque>
#include   "CountStat.h"
//#include   "Alea.h"
#include   <SDGString.h>
#include   <SDGMemBioSeq.h>

class FastExtAlign
{

  SDGString name1,name2;    /* names of sequence files    */
  char  *A,  *B;			/* Storing two sequences      */
  unsigned long  lenA, lenB;			/* Sequence lengths     */
  unsigned long  M, N;			/* Sequence lengths     */


  int match, mismh;	       /* max and min substitution weights */
  int v[128][128];	       /* substitution scores */
  int  q, r;                   /* gap penalties */

  unsigned endi,endj,starti,startj;
  int score;   		         /* the max score */

  int whosmax;
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

public:

 FastExtAlign()
   {
     setMismatch(1,2);
     setGap(3,1);
     A=B=NULL;
     M=N=0;
     lenA=lenB=0;
   };

 virtual ~FastExtAlign()
   {
     reset_align();
   };

 void reset_align(void)
   {
     if(A!=NULL)
       delete[] A;
     if(B!=NULL)
       delete[] B;
     score=0;
     starti=startj=endi=endj=0;
     N=M=lenA=lenB=0;
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

 void setGap(int gapopen, int gapextend)
   {
     q=gapopen;
     r=gapextend;
   };

 void setSeq(const SDGBioSeq& s1,const SDGBioSeq& s2)
   {
     reset_align();
     lenA=s1.length();
     A=new char[lenA+2];
     strcpy(&A[1],s1.toString().start());

     lenB=s2.length();
     B=new char[lenB+2];     
     strcpy(&B[1],s2.toString().start());
   };

 void setSeq1(const SDGBioSeq& s1)
   {
     if(A!=NULL)
       delete[] A;
     lenA=s1.length();
     A=new char[lenA+2];
     strcpy(&A[1],s1.toString().start());
   };

 void setSeq2(const SDGBioSeq& s2)
   {
     if(B!=NULL)
       delete[] B;
     lenB=s2.length();
     B=new char[lenB+2];     
     strcpy(&B[1],s2.toString().start());
   };

 void setStart(unsigned i, unsigned j, unsigned l)
   {
     starti=i;
     startj=j;
     if(l>=lenA || l>=lenB)
       throw SDGException(NULL,"FastExtAlign.setStart: extension length too long",-1);
     M=N=l;
   };

 int extend_dir(int score);
 int extend_rev(int score);

 int getStartSeq1(void){return starti;};
 int getStartSeq2(void){return startj;};
 int getEndSeq1(void){return endi;};
 int getEndSeq2(void){return endj;};

 long getExtLenSeq1(void){return std::abs((long)starti-endi);};
 long getExtLenSeq2(void){return std::abs((long)startj-endj);};
     
 int getGapOpen(void){return q;};
 int getGapExtend(void){return r;};

 SDGString getNameSeq1(){return name1;};
 SDGString getNameSeq2(){return name2;};
  

};

#endif








