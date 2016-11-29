//#include "Alea.h"
#include <stdio.h>
#include "CountStat.h"
#include "Galign.h"

//--------------------------------------------------------------------------
void Galign::align(void)
{ 
  reset_align();
  change = false;
  if ( M > N )
    {
      change=true;
      char* T=A;
      A=B;
      B=T;
      int L=M;
      M=N;
      N=L;
    }
  int  i, j;			/* row and column indices */
	
  /* allocate space for all vectors */
  j = (N + 1);
  CC = new int[j];
  DD = new int[j];
  RR = new int[j];
  SS = new int[j];
  i = (M + 1);
  S = new int[i + j];
  
  sapp = S;
  last = 0;
  al_len = 0;
  no_mat = 0;
  no_mis = 0;
  score = diff2(A,B,M,N,q,q,0,0,0,0);
  no_prev_align=false;
};
//--------------------------------------------------------------------------
void Galign::view(void)
{
   printf("Max Match   Min Mismatch   Gap-Open Penalty   Gap-Extension Penalty\n");
   printf("   %d          %d              %d                  %d\n\n", match, mismh, q, r);

  if ( change )
    {
      std::cout<<"                 Upper Sequence : "<< name2<<std::endl;
      printf("                         Length : %ld\n", M);
    }
  else
    {
      std::cout<<"                 Upper Sequence : "<< name1<<std::endl;
      printf("                         Length : %ld\n", M);
    }

  if ( change )
    {
      std::cout<<"                 Lower Sequence : "<<name1<<std::endl;
      printf("                         Length : %ld\n\n", N);
    }
  else
    {
      std::cout<<"                 Lower Sequence : "<<name2<<std::endl;
      printf("                         Length : %ld\n\n", N);
    }  
  /* Output the best alignment */
   printf("      Best Alignment with Constant Cost for Long Insertions\n");
   printf("      Similarity Score : %ld\n",score);
   printf("      Match Percentage : %ld%%\n", (100*no_mat)/al_len);
   printf("      Number of Matches : %ld\n", no_mat);
   printf("      Number of Mismatches : %ld\n", no_mis);
   printf("      Total Length of Gaps : %ld\n", al_len-no_mat-no_mis);
   printf("      Note that terminal gaps are not penalized\n");
   display(A,B,M,N,S,1,1);
};

//--------------------------------------------------------------------------
/* diff2(A,B,M,N,tb,te,sc,sr,ec,er) returns the score of an optimum conversion
   between A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script. If sc = 0, then
   the beginning deletion is not penalized; if sr = 0, the beginning insertion is
   not penalized; if ec = 0, the ending deletion is not charged; if er = 0;
   then the ending insertion is not charged. Any insertion of length at least
   gaplen is given a constant cost */

int Galign::diff2(char *A,char *B,int M,int N,int tb,int te,int sc,int sr,int ec,int er)
{ int   midi, midj, type;	/* Midpoint, type, and cost */
  int midc;
  int  ss,cc;

{ register int   i, j;
  register int c, e, d, s;
           int t, *va;
	   int  g, temp;

/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) DEL(M);
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
      va = v[(unsigned)A[1]];
      for (j = 1; j <= N; j++)
	{ c = va[(unsigned)B[j]] - ( sr * gap2(j-1) + er * gap2(N-j) );
          if (c > midc)
           { midc = c;
             midj = j;
           }
	}
      if (midj == -1)
        { DEL(1); INS(N); }
      else
      if (midj == 0)
        { INS(N); DEL(1); }
      else
        { if (midj > 1) INS(midj-1);
	REP();
	  if ( A[1] == B[midj] )
	     no_mat += 1;
	  else
	    if (A[1]!='N' && B[midj]!='N') // modif hadi 27/10/06
	      no_mis += 1;
          if (midj < N) INS(N-midj);
        }
      return midc;
    }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;			/* Forward phase:                          */
  CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
  t = - q * sr;
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
      va = v[(unsigned)A[i]];
      for (j = 1; j <= N; j++)
        { if ((c = c - qr) > (e = e - r)) e = c;
	  if ( j == N && !ec )
            { if ((c = CC[j] ) > (d = DD[j] )) d = c;}
	  else
            if ((c = CC[j] - qr) > (d = DD[j] - r)) d = c;
	  c = s+va[(unsigned)B[j]];
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
      va = v[(unsigned)A[i+1]];
      for (j = N-1; j >= 0; j--)
        { if ((c = c - qr) > (e = e - r)) e = c;
	  if ( !j && !sc )
            { if ((c = RR[j] ) > (d = SS[j] )) d = c;}
	  else
            if ((c = RR[j] - qr) > (d = SS[j] - r)) d = c;
	  c =  s+va[(unsigned)B[j+1]];
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
      if (c > midc || (CC[j] != DD[j] && RR[j] == SS[j]))
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
    {  diff2(A,B,midi,midj,tb,q,sc,sr,cc,1);
       diff2(A+midi,B+midj,M-midi,N-midj,q,te,ss,1,ec,er);
    }
  else
    {  diff2(A,B,midi-1,midj,tb,0,sc,sr,cc,1);
      DEL(2);
       diff2(A+midi+1,B+midj,M-midi-1,N-midj,0,te,ss,1,ec,er);
    }
  return midc;
}
//--------------------------------------------------------------------------
/* Alignment display routine */
void Galign::display(char A[], char B[], int M, int N, int S[], int AP,int BP)
{
  /* Alignment display routine */
  char ALINE[51], BLINE[51], CLINE[51];

  register char *a, *b, *c;
  register int   i,  j, op;
           int   lines, ap, bp;

  i = j = op = lines = 0;
  ap = AP;
  bp = BP;
  a = ALINE;
  b = BLINE;
  c = CLINE;
  while (i < M || j < N)
    { 
      if (op == 0 && *S == 0)
        { 
	  op = *S++;
          *a = A[++i];
          *b = B[++j];
          *c++ = (*a++ == *b++) ? '|' : ' ';
        }
      else
        { if (op == 0)
            op = *S++;
          if (op > 0)
            { 
	      *a++ = '-';
              *b++ = B[++j];
              op--;
            }
          else
            { 
	      *a++ = A[++i];
              *b++ = '-';
              op++;
            }
          *c++ = ' ';
        }
      if (a >= ALINE+50 || (i >= M && j >= N))
        { 
	  *a = *b = *c = '\0';
           printf("\n%5d ",50*lines++);
          for (b = ALINE+10; b <= a; b += 10)
             printf("    .    :");
          if (b <= a+5)
             printf("    .");
           printf("\n%5d %s\n      %s\n%5d %s\n",ap,ALINE,CLINE,bp,BLINE);

	  ap = AP + i;
	  bp = BP + j;
          a = ALINE;
          b = BLINE;
          c = CLINE;
        }
    }
}
//--------------------------------------------------------------------------
void Galign::getAlignedSeq(SDGString& seq1, SDGString& seq2)
{
  getAlignedStr(A,B,M,N,S);
  if(change)
    {
      seq1=alignedStr2;
      seq2=alignedStr1;
    }
  else
    {
      seq1=alignedStr1;
      seq2=alignedStr2;
    }
};
//--------------------------------------------------------------------------
// double Galign::z_score(int nbrep)
// {
//   CountStat stat;
//   char *tmp=A;
  
//   int score_ori=getScore();
//   A=new char[M+2];

//   for(int rep=0;rep<nbrep;rep++)
//     {
//       for(unsigned idx=1;idx<M;idx++)
// 	A[idx]=tmp[alea.rndLoHi(1,M)];
//       align();
//       stat.add(getScore());
//     }
//   delete [] A;

//   A=tmp;

//   return (score_ori-stat.mean())/stat.sd();
// };
//--------------------------------------------------------------------------
/* Alignment return routine */
void Galign::getAlignedStr(char A[], char B[], int M, int N, int S[])
{
  /* Alignment display routine */
  char ALINE[51], BLINE[51];

  alignedStr1="";
  alignedStr2="";

  register char *a, *b;
  register int   i,  j, op;

  i = j = op = 0;
  a = ALINE;
  b = BLINE;
  while (i < M || j < N)
    { 
      if (op == 0 && *S == 0)
        { 
		  op = *S++;
		  *a = A[++i];
		  *b = B[++j];
		  *a++; *b++; // gcc says value not used! true?
        }
      else
        { if (op == 0)
            op = *S++;
          if (op > 0)
            { 
        	  *a++ = '-';
              *b++ = B[++j];
              op--;
            }
          else
            { 
        	  *a++ = A[++i];
              *b++ = '-';
              op++;
            }
        }
      if (a >= ALINE+50 || (i >= M && j >= N))
        { 
		  *a = *b = '\0';

		  alignedStr1+=ALINE;
		  alignedStr2+=BLINE;

          a = ALINE;
          b = BLINE;
        }
    }
}

