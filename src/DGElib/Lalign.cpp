#include <stdio.h>
#include <string.h>
//#include "Alea.h"
#include "CountStat.h"
#include "Lalign.h"

//--------------------------------------------------------------------------
void Lalign::align()
{ 
  int  i, j;			/* row and column indices */
  
  reset_align();

  /* allocate space for all vectors */
  j = (N + 1);
  CC = new int[j];
  DD = new int[j];
  RR = new int[j];
  SS = new int[j];
  EE = new int[j];
  FF = new int[j];
  i = (M + 1);
  HH = new int[i];
  WW = new int[i];
  II = new int[i];
  JJ = new int[i];
  XX = new int[i];
  YY = new int[i];
  S = new int[i + j];
  row = new pairptr[(M + 1) * sizeof(pairptr)];
  
  /* set up list for each row */
  for ( i = 1; i <= M; i++ )
    if ( nseq == 2 )
      row[i] = 0;
    else
      { 
	row[i] = z = new pair[(int)sizeof(pair)];
	z->COL = i;			
	z->NEXT = 0;
      }
  
  
  LIST = new vertexptr[max_align * sizeof(vertexptr)];
  for ( i = 0; i < max_align ; i++ )
    LIST[i] = new vertex[(int) sizeof(vertex)];
  
  numnode = min = 0;
  big_pass();
  count=max_align;
  no_prev_align=false;
  findNext();

};
//--------------------------------------------------------------------------
void Lalign::findNext()
{
  vertexptr cur; 		  /* temporary pointer to current alignment  */

  if ( numnode == 0 || count<0)
    {
      std::cout<<"There are no more alignments between two sequences"<<std::endl;
      return;
    }

  if ( count && count!=max_align )
    {
      flag = false;
      locate();
      if ( flag )
	small_pass(count);
    }

  cur = findmax();	/* Return a pointer to a node with max score*/
  score = cur->SCORE;
  stari = ++cur->STARI;
  starj = ++cur->STARJ;
  endi = cur->ENDI;
  endj = cur->ENDJ;
  m1 = cur->TOP;
  mm = cur->BOT;
  n1 = cur->LEFT;
  nn = cur->RIGHT;
  rl = endi - stari + 1;
  cl = endj - starj + 1;
  I = stari - 1;
  J = starj - 1;
  sapp = S;
  last = 0;
  al_len = 0;
  no_mat = 0;
  no_mis = 0;
  diff(&A[stari]-1, &B[starj]-1,rl,cl,q,q);
  count--;
};

//--------------------------------------------------------------------------
void Lalign::view(void)
{
  // Show the current best alignment.
  if ( numnode == 0 || count<0)
    {
      std::cout<<"There are no more alignments between two sequences"<<std::endl;
      return;
    }

  /* Output the best alignment */
  (void) printf("\n*********************************************************\n");
  (void) printf("      Number %d Local Alignment\n", max_align - count);
  (void) printf("      Similarity Score : %d\n",score);
  (void) printf("      Match Percentage : %d%%\n", (100*no_mat)/al_len);
  (void) printf("      Number of Matches : %d\n", no_mat);
  (void) printf("      Number of Mismatches : %d\n", no_mis);
  (void) printf("      Total Length of Gaps : %d\n", al_len-no_mat-no_mis);
  (void) printf("      Begins at (%d, %d) and Ends at (%d, %d)\n",
		stari,starj, endi,endj);
  display(&A[stari]-1,&B[starj]-1,rl,cl,S,stari,starj);
  (void) fflush(stdout);
};
//--------------------------------------------------------------------------
void Lalign::write_map(std::ostream& out)
{
  unsigned i=max_align - count;
  out<<"alig"<<i
      <<"\t"<<getNameSeq1()
      <<"\t"<<getStartSeq1()
      <<"\t"<<getEndSeq1()<<std::endl;
  if(nseq==1)
    out<<"alig"<<i<<"\t"<<getNameSeq1();
  else
    out<<"alig"<<i<<"\t"<<getNameSeq2();
  out<<"\t"<<getStartSeq2()
      <<"\t"<<getEndSeq2()<<std::endl;
};

//--------------------------------------------------------------------------
/*
void Lalign::getAlignedSeq(SDGAlignedBioSeq_hdl& seq1, SDGAlignedBioSeq_hdl& seq2)
{
  SDGString s1,s2;
  getAlignedStr(&A[stari]-1,&B[starj]-1,rl,cl,S,s1,s2);
  seq1=new SDGAlignedBioSeq(s1);
  seq2=new SDGAlignedBioSeq(s2);
};
*/

void Lalign::getAlignedSeq(SDGString& seq1, SDGString& seq2)
{
  getAlignedStr(&A[stari]-1,&B[starj]-1,rl,cl,S,seq1,seq2);
};


//--------------------------------------------------------------------------
/* Alignment return routine */
void Lalign::getAlignedStr(char A[], char B[], int M, int N, int S[], SDGString& s1,SDGString& s2)
{
  
  s1="";
  s2="";

  char ALINE[51], BLINE[51];

  char *a, *b;
  int   i,  j, op;

  i = j = op = 0;
  a = ALINE;
  b = BLINE;
  while (i < M || j < N)
    { 
      if (op == 0 && *S == 0)
        { 
	  op = *S++;
          *a++ = A[++i];
          *b++ = B[++j];
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
	  s1+=ALINE;
	  s2+=BLINE;

          a = ALINE;
          b = BLINE;
        }
    }
}
//--------------------------------------------------------------------------
/* Alignment display routine */
void Lalign::display(char A[], char B[], int M, int N, int S[], int AP,int BP)
{
  /* Alignment display routine */
  char ALINE[51], BLINE[51], CLINE[51];

  char *a, *b, *c;
  int   i,  j, op;
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
	      *a++ = ' ';
              *b++ = B[++j];
              op--;
            }
          else
            { 
	      *a++ = A[++i];
              *b++ = ' ';
              op++;
            }
          *c++ = '-';
        }
      if (a >= ALINE+50 || (i >= M && j >= N))
        { 
	  *a = *b = *c = '\0';
          (void) printf("\n%5d ",50*lines++);
          for (b = ALINE+10; b <= a; b += 10)
            (void) printf("    .    :");
          if (b <= a+5)
            (void) printf("    .");
          (void) printf("\n%5d %s\n      %s\n%5d %s\n",ap,ALINE,CLINE,bp,BLINE);
	  ap = AP + i;
	  bp = BP + j;
          a = ALINE;
          b = BLINE;
          c = CLINE;
        }
    }
}
//--------------------------------------------------------------------------
/* A big pass to compute n best classes */
void Lalign::big_pass()
{ 
  int  i, j;			/* row and column indices */
  int  c;			/* best score at current point */
  int  f;			/* best score ending with insertion */
  int  d;			/* best score ending with deletion */
  int  p;			/* best score at (i-1, j-1) */
  int  ci, cj;		/* end-point associated with c */
  int  di, dj;		/* end-point associated with d */
  int  fi, fj;		/* end-point associated with f */
  int  pi, pj;		/* end-point associated with p */
  int  *va;				/* pointer to v(A[i], B[j]) */
	
	/* Compute the matrix and save the top K best scores in LIST
	   CC : the scores of the current row
	   RR and EE : the starting point that leads to score CC
	   DD : the scores of the current row, ending with deletion
	   SS and FF : the starting point that leads to score DD        */
 	/* Initialize the 0 th row */
	for ( j = 1; j <= N ; j++ )
	  {
	     CC[j] = 0;
	     RR[j] = 0;
	     EE[j] = j;
	     DD[j] = - (q);
	     SS[j] = 0;
	     FF[j] = j;
	  }
	for ( i = 1; i <= M; i++) 
	  {  c = 0;				/* Initialize column 0 */
	     f = - (q);
	     ci = fi = i;
	     va = v[(unsigned)A[i]];
	     if ( nseq == 2 )
	       { p = 0;
	         pi = i - 1;
	         cj = fj = pj = 0;
	       }
	     else
	       { p = CC[i];
		 pi = RR[i];
		 pj = EE[i];
	         cj = fj = i;
	       }
	     for ( j = (nseq == 2 ? 1 : (i+1)) ; j <= N ; j++ )  
	       {  
		 f = f - r;
		 c = c - qr;
		 ORDER(f, fi, fj, c, ci, cj);
		 c = CC[j] - qr; 
		 ci = RR[j];
		 cj = EE[j];
		 d = DD[j] - r;
		 di = SS[j];
		 dj = FF[j];
		 ORDER(d, di, dj, c, ci, cj);
		 c = 0;
		 DIAG(i, j, c, p+va[(unsigned)B[j]]);		/* diagonal */
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj);
		  ORDER(c, ci, cj, f, fi, fj);
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > min )	/* add the score into list */
		    min = addnode(c, ci, cj, i, j, max_align, min);
	        }
	  }
}

//--------------------------------------------------------------------------
/* Determine the left and top boundaries of the recomputed area */
void Lalign::locate()
{ 
  int  i, j;			/* row and column indices */
  int  c;			/* best score at current point */
  int  f;			/* best score ending with insertion */
  int  d;			/* best score ending with deletion */
  int  p;			/* best score at (i-1, j-1) */
  int  ci, cj;		/* end-point associated with c */
  int  di, dj;		/* end-point associated with d */
  int  fi, fj;		/* end-point associated with f */
  int  pi, pj;		/* end-point associated with p */
  short  cflag, rflag;			/* for recomputation */
  int  *va;				/* pointer to v(A[i], B[j]) */
  int  limit;				/* the bound on j */

	/* Reverse pass
	   rows
	   CC : the scores on the current row
	   RR and EE : the endpoints that lead to CC
	   DD : the deletion scores 
	   SS and FF : the endpoints that lead to DD

	   columns
	   HH : the scores on the current columns
	   II and JJ : the endpoints that lead to HH
	   WW : the deletion scores
	   XX and YY : the endpoints that lead to WW
	*/
	for ( j = nn; j >= n1 ; j-- )
          {  CC[j] = 0;
	     EE[j] = j;
	     DD[j] = - (q);
	     FF[j] = j;
	     if ( nseq == 2 || j > mm )
                RR[j] = SS[j] = mm + 1;
	     else
                RR[j] = SS[j] = j;
	  }

        for ( i = mm; i >= m1; i-- )
	  {  c = p = 0;
	     f = - (q);
	     ci = fi = i;
	     pi = i + 1;
	     cj = fj = pj = nn + 1;
	     va = v[(unsigned)A[i]];
	     if ( nseq == 2 || n1 > i )
		limit = n1;
	     else
		limit = i + 1;
	     for ( j = nn; j >= limit ; j-- )  
	       {  f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj);
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj);
		  c = 0;
		  DIAG(i, j, c, p+va[(unsigned)B[j]]);	/* diagonal */
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj);
		  ORDER(c, ci, cj, f, fi, fj);
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > min )
		    flag = true;
	        }
	     if ( nseq == 2 || i < n1 )
	       { HH[i] = CC[n1];
	         II[i] = RR[n1];
	         JJ[i] = EE[n1];
	         WW[i] = f;
	         XX[i] = fi;
	         YY[i] = fj;
	       }
	  }
      
  for ( rl = m1, cl = n1; ; )
    { for ( rflag = cflag = 1; ( rflag && m1 > 1 ) || ( cflag && n1 > 1 ) ;  )
        { if ( rflag && m1 > 1 )	/* Compute one row */
            { rflag = 0;
	      m1--;
      	      c = p = 0;
	      f = - (q);
	      ci = fi = m1;
	      pi = m1 + 1;
	      cj = fj = pj = nn + 1;
	      va = v[(unsigned)A[m1]];
	      for ( j = nn; j >= n1 ; j-- )  
	        { f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj);
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj);
		  c = 0;
		  DIAG(m1, j, c, p+va[(unsigned)B[j]]);		/* diagonal */
		  if ( c <= 0 )
		    { c = 0; ci = m1; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj);
		  ORDER(c, ci, cj, f, fi, fj);
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > min )
		     flag = true;
		  if ( ! rflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
	 		                            || (fi > rl && fj > cl )) )
		      rflag = 1;
	        }
	      HH[m1] = CC[n1];
	      II[m1] = RR[n1];
	      JJ[m1] = EE[n1];
	      WW[m1] = f;
	      XX[m1] = fi;
	      YY[m1] = fj;
	      if ( ! cflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
			     || (fi > rl && fj > cl )) )
	         cflag = 1;
	    }

	  if ( nseq == 1 && n1 == (m1 + 1) && ! rflag )
	     cflag = 0;
	  if ( cflag && n1 > 1 )	/* Compute one column */
	    { cflag = 0;
	      n1--;
	      c = 0;
	      f = - (q);
	      cj = fj = n1;
	      va = v[(unsigned)B[n1]];
	      if ( nseq == 2 || mm < n1 )
		{ p = 0;
	          ci = fi = pi = mm + 1;
	          pj = n1 + 1;
		  limit = mm;
		}
	      else
		{ p = HH[n1];
		  pi = II[n1];
		  pj = JJ[n1];
	          ci = fi = n1;
		  limit = n1 - 1;
		}
	      for ( i = limit; i >= m1 ; i-- )  
	        { f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj);
		  c = HH[i] - qr; 
		  ci = II[i];
		  cj = JJ[i];
		  d = WW[i] - r;
		  di = XX[i];
		  dj = YY[i];
		  ORDER(d, di, dj, c, ci, cj);
		  c = 0;
	          DIAG(i, n1, c, p+va[(unsigned)A[i]]);
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = n1; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj);
		  ORDER(c, ci, cj, f, fi, fj);
		  p = HH[i];
		  HH[i] = c;
		  pi = II[i];
		  pj = JJ[i];
		  II[i] = ci;
		  JJ[i] = cj;
		  WW[i] = d;
		  XX[i] = di;
		  YY[i] = dj;
		  if ( c > min )
		     flag = true;
	          if ( ! cflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
		               || (fi > rl && fj > cl) ) )
		     cflag = 1;
	        }
	      CC[n1] = HH[m1];
	      RR[n1] = II[m1];
	      EE[n1] = JJ[m1];
	      DD[n1] = f;
	      SS[n1] = fi;
	      FF[n1] = fj;
	      if ( ! rflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
		                                 || (fi > rl && fj > cl) ) )
	         rflag = 1;
	    }
	}
      if ( (m1 == 1 && n1 == 1) || no_cross() )
	 break;
   }
  m1--;
  n1--;
}

//--------------------------------------------------------------------------
/* recompute the area on forward pass */
void Lalign::small_pass(int count)
{ int  i, j;			/* row and column indices */
  int  c;			/* best score at current point */
  int  f;			/* best score ending with insertion */
  int  d;			/* best score ending with deletion */
  int  p;			/* best score at (i-1, j-1) */
  int  ci, cj;		/* end-point associated with c */
  int  di, dj;		/* end-point associated with d */
  int  fi, fj;		/* end-point associated with f */
  int  pi, pj;		/* end-point associated with p */
  int  *va;				/* pointer to v(A[i], B[j]) */
  int  limit;				/* lower bound on j */

	for ( j = n1 + 1; j <= nn ; j++ )
	  {  CC[j] = 0;
	     RR[j] = m1;
	     EE[j] = j;
	     DD[j] = - (q);
	     SS[j] = m1;
	     FF[j] = j;
	  }
	for ( i = m1 + 1; i <= mm; i++) 
	  {  c = 0;				/* Initialize column 0 */
	     f = - (q);
	     ci = fi = i;
	     va = v[(unsigned)A[i]];
	     if ( nseq == 2 || i <= n1 )
	       { p = 0;
	         pi = i - 1;
	         cj = fj = pj = n1;
		 limit = n1 + 1;
	       }
	     else
	       { p = CC[i];
		 pi = RR[i];
		 pj = EE[i];
	         cj = fj = i;
		 limit = i + 1;
	       }
	     for ( j = limit ; j <= nn ; j++ )  
	       {  f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj);
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj);
		  c = 0;
		  DIAG(i, j, c, p+va[(unsigned)B[j]]);		/* diagonal */
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj);
		  ORDER(c, ci, cj, f, fi, fj);
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > min )	/* add the score into list */
		    min = addnode(c, ci, cj, i, j, count, min);
	        }
	  }
}

//--------------------------------------------------------------------------
/* Add a new node into list.  */
int Lalign::addnode(int c, int ci, int cj, int i, int j, int K, int cost)
{ 
  short found;				/* 1 if the node is in LIST */
  int d;

  found = 0;
  if ( most != 0 && most->STARI == ci && most->STARJ == cj )
    found = 1;
  else
     for ( d = 0; d < numnode ; d++ )
	{ most = LIST[d];
	  if ( most->STARI == ci && most->STARJ == cj )
	    { found = 1;
	      break;
	    }
        }
  if ( found )
    { if ( most->SCORE < c )
        { most->SCORE = c;
          most->ENDI = i;
          most->ENDJ = j;
        }
      if ( most->TOP > i ) most->TOP = i;
      if ( most->BOT < i ) most->BOT = i;
      if ( most->LEFT > j ) most->LEFT = j;
      if ( most->RIGHT < j ) most->RIGHT = j;
    }
  else
    { if ( numnode == K )	/* list full */
	 most = low;
      else
         most = LIST[numnode++];
      most->SCORE = c;
      most->STARI = ci;
      most->STARJ = cj;
      most->ENDI = i;
      most->ENDJ = j;
      most->TOP = most->BOT = i;
      most->LEFT = most->RIGHT = j;
    }
  if ( numnode == K )
    { if ( low == most || ! low ) 
        { for ( low = LIST[0], d = 1; d < numnode ; d++ )
            if ( LIST[d]->SCORE < low->SCORE )
              low = LIST[d];
	}
      return ( low->SCORE ) ;
    }
  else
    return cost;
}

//--------------------------------------------------------------------------
/* Find and remove_self_hits the largest score in list */
Lalign::vertexptr Lalign::findmax()
{ 
  vertexptr  cur;
  int i, j;

  for ( j = 0, i = 1; i < numnode ; i++ )
    if ( LIST[i]->SCORE > LIST[j]->SCORE )
       j = i;
  cur = LIST[j];
  if ( j != --numnode )
    { LIST[j] = LIST[numnode];
      LIST[numnode] =  cur;
    }
  most = LIST[0];
  if ( low == cur ) low = LIST[0];
  return ( cur );
}

//--------------------------------------------------------------------------
/* return 1 if no node in LIST share vertices with the area */
int Lalign::no_cross()
{ 
  vertexptr  cur;
  int i;

      for ( i = 0; i < numnode; i++ )
	{ cur = LIST[i];
	  if ( cur->STARI <= mm && cur->STARJ <= nn && cur->BOT >= m1-1 && 
	       cur->RIGHT >= n1-1 && ( cur->STARI < rl || cur->STARJ < cl ))
	     { if ( cur->STARI < rl ) rl = cur->STARI;
	       if ( cur->STARJ < cl ) cl = cur->STARJ;
	       flag = true;
	       break;
	     }
	}
      if ( i == numnode )
	return 1;
      else
	return 0;
}

//--------------------------------------------------------------------------
/* diff(tb,te) returns the score of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */

int Lalign::diff(char *A, char *B, int M, int N, int tb, int te)
{ 
  int   midi, midj, type;	/* Midpoint, type, and cost */
  int midc;

{ int   i, j;
  int c, e, d, s;
           int t, *va;

/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) DEL(M);
      return - gap(M);
    }
  if (M <= 1)
    { if (M <= 0)
        { INS(N);
          return - gap(N);
        }
      if (tb > te) tb = te;
      midc = - (tb + r + gap(N) );
      midj = 0;
      va = v[(unsigned)A[1]];
      for (j = 1; j <= N; j++)
        {  for ( tt = 1, z = row[I+1]; z != 0; z = z->NEXT )	
              if ( z->COL == j+J )			
	         { tt = 0; break; }		
           if ( tt )			
            { c = va[(unsigned)B[j]] - ( gap(j-1) + gap(N-j) );
              if (c > midc)
               { midc = c;
                 midj = j;
               }
	    }
	}
      if (midj == 0)
        { INS(N); DEL(1); }
      else
        { if (midj > 1) INS(midj-1);
	REP();
	  if ( A[1] == B[midj] )
	     no_mat += 1;
	  else
	     no_mis += 1;
	  /* mark (A[I],B[J]) as used: put J into list row[I] */	
          I++; J++;
	  //z = ( pairptr ) ckalloc( (int) sizeof(pair));
	  z = new pair[(int) sizeof(pair)];
          z->COL = J;			
          z->NEXT = row[I];				
	  row[I] = z;
          if (midj < N) INS(N-midj);
        }
      return midc;
    }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;			/* Forward phase:                          */
  CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
  t = -q;
  for (j = 1; j <= N; j++)
    { CC[j] = t = t-r;
      DD[j] = t-q;
    }
  t = -tb;
  for (i = 1; i <= midi; i++)
    { s = CC[0];
      CC[0] = c = t = t-r;
      e = t-q;
      va = v[(unsigned)A[i]];
      for (j = 1; j <= N; j++)
        { if ((c = c - qr) > (e = e - r)) e = c;
          if ((c = CC[j] - qr) > (d = DD[j] - r)) d = c;
	  DIAG(i+I, j+J, c, s+va[(unsigned)B[j]]);
          if (c < d) c = d;
          if (c < e) c = e;
          s = CC[j];
          CC[j] = c;
          DD[j] = d;
        }
    }
  DD[0] = CC[0];

  RR[N] = 0;			/* Reverse phase:                          */
  t = -q;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
  for (j = N-1; j >= 0; j--)
    { RR[j] = t = t-r;
      SS[j] = t-q;
    }
  t = -te;
  for (i = M-1; i >= midi; i--)
    { s = RR[N];
      RR[N] = c = t = t-r;
      e = t-q;
      va = v[(unsigned)A[i+1]];
      for (j = N-1; j >= 0; j--)
        { if ((c = c - qr) > (e = e - r)) e = c;
          if ((c = RR[j] - qr) > (d = SS[j] - r)) d = c;
	  DIAG(i+1+I, j+1+J, c, s+va[(unsigned)B[j+1]]);
          if (c < d) c = d;
          if (c < e) c = e;
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
    if ((c = DD[j] + SS[j] + q) > midc)
      { midc = c;
        midj = j;
        type = 2;
      }
}

/* Conquer: recursively around midpoint */

  if (type == 1)
    { (void) diff(A,B,midi,midj,tb,q);
      (void) diff(A+midi,B+midj,M-midi,N-midj,q,te);
    }
  else
    { (void) diff(A,B,midi-1,midj,tb,0);
      DEL(2);
      (void) diff(A+midi+1,B+midj,M-midi-1,N-midj,0,te);
    }
  return midc;
}
//--------------------------------------------------------------------------
// double Lalign::z_score(int nbrep)
// {
//   CountStat stat;
//   char *tmp=A;
  
//   int score_ori=getScore();
//   int nb_sub=count;

//   A=new char[M+2];

//   for(int rep=0;rep<nbrep;rep++)
//     {
//       for(int idx=1;idx<M;idx++)
// 	A[idx]=tmp[alea.rndLoHi(1,M)];
//       align();
//       for(int al=1;al<nb_sub;al++)
// 	findNext();
//       stat.add(getScore());
//     }
//   delete [] A;

//   A=tmp;

//   return (score_ori-stat.mean())/stat.sd();
// };
//--------------------------------------------------------------------------





