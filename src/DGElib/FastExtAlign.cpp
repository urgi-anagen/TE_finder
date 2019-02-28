#include <limits.h>
#include "FastExtAlign.h"

//--------------------------------------------------------------------------
int FastExtAlign::extend_dir(int start_score)
{
  int gapAscores,vtCellBestScore,diagCellBestScore;  
  unsigned i,j;

  if(starti+M>=lenA)
    M=lenA-starti;
  if(startj+N>=lenB)
    N=lenB-startj;

  int *hzCellBestScores=new int[N+1];
  int* gapBscores=new int[N+1];

  score=start_score;
  endi=starti;
  endj=startj;
  int t=start_score+r;
  for(j=1;j<=N;j++)
    {
      t=t-r;
      hzCellBestScores[j]=t;
      gapBscores[j]=t-q;
    }
  t=start_score+r;
  for(i=1;i<=M;i++)
    {
      t=t-r;
      vtCellBestScore=diagCellBestScore=t;
      gapAscores=t-q;
      for(j=1;j<=N;j++)
	{
	  //score with a gap in seq A
	  if(gapAscores>=vtCellBestScore-q)
	    {
	      gapAscores=gapAscores-r;
	    }
	  else
	    {
	      gapAscores=vtCellBestScore-q-r;
	    }

	  //score with a gap in seq B
	  if(gapBscores[j]>=hzCellBestScores[j]-q)
	    {
	      gapBscores[j]=gapBscores[j]-r;
	    }
	  else
	    {
	      gapBscores[j]=hzCellBestScores[j]-q-r;
	    }
	  //max at (i,j)
	  vtCellBestScore=search_max(diagCellBestScore+v[(unsigned)A[starti+i]][(unsigned)B[startj+j]],
			      gapBscores[j],gapAscores);

	  diagCellBestScore=hzCellBestScores[j];
	  hzCellBestScores[j]=vtCellBestScore;
	  if(vtCellBestScore>score)
	    {
	      score=vtCellBestScore;
	      endi=starti+i;
	      endj=startj+j;
	    }
	}
    }
  delete [] gapBscores;
  delete [] hzCellBestScores;

  return score;
};
//--------------------------------------------------------------------------
int FastExtAlign::extend_rev(int start_score)
{
  int gapAscores,vtCellBestScore,diagCellBestScore;  
  unsigned i,j;

  if(starti<=M)
    M=starti-1;
  if(startj<=N)
    N=startj-1;

  int *hzCellBestScores=new int[N+1];
  int* gapBscores=new int[N+1];

  score=start_score;
  endi=starti;
  endj=startj;
  int t=start_score+r;
  for(j=1;j<=N;j++)
    {
      t=t-r;
      hzCellBestScores[j]=t;
      gapBscores[j]=t-q;
    }
  t=start_score+r;
  for(i=1;i<=M;i++)
    {
      t=t-r;
      vtCellBestScore=diagCellBestScore=t;
      gapAscores=t-q;
      for(j=1;j<=N;j++)
	{
	  //score with a gap in seq A
	  if(gapAscores>=vtCellBestScore-q)
	    {
	      gapAscores=gapAscores-r;
	    }
	  else
	    {
	      gapAscores=vtCellBestScore-q-r;
	    }

	  //score with a gap in seq B
	  if(gapBscores[j]>=hzCellBestScores[j]-q)
	    {
	      gapBscores[j]=gapBscores[j]-r;
	    }
	  else
	    {
	      gapBscores[j]=hzCellBestScores[j]-q-r;
	    }
	  //max at (i,j)
	  vtCellBestScore=search_max(diagCellBestScore+v[(unsigned)A[starti-i]][(unsigned)B[startj-j]],
			      gapBscores[j],gapAscores);

	  diagCellBestScore=hzCellBestScores[j];
	  hzCellBestScores[j]=vtCellBestScore;
	  if(vtCellBestScore>score)
	    {
	      score=vtCellBestScore;
	      endi=starti-i;
	      endj=startj-j;
	    }
	}
    }
  delete [] gapBscores;
  delete [] hzCellBestScores;

  return score;
};


























