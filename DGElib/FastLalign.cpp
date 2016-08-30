#include <limits.h>
#include "FastLalign.h"

//--------------------------------------------------------------------------
int FastLalign::score_pass_no_long_gap(void)
{
  int gapAscores,vtCellBestScore,diagCellBestScore;  
  int i,j;

  int *hzCellBestScores=new int[N+1];
  int* gapBscores=new int[N+1];

  score=endi=endj=0;
  for(j=1;j<=N;j++)
    {
      hzCellBestScores[j]=0;
      gapBscores[j]=-q;
    }
  for(i=1;i<=M;i++)
    {
      vtCellBestScore=diagCellBestScore=0;
      gapAscores=-q;
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
	  vtCellBestScore=search_max(0,diagCellBestScore+v[(unsigned)A[i]][(unsigned)B[j]],
			      gapBscores[j],gapAscores);

	  diagCellBestScore=hzCellBestScores[j];
	  hzCellBestScores[j]=vtCellBestScore;
	  if(vtCellBestScore>score)
	    {
	      score=vtCellBestScore;
	      endi=i;
	      endj=j;
	    }
	}
    }
  delete [] gapBscores;
  delete [] hzCellBestScores;

  return score;
};
//--------------------------------------------------------------------------
void FastLalign::align_pass_no_long_gap(void)
{
  int gapAscores,vtCellBestScore,diagCellBestScore;  
  int i,j;
  int *hzCellBestScores=new int[N+1];
  int* gapBscores=new int[N+1];

  int posJumpA=0;
  int* posJumpB=new int[N+1];

  score=endi=endj=0;
  path[0]=0;
  for(j=1;j<=N;j++)
    {
      hzCellBestScores[j]=0;
      gapBscores[j]=-q;
      path[0+(M+1)*j]=0;
    }
  for(i=1;i<=M;i++)
    {
      vtCellBestScore=diagCellBestScore=0;
      gapAscores=-q;
      path[i+(M+1)*0]=0;

      for(j=1;j<=N;j++)
	{
	  //score with a gap in seq A
	  if(gapAscores>=vtCellBestScore-q)
	    {
	      gapAscores-=r;
	    }
	  else
	    {
	      gapAscores=vtCellBestScore-q-r;
	      posJumpA=j-1;
	    }
	  
	  //score with a gap in seq B
	  if(gapBscores[j]>=hzCellBestScores[j]-q)
	    {
	      gapBscores[j]-=r;
	    }
	  else
	    {
	      gapBscores[j]=hzCellBestScores[j]-q-r;
	      posJumpB[j]=i-1;
	    }

	  //max at (i,j)
	  vtCellBestScore=search_max(0,diagCellBestScore+v[(unsigned)A[i]][(unsigned)B[j]],
			      gapBscores[j],gapAscores);

	  switch(whosmax)
	    {
	    case 1:
	      { 
		path[i+(M+1)*j]=0; 
		break;
	      } //stop
	    case 2: //no gap
	      { 
		path[i+(M+1)*j]=1;
		break;
	      } 
	    case 3: //gap in seq B 
	      {
		  path[i+(M+1)*j]=-((i-posJumpB[j])+1); 
		  break;
	      }
	    case 4: //gap in seq A
	      {
		path[i+(M+1)*j]=(j-posJumpA)+1;
		break;
	      }
	    }
	  
	  diagCellBestScore=hzCellBestScores[j];
	  hzCellBestScores[j]=vtCellBestScore;
	  if(vtCellBestScore>score)
	    {
	      score=vtCellBestScore;
	      endi=i;
	      endj=j;
	    }
	}
    }
  delete [] gapBscores;
  delete [] hzCellBestScores;
  delete [] posJumpB;
};
//--------------------------------------------------------------------------
int FastLalign::score_pass_long_gap(void)
{
  int gapAscores,vtCellBestScore,diagCellBestScore;  
  int i,j;
  int *hzCellBestScores=new int[N+1];
  int* gapBscores=new int[N+1];

  int scoreLongGapA;

  int* scoreLongGapB=new int[N+1];
  std::deque<int> kPrevIscores;
  std::deque<int>* kPrevJscores=new std::deque<int>[N+1];

  score=endi=endj=0;
  for(j=1;j<=N;j++)
    {
      hzCellBestScores[j]=0;
      gapBscores[j]=-q;
      scoreLongGapB[j]=0;
    }
  for(i=1;i<=M;i++)
    {
      vtCellBestScore=diagCellBestScore=0;
      gapAscores=-q;
      scoreLongGapA=0;
      kPrevIscores.clear();
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
	  if(j>gaplen)
	    { 
	      int tmpscore=kPrevIscores.front()-long_gap_pen;
	      if(tmpscore>scoreLongGapA)
		{
		  scoreLongGapA=tmpscore;
		}
	    }
	  int gapA=(gapAscores>scoreLongGapA)? gapAscores : scoreLongGapA;

	  //score with a gap in seq B
	  if(gapBscores[j]>=hzCellBestScores[j]-q)
	    {
	      gapBscores[j]=gapBscores[j]-r;
	    }
	  else
	    {
	      gapBscores[j]=hzCellBestScores[j]-q-r;
	    }
	  if(i>gaplen)
	    {
	      int tmpscore=kPrevJscores[j].front()-long_gap_pen;
	      if(tmpscore>scoreLongGapB[j])
		{
		  scoreLongGapB[j]=tmpscore;
		}
	    }

	  int gapB=(gapBscores[j]>scoreLongGapB[j])? 
	    gapBscores[j] : scoreLongGapB[j];

	  //max at (i,j)
	  vtCellBestScore=search_max(0,diagCellBestScore+v[(unsigned)A[i]][(unsigned)B[j]],gapB,gapA);

	  kPrevJscores[j].push_back(vtCellBestScore);
	  if(i>gaplen) kPrevJscores[j].pop_front();

	  kPrevIscores.push_back(vtCellBestScore);
	  if(j>gaplen) kPrevIscores.pop_front();

	  diagCellBestScore=hzCellBestScores[j];
	  hzCellBestScores[j]=vtCellBestScore;
	  if(vtCellBestScore>score)
	    {
	      score=vtCellBestScore;
	      endi=i;
	      endj=j;
	    }
	}
    }
  delete [] scoreLongGapB;
  delete [] kPrevJscores;
  delete [] gapBscores;
  delete [] hzCellBestScores;
  return score;
};
//--------------------------------------------------------------------------
void FastLalign::align_pass_long_gap(void)
{
  int gapAscores,vtCellBestScore,diagCellBestScore;  
  int i,j;
  int *hzCellBestScores=new int[N+1];
  int* gapBscores=new int[N+1];


  int posJumpA=0,scoreLongGapA=0,posLongGapA=0;

  int* posJumpB=new int[N+1];
  int* scoreLongGapB=new int[N+1];
  int* posLongGapB=new int[N+1];
  std::deque<int> kPrevIscores;
  std::deque<int>* kPrevJscores=new std::deque<int>[N+1];

  score=endi=endj=0;
  path[0]=0;
  for(j=1;j<=N;j++)
    {
      hzCellBestScores[j]=0;
      gapBscores[j]=-q;
      path[0+(M+1)*j]=0;
      scoreLongGapB[j]=0;
    }
  for(i=1;i<=M;i++)
    {
      vtCellBestScore=diagCellBestScore=0;
      gapAscores=-q;
      path[i+(M+1)*0]=0;
      scoreLongGapA=0;
      kPrevIscores.clear();
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
	      posJumpA=j-1;
	    }
	  if(j>gaplen)
	    { 
	      int tmpscore=kPrevIscores.front()-long_gap_pen;
	      if(tmpscore>scoreLongGapA)
		{
		  scoreLongGapA=tmpscore;
		  posLongGapA=j-gaplen;
		}
	    }
	  int gapA;
	  if(gapAscores>scoreLongGapA)
	    {
	      gapA=gapAscores;
	    }
	  else
	    {
	      gapA=scoreLongGapA;
	      posJumpA=posLongGapA;
	    }

	  //score with a gap in seq B
	  if(gapBscores[j]>=hzCellBestScores[j]-q)
	    {
	      gapBscores[j]=gapBscores[j]-r;
	    }
	  else
	    {
	      gapBscores[j]=hzCellBestScores[j]-q-r;
	      posJumpB[j]=i-1;
	    }
	  if(i>gaplen)
	    {
	      int tmpscore=kPrevJscores[j].front()-long_gap_pen;
	      if(tmpscore>scoreLongGapB[j])
		{
		  scoreLongGapB[j]=tmpscore;
		  posLongGapB[j]=i-gaplen;
		}
	    }

	  int gapB;
	  if(gapBscores[j]>scoreLongGapB[j])
	    {
	      gapB=gapBscores[j];
	    }
	  else
	    {
	      gapB=scoreLongGapB[j];
	      posJumpB[j]=posLongGapB[j];
	    }

	  //max at (i,j)
	  vtCellBestScore=search_max(0,diagCellBestScore+v[(unsigned)A[i]][(unsigned)B[j]],gapB,gapA);

	  kPrevJscores[j].push_back(vtCellBestScore);
	  if(i>gaplen) kPrevJscores[j].pop_front();

	  kPrevIscores.push_back(vtCellBestScore);
	  if(j>gaplen) kPrevIscores.pop_front();

	  switch(whosmax)
	    {
	    case 1:
	      { 
		path[i+(M+1)*j]=0; 
		break;
	      } //stop
	    case 2: //no gap
	      { 
		path[i+(M+1)*j]=1;
		break;
	      } 
	    case 3: //gap in seq B 
	      {
		  path[i+(M+1)*j]=-((i-posJumpB[j])+1); 
		  break;
	      }
	    case 4: //gap in seq A
	      {
		path[i+(M+1)*j]=(j-posJumpA)+1;
		break;
	      }
	    }
	  
	  diagCellBestScore=hzCellBestScores[j];
	  hzCellBestScores[j]=vtCellBestScore;
	  if(vtCellBestScore>score)
	    {
	      score=vtCellBestScore;
	      endi=i;
	      endj=j;
	    }
	}
    }
  delete [] gapBscores;
  delete [] hzCellBestScores;
  delete [] posJumpB;
  delete [] scoreLongGapB;
  delete [] posLongGapB;
  delete [] kPrevJscores;
};
//--------------------------------------------------------------------------
void FastLalign::traceback(void)
{
  starti=endi;
  startj=endj;
  alignment.clear();
  if(score==0)
    {
      startj=endj=starti=endi=0;;
      return;
    }
  while(starti>=0 && startj>=0)
      {
	if(path[starti+(M+1)*startj]==0)
	  {
	    int d=alignment.back();
	    if(d==0)
	      {
		starti++;
		startj++;
	      }
	    else
		std::cout<<"backtrack error!!"<<std::endl;
	    return;
	  }
	if(path[starti+(M+1)*startj]==1) 
	  {
	    alignment.push_back(0);
	    starti--; startj--;
	  }
	else
	  {
	    if(path[starti+(M+1)*startj]>1) // gap  in seq A
	      {
		int size_gap=path[starti+(M+1)*startj]-1;
		alignment.push_back(size_gap);
		startj-=size_gap;
	      }
	    if(path[starti+(M+1)*startj]<-1) // gap  in seq B
	      {
		int size_gap=-path[starti+(M+1)*startj]-1;
		alignment.push_back(-size_gap);
		starti-=size_gap;
	      }
	  }
      }
};
//--------------------------------------------------------------------------
// double FastLalign::z_score(int nbrep)
// {
//   CountStat stat;
//   char *tmp=A;
  
//   int score_ori=getScore();
//   A=new char[M+2];

//   for(int rep=0;rep<nbrep;rep++)
//     {
//       for(int idx=1;idx<=M;idx++)
// 	A[idx]=tmp[alea.rndLoHi(1,M)];
//       if(gaplen==0) stat.add(score_pass_no_long_gap());
//       else stat.add(score_pass_long_gap());
//     }
//   delete [] A;

//   A=tmp;

//   double zval=0;
//   if(stat.sd()!=0)
//       zval=(score_ori-stat.mean())/stat.sd();
//   else
//     zval=(score_ori-stat.mean())/(1e-100);
//   return zval;
// };
//--------------------------------------------------------------------------
double FastLalign::getIdentity(void)
{
  if(identity!=0) return identity;
  if(alignment.empty()) traceback();
  if(score==0)
    {
      identity=0;
      return identity;
    }

  int   i=starti,  j=startj;
  unsigned short count_match=0,count=0;
 
  while (i <= endi || j <= endj)
    { 
      int dir=alignment.back();
      if(dir==0)
	{
	  if(A[i++]==B[j++]) count_match++;
	  count++;
	}
      else
	if(dir>0)
	  {
	    j++;
	  }
	else
	  {
	    i++;
	  }
      if(!alignment.empty()) alignment.pop_back();
    }
  identity=(double)count_match/count;
  return identity;
}
//--------------------------------------------------------------------------
void FastLalign::view()
{
  /* Alignment display routine */
  // Show the current best alignment.
  if(alignment.empty()) traceback();

  if(score==0)
    {
      std::cout<<" no alignment!"<<std::endl;
      return;
    }

  char aline[51],bline[51],cline[51];
  int   i=starti,  j=startj;
  unsigned short count=0;

  unsigned short count_match=0,count_diag=0;

  int recomp_score=0;

  int as=starti;
  int bs=startj;
  while (i <= endi || j <= endj)
    { 
      int dir=alignment.back();
      if(dir==0)
	{
	  aline[count]=A[i];
	  bline[count]=B[j];
	  recomp_score+=v[(unsigned)A[i]][(unsigned)B[j]];
	  if(A[i] ==B[j]) count_match++;
	  count_diag++;
	  cline[count]=(A[i++] ==B[j++] ) ? '|' : ' ';
	  count++;
	  if(count==50)
	    {
	      aline[count]='\0';
	      bline[count]='\0';
	      cline[count]='\0';
	      count=0;
	      std::cout<<std::setw(7)<<"   "<<" ";
	      for(int p=1;p<=50;p++)
		{
		  if(p%10==0) std::cout<<":";
		  else if(p%5==0) std::cout<<"+";
		  else std::cout<<".";
		}
	      std::cout<<std::endl;
	      std::cout<<std::setw(7)<<as<<" "<<aline<<" "<<i-1<<std::endl;
	      std::cout<<std::setw(7)<<"   "<<" "<<cline<<std::endl;
	      std::cout<<std::setw(7)<<bs<<" "<<bline<<" "<<j-1<<std::endl;
	      std::cout<<std::endl;
	      as=i;
	      bs=j;
	    };
	}
      else
	if(dir>0)
	  {
	    int size_gap=dir;
	    if(size_gap<=gaplen)
	      recomp_score-=(q+(size_gap)*r);
	    else
	       recomp_score-=long_gap_pen;
	    for(int p=0;p<size_gap;p++)
	    {
	      aline[count]='-';
	      bline[count]=B[j++];
	      cline[count]=' ';
	      count++;
	      if(count==50)
		{
		  aline[count]='\0';
		  bline[count]='\0';
		  cline[count]='\0';
		  count=0;
		  std::cout<<std::setw(7)<<"   "<<" ";
		  for(int p=1;p<=50;p++)
		    {
		      if(p%10==0) std::cout<<":";
		      else if(p%5==0) std::cout<<"+";
		      else std::cout<<".";
		    }
		  std::cout<<std::endl;
		  std::cout<<std::setw(7)<<as<<" "<<aline<<" "<<i-1<<std::endl;
		  std::cout<<std::setw(7)<<"   "<<" "<<cline<<std::endl;
		  std::cout<<std::setw(7)<<bs<<" "<<bline<<" "<<j-1<<std::endl;
		  std::cout<<std::endl;
		  as=i;
		  bs=j;
		};
	    }		  
	  }
	else
	  {
	    int size_gap=-dir;
	    if(size_gap<=gaplen)
	      recomp_score-=(q+(size_gap)*r);
	    else
	       recomp_score-=long_gap_pen;
	    for(int p=0;p<size_gap;p++)
	    {
	      aline[count]=A[i++];
	      bline[count]='-';
	      cline[count]=' ';
	      count++;
	      if(count==50)
		{
		  aline[count]='\0';
		  bline[count]='\0';
		  cline[count]='\0';
		  count=0;
		  std::cout<<std::setw(7)<<"   "<<" ";
		  for(int p=1;p<=50;p++)
		    {
		      if(p%10==0) std::cout<<":";
		      else if(p%5==0) std::cout<<"+";
		      else std::cout<<".";
		    }
		  std::cout<<std::endl;
		  std::cout<<std::setw(7)<<as<<" "<<aline<<" "<<i-1<<std::endl;
		  std::cout<<std::setw(7)<<"   "<<" "<<cline<<std::endl;
		  std::cout<<std::setw(7)<<bs<<" "<<bline<<" "<<j-1<<std::endl;
		  std::cout<<std::endl;
		  as=i;
		  bs=j;
		};
	    }	
	  }
      if(!alignment.empty())alignment.pop_back();
    }
  if(count!=0)
    {
      aline[count]='\0';
      bline[count]='\0';
      cline[count]='\0';
      std::cout<<std::setw(7)<<"   "<<" ";
      for(int p=1;p<=count;p++)
	{
	  if(p%10==0) std::cout<<":";
	  else if(p%5==0) std::cout<<"+";
	  else std::cout<<".";
	}
      std::cout<<std::endl;
      std::cout<<std::setw(7)<<as<<" "<<aline<<" "<<i-1<<std::endl;
      std::cout<<std::setw(7)<<"   "<<" "<<cline<<std::endl;
      std::cout<<std::setw(7)<<bs<<" "<<bline<<" "<<j-1<<std::endl;
      std::cout<<std::endl;
    }

  identity=(double)count_match/count_diag;
  std::cout<<"match="<<match<<" mismatch="<<mismh
      <<" gap_open="<<q<<" gap_extend="<<r
      <<" max gap length penalized="<<gaplen<<std::endl;
  std::cout<<"      Similarity Score : "<<recomp_score<<std::endl;
  std::cout<<"   Identity percentage : "<<identity<<std::endl;
  std::cout<<"      Begins at ("<<starti<<","<<startj<<") and Ends at ("
      <<endi<<","<<endj<<")\n"<<std::flush;
}
//--------------------------------------------------------------------------
void FastLalign::getAlignedSeq(SDGString& seq1, SDGString& seq2)
{
  if(alignment.empty()) traceback();
  if(score==0)
    {
      std::cout<<" no alignment!"<<std::endl;
      return;
    }

  int   i=starti,  j=startj;
  unsigned short count_match=0,count_diag=0;

  while (i <= endi || j <= endj)
    { 
      int dir=alignment.back();
      if(dir==0)
	{
	  seq1+=A[i];
	  seq2+=B[j];
	  if(A[i++] ==B[j++]) count_match++;
	  count_diag++;
	}
      else
	if(dir>0)
	  {
	    int size_gap=dir;
	    for(int p=0;p<size_gap;p++)
	    {
	      seq1+='-';
	      seq2+=B[j++];
	    }		  
	  }
	else
	  {
	    int size_gap=-dir;
	    for(int p=0;p<size_gap;p++)
	    {
	      seq1+=A[i++];
	      seq2+='-';
	    }	
	  }
      if(!alignment.empty())alignment.pop_back();
    }

  identity=(double)count_match/count_diag;
}

























