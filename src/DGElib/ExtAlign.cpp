#include <limits.h>
#include "ExtAlign.h"

//--------------------------------------------------------------------------
int ExtAlign::score_pass(void)
{
  int gapAscores,vtCellBestScore,diagCellBestScore;  
  int i,j;

  int *hzCellBestScores=new int[N+1];
  int* gapBscores=new int[N+1];

  score=endi=endj=0;
  int t=q;
  for(j=1;j<=N;j++)
    {
      t=t+r;
      hzCellBestScores[j]=-(t);
      gapBscores[j]=-(t+q);
    }
  t=q;
  for(i=1;i<=M;i++)
    {
      t=t+r;
      vtCellBestScore=diagCellBestScore=-(t);
      gapAscores=-(t+q);
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
	  vtCellBestScore=search_max(diagCellBestScore+v[(unsigned)A[i]][(unsigned)B[j]],
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
void ExtAlign::align_pass(void)
{
  int gapAscores,vtCellBestScore,diagCellBestScore;  
  int i,j;
  int *hzCellBestScores=new int[N+1];
  int* gapBscores=new int[N+1];

  int posJumpA=0;
  int* posJumpB=new int[N+1];

  score=endi=endj=0;
  path[0]=0;
  int t=q;
  for(j=1;j<=N;j++)
    {
      t=t+r;
      hzCellBestScores[j]=-(t);
      gapBscores[j]=-(t+q);
      path[0+(M+1)*j]=(j+1); 
    }
  t=q;
  for(i=1;i<=M ;i++)
    {
      t=t+r;
      vtCellBestScore=diagCellBestScore=-(t);
      gapAscores=-(t+q);
      path[i+(M+1)*0]=-(i+1);
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

	  vtCellBestScore=search_max(diagCellBestScore+v[(unsigned)A[i]][(unsigned)B[j]],
				     gapBscores[j],gapAscores);
	      switch(whosmax)
		{
		case 1: //no gap
		  { 
		    path[i+(M+1)*j]=1;
		    break;
		  } 
		case 2: //gap in seq B 
		  {
		    path[i+(M+1)*j]=-((i-posJumpB[j])+1); 
		    break;
		  }
		case 3: //gap in seq A
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
void ExtAlign::traceback(void)
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
	    starti++;
	    startj++;
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
// double ExtAlign::z_score(int nbrep)
// {
//   CountStat stat;
//   char *tmp=A;
  
//   int score_ori=getScore();
//   A=new char[M+2];

//   for(int rep=0;rep<nbrep;rep++)
//     {
//       for(int idx=1;idx<=M;idx++)
// 	A[idx]=tmp[alea.rndLoHi(1,M)];
//       stat.add(score_pass());
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
double ExtAlign::getIdentity(void)
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
 
  while (i <= endi && j <= endj)
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
void ExtAlign::view()
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
	    recomp_score-=(q+(size_gap)*r);
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
	    recomp_score-=(q+(size_gap)*r);
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
      <<" gap_open="<<q<<" gap_extend="<<r<<std::endl;
  std::cout<<"      Similarity Score : "<<recomp_score<<std::endl;
  std::cout<<"   Identity percentage : "<<identity<<std::endl;
  std::cout<<"      Begins at ("<<starti<<","<<startj<<") and Ends at ("
      <<endi<<","<<endj<<")\n"<<std::flush;
}


























