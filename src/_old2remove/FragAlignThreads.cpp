#include "FragAlignThreads.h"

//------------------------------------------------------------------------
void FragAlignThreads::join(std::list<RangePair>& l, std::list<RangePairSet>& finished_list)
{


	//start joining
	 std::list<std::list<RangePair>::iterator>
		 path_dd, path_cd, path_dc, path_cc;

	 finished_list.clear(); //to be returned
	 while(!l.empty())
	  {
		alignDirectDirect(l,path_dd);
		if(!path_dd.empty())
		  {
		RangePairSet newRange=*(path_dd.front());
		std::list<RangePair> a_path;
		for(std::list<std::list<RangePair>::iterator>::iterator r=path_dd.begin();
			r!=path_dd.end();r++)
		  {
			a_path.push_back(*(*r));
			l.erase(*r);
		  }
		newRange.setPath(a_path);
		newRange.computeScoreWithDynaProg(mism_pen,gapo_pen,gape_pen);
		finished_list.push_back(newRange);
		  }
		if(l.empty()) break;

		alignDirectCompl(l,path_dc);
		if(!path_dc.empty())
		  {
		RangePairSet newRange=*(path_dc.front());

		std::list<RangePair> a_path;
		for(std::list<std::list<RangePair>::iterator>::iterator r=path_dc.begin();
			r!=path_dc.end();r++)
		  {
			a_path.push_back(*(*r));
			l.erase(*r);
		  }
		newRange.setPath(a_path);
		newRange.computeScoreWithDynaProg(mism_pen,gapo_pen,gape_pen);
		finished_list.push_back(newRange);
		  }

		if(l.empty()) break;
		alignComplDirect(l,path_cd);
		if(!path_cd.empty())
		  {
		RangePairSet newRange=*(path_cd.front());

		std::list<RangePair> a_path;
		for(std::list<std::list<RangePair>::iterator>::iterator r=path_cd.begin();
			r!=path_cd.end();r++)
		  {
			a_path.push_back(*(*r));
			l.erase(*r);
		  }
		newRange.setPath(a_path);
		newRange.computeScoreWithDynaProg(mism_pen,gapo_pen,gape_pen);
		finished_list.push_back(newRange);
		  }

		if(l.empty()) break;
		alignComplCompl(l,path_cc);
		if(!path_cc.empty())
		  {
		RangePairSet newRange=*(path_cc.front());

		std::list<RangePair> a_path;
		for(std::list<std::list<RangePair>::iterator>::iterator r=path_cc.begin();
			r!=path_cc.end();r++)
		  {
			a_path.push_back(*(*r));
			l.erase(*r);
		  }
		newRange.setPath(a_path);
		newRange.computeScoreWithDynaProg(mism_pen,gapo_pen,gape_pen);
		finished_list.push_back(newRange);
		  }
	  } //end while


}
//------------------------------------------------------------------------
void FragAlignThreads::alignDirectDirect(std::list<RangePair>& l, std::list<std::list<RangePair>::iterator>& path_result)
{
	I.clear();
	V.clear();
	unsigned nb_frag=0;
	V.push_back(0);
	unsigned k=0;
	path_result.clear();

  for(std::list<RangePair>::iterator i=l.begin();i!=l.end();i++)
    {
      if(!i->getRangeQ().isPlusStrand() || !i->getRangeS().isPlusStrand() ) 
	continue;
      nb_frag++;
      bound b;
      unsigned o=over;
      if(i->getRangeQ().getLength()/2<over 
	 || i->getRangeS().getLength()/2<over)
	o=std::min(i->getRangeQ().getLength()/2,i->getRangeS().getLength()/2);
      b.end_point=i->getRangeQ().getStart()+o;
      b.y_high=i->getRangeS().getStart()+o;
      b.y_low=i->getRangeS().getEnd()-o;
      b.left=true;
      b.iter_range_pair=i;  
      b.rect=++k;  
      I.push_back(b);
      
      b.end_point=i->getRangeQ().getEnd()-o;
      b.y_high=i->getRangeS().getStart()+o;
      b.y_low=i->getRangeS().getEnd()-o;
      b.left=false;
      b.iter_range_pair=i;  
      b.rect=k;  
      I.push_back(b);
      
      V.push_back(0);
    }

  if(nb_frag>0)
  {
	  align(gapo_pen,gape_pen,mism_pen,connect_dist_limit,I,V,nb_frag,path_result);
  }
}
//------------------------------------------------------------------------
void FragAlignThreads::alignDirectCompl(std::list<RangePair>& l, std::list<std::list<RangePair>::iterator>& path_result)
{
	I.clear();
	V.clear();
	unsigned nb_frag=0;
	V.push_back(0);
	unsigned k=0;
	path_result.clear();

  unsigned long max_coord=0;
  for(std::list<RangePair>::iterator i=l.begin();i!=l.end();i++)
    {
      unsigned long m=std::max(i->getRangeS().getStart(),i->getRangeS().getEnd());
      max_coord=m>max_coord?m:max_coord;
    }
  max_coord++;
  for(std::list<RangePair>::iterator i=l.begin();i!=l.end();i++)
    {
      if(!i->getRangeQ().isPlusStrand() || i->getRangeS().isPlusStrand())
	continue;
      nb_frag++;
      bound b;
      unsigned o=over;
      if(i->getRangeQ().getLength()/2<over 
	 || i->getRangeS().getLength()/2<over)
	o=std::min(i->getRangeQ().getLength()/2,i->getRangeS().getLength()/2);
      b.end_point=i->getRangeQ().getStart()+o;
      b.y_high=max_coord-i->getRangeS().getStart()+o;
      b.y_low=max_coord-i->getRangeS().getEnd()-o;
      b.left=true;
      b.iter_range_pair=i;  
      b.rect=++k;  
      I.push_back(b);
      
      b.end_point=i->getRangeQ().getEnd()-o;
      b.y_high=max_coord-i->getRangeS().getStart()+o;
      b.y_low=max_coord-i->getRangeS().getEnd()-o;
      b.left=false;
      b.iter_range_pair=i;  
      b.rect=k;  
      I.push_back(b);
      
      V.push_back(0);
    }

  if(nb_frag>0)
  {
	  align(gapo_pen,gape_pen,mism_pen,connect_dist_limit,I,V,nb_frag,path_result);
  }

}
//------------------------------------------------------------------------
void FragAlignThreads::alignComplDirect(std::list<RangePair>& l, std::list<std::list<RangePair>::iterator>& path_result)
{
	I.clear();
	V.clear();
	unsigned nb_frag=0;
	V.push_back(0);
	unsigned k=0;
	path_result.clear();

	unsigned long max_coord=0;
  for(std::list<RangePair>::iterator i=l.begin();i!=l.end();i++)
    {
      unsigned long m=std::max(i->getRangeQ().getStart(),i->getRangeQ().getEnd());
      max_coord=m>max_coord?m:max_coord;
    }
  max_coord++;
  for(std::list<RangePair>::iterator i=l.begin();i!=l.end();i++)
    {
      if(i->getRangeQ().isPlusStrand() || !i->getRangeS().isPlusStrand())
	continue;
      nb_frag++;
      bound b;
      unsigned o=over;
      if(i->getRangeQ().getLength()/2<over 
	 || i->getRangeS().getLength()/2<over)
	o=std::min(i->getRangeQ().getLength()/2,i->getRangeS().getLength()/2);
      b.end_point=max_coord-i->getRangeQ().getStart()+o;
      b.y_high=i->getRangeS().getStart()+o;
      b.y_low=i->getRangeS().getEnd()-o;
      b.left=true;
      b.iter_range_pair=i;  
      b.rect=++k;  
      I.push_back(b);
      
      b.end_point=max_coord-i->getRangeQ().getEnd()-o;
      b.y_high=i->getRangeS().getStart()+o;
      b.y_low=i->getRangeS().getEnd()-o;
      b.left=false;
      b.iter_range_pair=i;  
      b.rect=k;  
      I.push_back(b);
      
      V.push_back(0);
    }

  if(nb_frag>0)
  {
	  align(gapo_pen,gape_pen,mism_pen,connect_dist_limit,I,V,nb_frag,path_result);
  }
}
//------------------------------------------------------------------------
void  FragAlignThreads::alignComplCompl(std::list<RangePair>& l, std::list<std::list<RangePair>::iterator>& path_result)
{
	I.clear();
	V.clear();
	unsigned nb_frag=0;
	V.push_back(0);
	unsigned k=0;
	path_result.clear();

	unsigned long max_coordQ=0;
	  for(std::list<RangePair>::iterator i=l.begin();i!=l.end();i++)
		{
		  unsigned long m=std::max(i->getRangeQ().getStart(),i->getRangeQ().getEnd());
		  max_coordQ=m>max_coordQ?m:max_coordQ;
		}
	  max_coordQ++;
	  unsigned long max_coordS=0;
	  for(std::list<RangePair>::iterator i=l.begin();i!=l.end();i++)
		{
		  unsigned long m=std::max(i->getRangeQ().getStart(),i->getRangeQ().getEnd());
		  max_coordS=m>max_coordS?m:max_coordS;
		}
	  max_coordS++;
	  for(std::list<RangePair>::iterator i=l.begin();i!=l.end();i++)
		{
		  if(i->getRangeQ().isPlusStrand() || i->getRangeS().isPlusStrand())
		continue;
		  nb_frag++;
		  bound b;
		  unsigned o=over;
		  if(i->getRangeQ().getLength()/2<over
		 || i->getRangeS().getLength()/2<over)
		o=std::min(i->getRangeQ().getLength()/2,i->getRangeS().getLength()/2);
		  b.end_point=max_coordQ-i->getRangeQ().getStart()+o;
		  b.y_high=max_coordS-i->getRangeS().getStart()+o;
		  b.y_low=max_coordS-i->getRangeS().getEnd()-o;
		  b.left=true;
		  b.iter_range_pair=i;
		  b.rect=++k;
		  I.push_back(b);

		  b.end_point=max_coordQ-i->getRangeQ().getEnd()-o;
		  b.y_high=max_coordS-i->getRangeS().getStart()+o;
		  b.y_low=max_coordS-i->getRangeS().getEnd()-o;
		  b.left=false;
		  b.iter_range_pair=i;
		  b.rect=k;
		  I.push_back(b);

		  V.push_back(0);
		}
	  if(nb_frag>0)
	  {
		  align(gapo_pen,gape_pen,mism_pen,connect_dist_limit,I,V,nb_frag,path_result);
	  }
}
//------------------------------------------------------------------------
void FragAlignThreads::align(double gapo_pen,double gape_pen, double mism_pen, unsigned connect_dist_limit, std::list<bound>& I,
		  std::vector<long long>& V, unsigned nb_frag,
		  std::list<std::list<RangePair>::iterator> &path_result)
{
		//Gusfield, page 328

  std::vector< std::list<std::list<RangePair>::iterator> > path((nb_frag)+1);
  unsigned end_best_path=0;
  std::list<triple> L; //list of vertical rectangle coordinates (kept sorted)
  V[0]=0;
  I.sort(bound::less);
  std::list<bound>::iterator i=I.begin();
  while(i!=I.end()) // for each hz rect coordinate
    {
      unsigned k=i->rect;
      if(i->left) //left coordinate: compute the score and chain 
	{
	  unsigned long hk=i->y_high;

	  long long max_score=std::numeric_limits<long long>::min();
	  std::list<triple>::iterator j=L.begin();
	  std::list<triple>::iterator max_j=L.end();
	  while(j!=L.end() && j->lower<hk) // modif Hadi: search for ALL rectangle above the current
	    {
	      long long diag1=
		(long long)j->iter_range_pair->getRangeQ().getEnd()
		-(long long)j->iter_range_pair->getRangeS().getEnd();
	      long long diag2=
		(long long)i->iter_range_pair->getRangeQ().getStart()
		-(long long)i->iter_range_pair->getRangeS().getStart();
	      long long gap=llabs(diag1-diag2);
		  
	      long long mismatch; 
	      if(diag1>diag2)
		mismatch=
		  (long long)i->iter_range_pair->getRangeQ().getStart()
		  -(long long)j->iter_range_pair->getRangeQ().getEnd();
	      else
		mismatch=
		  (long long)i->iter_range_pair->getRangeS().getStart()
		  -(long long)j->iter_range_pair->getRangeS().getEnd();
	      
	      mismatch=mismatch>0?mismatch:0;
	      
	      long long score=i->iter_range_pair->getScore()
		+V[j->rect]
		-(long long)rint(gapo_pen+gape_pen*gap)
		-(long long)rint(mismatch*mism_pen);
	      
	      if(score>max_score)
		{
		  max_score=score;
		  max_j=j;
		}
	      j++;
	    }


	  if(max_j==L.end())
	    {
		  V[k]=i->iter_range_pair->getScore();
	      path[k].push_back(i->iter_range_pair);
	    }
	  else
	    {
	      if(max_score>=i->iter_range_pair->getScore())
	      {
	    	  V[k]=max_score;
	    	  path[k].insert(path[k].end(),path[max_j->rect].begin(),
	    			  path[max_j->rect].end());
	    	  path[k].push_back(i->iter_range_pair);
	      }
	      else
	      {
	    	  V[k]=i->iter_range_pair->getScore();
	    	  path[k].push_back(i->iter_range_pair);
	      }
	    }

	  if(V[k]>V[end_best_path]) end_best_path=k;
	}
      else // right coordinate: add vertical coord as a triple
	{
	  unsigned long lk=i->y_low;
	  std::list<triple>::iterator j=std::lower_bound(L.begin(),L.end(),lk,
							 triple::lessI);
	  if(j!=L.end()) j--;
	  long long Vj=0;
	  unsigned long lj=0;
	  if(j!=L.end())
	    {
	      Vj=V[j->rect];
	      lj=j->lower;
	    }
	  if(lj<lk || (lj==lk && V[k]>Vj))
	    {
	      triple t;
	      t.lower=i->y_low;
	      t.path_val=V[k];
	      t.iter_range_pair=i->iter_range_pair;
	      t.rect=k;
	      L.insert(std::lower_bound(L.begin(),L.end(),t,triple::less),t);
	    }
//	  std::list<triple>::iterator n=L.begin();
//  	  while(n!=L.end())
//  	    if(V[k]>V[n->rect] && n->lower<=lk) 
//  		n=L.erase(n);
//  	    else n++;
	  std::list<triple>::iterator n=L.begin();
 	  while(n!=L.end())
 	    if(V[k]>V[n->rect] && n->lower<=lk
	       && (i->iter_range_pair->getRangeQ().getStart()
		   -n->iter_range_pair->getRangeQ().getEnd()
		   > connect_dist_limit
		   ||
		   i->iter_range_pair->getRangeS().getStart()
		   -n->iter_range_pair->getRangeS().getEnd()
		   > connect_dist_limit
		   ) 
	       // modif Hadi to take into account above modifications 
	       //and set limits to search for above rectangles
	       )
 		n=L.erase(n);
 	    else n++;

	}
      i++;
    }
  path_result=path[end_best_path];
};
////------------------------------------------------------------------------
//long long FragAlignThreads::align(std::list<bound> *pI, std::vector<long long> *pV, std::list<std::list<RangePair>::iterator> *pPath, unsigned nb_frag )
//{
//  //Gusfield, page 328
//
//  static pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER;
//
//  std::vector< std::list<std::list<RangePair>::iterator> > path(nb_frag+1);
//  unsigned end_best_path=0;
//  std::list<triple> L; //list of vertical rectangle coordinates (kept sorted)
//  pthread_mutex_lock (&my_mutex);
//  pV->operator[](0)=0;
//  pI->sort(bound::less);
//  pthread_mutex_unlock (&my_mutex);
//  std::list<bound>::iterator i=pI->begin();
//  while(i!=pI->end()) // for each hz rect coordinate
//    {
//      unsigned k=i->rect;
//      if(i->left) //left coordinate: compute the score and chain
//	{
//	  unsigned long hk=i->y_high;
//
//	  long long max_score=std::numeric_limits<long long>::min();
//	  std::list<triple>::iterator j=L.begin();
//	  std::list<triple>::iterator max_j=L.end();
//	  while(j!=L.end() && j->lower<hk) // modif Hadi: search for ALL rectangle above the current
//	    {
//	      long long diag1=
//		(long long)j->iter_range_pair->getRangeQ().getEnd()
//		-(long long)j->iter_range_pair->getRangeS().getEnd();
//	      long long diag2=
//		(long long)i->iter_range_pair->getRangeQ().getStart()
//		-(long long)i->iter_range_pair->getRangeS().getStart();
//	      long long gap=llabs(diag1-diag2);
//
//	      long long mismatch;
//	      if(diag1>diag2)
//		mismatch=
//		  (long long)i->iter_range_pair->getRangeQ().getStart()
//		  -(long long)j->iter_range_pair->getRangeQ().getEnd();
//	      else
//		mismatch=
//		  (long long)i->iter_range_pair->getRangeS().getStart()
//		  -(long long)j->iter_range_pair->getRangeS().getEnd();
//
//	      mismatch=mismatch>0?mismatch:0;
//
//	      long long score=i->iter_range_pair->getScore()
//		+pV->operator[](j->rect)
//		-(long long)rint(gapo_pen+gape_pen*gap)
//		-(long long)rint(mismatch*mism_pen);
//
//	      if(score>max_score)
//		{
//		  max_score=score;
//		  max_j=j;
//		}
//	      j++;
//	    }
//
//
//	  if(max_j==L.end())
//	    {
//		  pthread_mutex_lock (&my_mutex);
//	      pV->operator[](k)=i->iter_range_pair->getScore();
//	      pthread_mutex_unlock (&my_mutex);
//	      path[k].push_back(i->iter_range_pair);
//	    }
//	  else
//	    {
//	      if(max_score>=i->iter_range_pair->getScore())
//	      {
//	    	  pthread_mutex_lock (&my_mutex);
//	    	  pV->operator[](k)=max_score;
//	    	  pthread_mutex_unlock (&my_mutex);
//	    	  path[k].insert(path[k].end(),path[max_j->rect].begin(),
//	    			  path[max_j->rect].end());
//	    	  path[k].push_back(i->iter_range_pair);
//	      }
//	      else
//	      {
//	    	  pthread_mutex_lock (&my_mutex);
//	    	  pV->operator[](k)=i->iter_range_pair->getScore();
//	    	  pthread_mutex_unlock (&my_mutex);
//	    	  path[k].push_back(i->iter_range_pair);
//	      }
//	    }
//
//	  if(pV->operator[](k)>pV->operator[](end_best_path)) end_best_path=k;
//	}
//      else // right coordinate: add vertical coord as a triple
//	{
//	  unsigned long lk=i->y_low;
//	  std::list<triple>::iterator j=std::lower_bound(L.begin(),L.end(),lk,
//							 triple::lessI);
//	  if(j!=L.end()) j--;
//	  long long Vj=0;
//	  unsigned long lj=0;
//	  if(j!=L.end())
//	    {
//	      Vj=pV->operator[](j->rect);
//	      lj=j->lower;
//	    }
//	  if(lj<lk || (lj==lk && pV->operator[](k)>Vj))
//	    {
//	      triple t;
//	      t.lower=i->y_low;
//	      t.path_val=pV->operator[](k);
//	      t.iter_range_pair=i->iter_range_pair;
//	      t.rect=k;
//	      L.insert(std::lower_bound(L.begin(),L.end(),t,triple::less),t);
//	    }
////	  std::list<triple>::iterator n=L.begin();
////  	  while(n!=L.end())
////  	    if(V[k]>V[n->rect] && n->lower<=lk)
////  		n=L.erase(n);
////  	    else n++;
//	  std::list<triple>::iterator n=L.begin();
// 	  while(n!=L.end())
// 	    if(pV->operator[](k)>pV->operator[](n->rect) && n->lower<=lk
//	       && (i->iter_range_pair->getRangeQ().getStart()
//		   -n->iter_range_pair->getRangeQ().getEnd()
//		   > connect_dist_limit
//		   ||
//		   i->iter_range_pair->getRangeS().getStart()
//		   -n->iter_range_pair->getRangeS().getEnd()
//		   > connect_dist_limit
//		   )
//	       // modif Hadi to take into account above modifications
//	       //and set limits to search for above rectangles
//	       )
// 		n=L.erase(n);
// 	    else n++;
//
//	}
//      i++;
//    }
//  pthread_mutex_lock (&my_mutex);
//  *pPath=path[end_best_path];
//  pthread_mutex_unlock (&my_mutex);
//  return pV->operator[](end_best_path);
//};









