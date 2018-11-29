
/**
 *
 * BLRMatchMap.cpp
 *
 **/
#include <stdlib.h>
#include <regex.h>
#include <fstream>
#include "FragAlign.h"
#include "BLRMatchMap.h"
#include <SDGError.h>
#include "BLRMatchMapLoader.h"
	// TODO CLEAN CODE (methods for test and comments ...)
	//----------------------------------------------------------------------------
	void BLRMatchMap::insert(RangePair& rangePair)
	{
	  //insert rangePair in the right place
	  std::list<RangePair>& al_list
	    =map_align[Key(rangePair.getRangeQ().getNumChr(),
			   rangePair.getRangeS().getNumChr())];

	  std::list<RangePair>::iterator r=std::lower_bound(al_list.begin(),
							    al_list.end(),rangePair);
	  if( rangePair.getRangeQ().getMin() != r->getRangeQ().getMin()
	      || rangePair.getRangeQ().getMax() != r->getRangeQ().getMax()
	      || rangePair.getRangeS().getMin() != r->getRangeS().getMin()
	      || rangePair.getRangeS().getMax() != r->getRangeS().getMax()
	      || rangePair.getE_value() != r->getE_value()
	      || rangePair.getScore() != r->getScore()
	      || rangePair.getIdentity() != r->getIdentity() )
	    al_list.insert(r,rangePair);
	}
	//----------------------------------------------------------------------------
	void BLRMatchMap::readAlign(std::istringstream& streamName, int verbose)
	//  read from a istringstream
	{
	   BLRMatchMapLoader blrmm = BLRMatchMapLoader(); 
	   blrmm.readAlign(*this, streamName, verbose);
	}
	//----------------------------------------------------------------------------
	// TODO rename load in loadAlgin
	void BLRMatchMap::load(SDGString filename, int verbose)
	//  load from txt file
	{
	   BLRMatchMapLoader blrmm = BLRMatchMapLoader(); 
	   blrmm.loadAlign(*this, filename, verbose);
	}
	//----------------------------------------------------------------------------
	void BLRMatchMap::readPath(std::istringstream& streamName, int verbose)
	{
	   BLRMatchMapLoader blrmm = BLRMatchMapLoader(); 
	   blrmm.readPath(*this, streamName, verbose);
	}
	//---------------------------------------------------------------------------
	// TODO verbosity > 0 change getNbMatchesInMapAlign
	void BLRMatchMap::loadPath(SDGString filename, int verbose)
	{
	   BLRMatchMapLoader blrmm = BLRMatchMapLoader();
	   blrmm.loadPath(*this, filename, verbose);
	}
	//----------------------------------------------------------------------------
	void BLRMatchMap::add_clean(std::list<RangePair>& rp_list,
			      std::list<RangePair>::iterator iter)
	//add a range pair and post process it removing conflicting subjects
	{
	  // search for a conflicting subject
	  bool found_over=false;
	  std::list<RangePair> lrp; //list of modified and cleaned RangePair
	  lrp.push_back(*iter);
	  for(MapAlign::iterator m=map_align.begin(); m!=map_align.end() ;m++)
	    if(m->first.first==iter->getRangeQ().getNumChr() &&
	       m->first.second!=iter->getRangeS().getNumChr())
	      {
		// check overlap only with a different subject

			  for(std::list<RangePair>::iterator lrp_it=lrp.begin();
				  lrp_it!=lrp.end();
				  lrp_it++)
				  {
						for(std::list<RangePair>::iterator iter_list=m->second.begin();
						iter_list!=m->second.end();
						iter_list++)
						  {
							if(lrp_it->getScore()<iter_list->getScore()
							   &&  lrp_it->overlapQ(*iter_list))
							  {
								found_over=true;
								RangePair rp=lrp_it->diffQ(*iter_list);
								if(!rp.empty()
								   && rp.getRangeQ().getLength()>para->getLenFilter())
								{
									lrp.push_back(rp);
								}
							  } //end if (...)
						} //end loop for
				  }//end loop for
	      } //end if

	  if(found_over) // RangePair found to overlap (conflicts!)
	    {
	      for(std::list<RangePair>::iterator lrp_it=lrp.begin();
		  lrp_it!=lrp.end();
		  lrp_it++)

			if(!lrp_it->empty()
			   && lrp_it->getRangeQ().getLength()>para->getLenFilter())
			  {
				std::list<RangePair>::iterator it
				  =std::lower_bound(iter,rp_list.end(),
						*lrp_it,
						RangePair::greaterScore); // search for the right place to insert
				while(it!=rp_list.end() && it==iter)
				  it++;
				rp_list.insert(it,*lrp_it);
		  }
	    }
	  else // already cleaned RangePair
	    if (!iter->empty() && iter->getRangeQ().getLength()>para->getLenFilter() && iter->getScore()>0)
		insert(*iter);
	}
	//----------------------------------------------------------------------------
	void BLRMatchMap::clean_conflicts(void)
	// removing conflicting subjects
	{
	  std::list<RangePair> rp_list;
	  for(MapAlign::iterator m=map_align.begin(); m!=map_align.end();m++)
	    {
	      while(!m->second.empty())
			{
			  RangePair rp=m->second.back();
			  m->second.pop_back();
			  if(para->getEvalFilter()>=rp.getE_value()
				 || para->getIdFilter()<=rp.getIdentity()
				 || para->getLenFilter()<=rp.getLength())
				{
				  rp_list.push_back(rp);
				}
			}
	    }
	  map_align.clear();

	  rp_list.sort(RangePair::greaterScore);
	  for(std::list<RangePair>::iterator i=rp_list.begin();
	      i!=rp_list.end();i++)
	      add_clean(rp_list,i);
	}
	//----------------------------------------------------------------------------
	void BLRMatchMap::insert_path(RangePairSet& rangePair)
	//insert a RangePairSet in map_path at the right place
	{
	    std::list<RangePairSet>& al_list
	    =map_path[Key(rangePair.getRangeQ().getNumChr(),
			  rangePair.getRangeS().getNumChr())];

	  std::list<RangePairSet>::iterator r=std::lower_bound(al_list.begin(),
							       al_list.end(),rangePair);
	  if( rangePair.getRangeQ().getMin() != r->getRangeQ().getMin()
	      || rangePair.getRangeQ().getMax() != r->getRangeQ().getMax()
	      || rangePair.getRangeS().getMin() != r->getRangeS().getMin()
	      || rangePair.getRangeS().getMax() != r->getRangeS().getMax()
	      || rangePair.getE_value() != r->getE_value()
	      || rangePair.getScore() != r->getScore()
	      || rangePair.getIdentity() != r->getIdentity() )
	    al_list.insert(r,rangePair);
	}
	//----------------------------------------------------------------------------
	void BLRMatchMap::add_clean_path_same_S(std::list<RangePairSet>& rp_list,
			      std::list<RangePairSet>::iterator iter)
	//add a range pair set and post process it removing conflicting subjects
	{
	  // search for a conflicting subject
	  bool found_over=false;
	  std::list<RangePairSet>& list
	    =map_path[Key(iter->getRangeQ().getNumChr(),
			  iter->getRangeS().getNumChr())];
	  for(std::list<RangePairSet>::iterator iter_list=list.begin();
	      iter_list!=list.end() ;iter_list++)
	    if(RangePair::greaterScore(*iter_list,*iter)
	       && iter->overlapQ(*iter_list))
	      if(iter->diffQ(*iter_list))
		found_over=true;

	  if(found_over)
	    {
	      if(!iter->empty()
		 && iter->getRangeQ().getLength()>para->getLenFilter())
		{
		  std::list<RangePairSet>::iterator it
		    =std::lower_bound(iter,rp_list.end(),
				      *iter,
				      RangePair::greaterScore);
		  if(it==iter)
		    it++;
		  rp_list.insert(it,*iter);
		}
	    }
	  else
	    if(iter->getRangeQ().getLength()>para->getLenFilter()
	       && iter->getScore()>0)
	      insert_path(*iter);

	}
	//----------------------------------------------------------------------------
	void BLRMatchMap::add_clean_path_all_S(std::list<RangePairSet>& rp_list,
					       std::list<RangePairSet>::iterator iter,
					       int verbose)
	//add a range pair set and post process it removing conflicting subjects
	{
	  bool found_over=false;
	  bool atLeastOneOverlap=false;
	  unsigned nbseqS=getNbSseq();
	  //unsigned s=1;
	  // start subject at -1 to consider merged data
	  long s=-1;
	  // for current query iter on each subject in map_path.
	  while(s<=nbseqS)
	    {
	      // list of subject      
	      std::list<RangePairSet>& list
		=map_path[Key(iter->getRangeQ().getNumChr(),s)];
	      for(std::list<RangePairSet>::iterator iter_list=list.begin();
		  iter_list!=list.end() ;iter_list++){
		    if( RangePair::greaterScore(*iter_list,*iter)
		    && iter->overlapQ(*iter_list) )
				if(iter->diffQ(*iter_list))
					found_over=true;
		}
	      if(found_over)
		{
		  atLeastOneOverlap=true;
		  if(!iter->empty()
		     && iter->getRangeQ().getLength()>para->getLenFilter())
		    {
		
		    std::list<RangePairSet>::iterator it
			=std::lower_bound(iter,rp_list.end(),
					*iter,
					RangePair::greaterScore);
		      if(it==iter)
			it++;
		      rp_list.insert(it,*iter);
	       

		    }
		}

	      if(!found_over) s++;else found_over=false;
	    }
	  if(!atLeastOneOverlap && iter->getRangeQ().getLength()>para->getLenFilter()
	     && iter->getScore()>0){
	    	insert_path(*iter);
	  }
	}
	/*
	//----------------------------------------------------------------------------
	void BLRMatchMap::add_clean_path_all_S(std::list<RangePairSet>& rp_list,
					       std::list<RangePairSet>::iterator iter,
					       int verbose)
	//add a range pair set and post process it removing conflicting subjects
	{
	  bool found_over=false;
	  bool atLeastOneOverlap=false;

	  unsigned nbseqS=getNbSseq();
	  unsigned s=1;

	  while(s<=nbseqS)
	    {
	      // for each query/subject in map_path. query is:iter->getRangeQ().getNumChr()  subjet: all subjects s 1->nbseqS
	      std::list<RangePairSet>& list
		=map_path[Key(iter->getRangeQ().getNumChr(),s)];
	      if(verbose>0)
		std::cout<<"qry  "<<numQ2name(iter->getRangeQ().getNumChr())<<" sbj "<<numS2name(s)<<": list.size="<<list.size()<<std::endl<<std::flush;
	      for(std::list<RangePairSet>::iterator iter_list=list.begin();
		  iter_list!=list.end() ;iter_list++)
		if( RangePair::greaterScore(*iter_list,*iter)
		    && iter->overlapQ(*iter_list) )
		  if(iter->diffQ(*iter_list))
		    found_over=true;

	      if(found_over)
		{
		  if(verbose>0)
		    std::cout<<"found overlap"<<std::endl;
		  atLeastOneOverlap=true;
		  if(!iter->empty()
		     && iter->getRangeQ().getLength()>para->getLenFilter())
		    {
		      std::list<RangePairSet>::iterator it
			=std::lower_bound(iter,rp_list.end(),
					*iter,
					RangePair::greaterScore);
		      if(it==iter)
			it++;
		      rp_list.insert(it,*iter);
		    }
		}
	      else
		if(verbose>0)
		  std::cout<<"no overlap"<<std::endl;

	      if(!found_over) s++;else found_over=false;
	    }

	  if(!atLeastOneOverlap && iter->getRangeQ().getLength()>para->getLenFilter()
	     && iter->getScore()>0)
	    insert_path(*iter);
	}
	/
	*/
	//----------------------------------------------------------------------------
	void BLRMatchMap::clean_path(bool same_S, int verbose)
		// removing conflicting subjects
	{
		std::list<RangePairSet> rp_list;
		for(MapPath::iterator m=map_path.begin(); m!=map_path.end();m++)
		{
			while(!m->second.empty())
			{
			  RangePairSet rp=m->second.back();
			  m->second.pop_back();
			  rp.computeScoreWithLength();
			  //rp.computeScoreWithDynaProg();
			  rp_list.push_back(rp);
			}
		}
	  	map_path.clear();
	  	rp_list.sort(RangePair::greaterScore);

	  if(same_S)
	    for(std::list<RangePairSet>::iterator i=rp_list.begin();
		i!=rp_list.end();i++)
	      add_clean_path_same_S(rp_list,i);
	  else
	    for(std::list<RangePairSet>::iterator i=rp_list.begin();
		i!=rp_list.end();i++)
		  {
			if(verbose>0)
				std::cout<<"Add "<<numQ2name(i->getRangeQ().getNumChr())<<"-"<<numS2name(i->getRangeS().getNumChr())
			<<"-"<<i->getScore()<<std::endl<<std::flush;
			add_clean_path_all_S( rp_list, i, verbose-1 );
			if(verbose>0)
			{
			  std::cout<<"nb of matches: "<<getNbMatchesInMapPath()<<std::endl;
			  std::cout<<"nb of paths: "<<getNbDistinctPaths()<<std::endl;
			}
		  }
	}
	//----------------------------------------------------------------------------
	void BLRMatchMap::add_split_path(std::list<RangePairSet>& rp_list, std::list<RangePairSet>::iterator iter)
	//add a range pair set and post process it removing split nest
	{
	  bool found_over=false;
	  std::list<RangePairSet> lrp;
	  lrp.push_back(*iter);

	    for(MapPath::iterator m=map_path.begin(); m!=map_path.end() ;m++)
	    {
		if(m->first.first==iter->getRangeQ().getNumChr())
		{
		  for(std::list<RangePairSet>::iterator lrp_it=lrp.begin();
		      lrp_it!=lrp.end();
		      lrp_it++)
		    for(std::list<RangePairSet>::iterator iter_list=m->second.begin();
		      iter_list!=m->second.end() ;iter_list++)
		      if(lrp_it->overlapQ(*iter_list))
			{
			  if(lrp_it->getLength()>=100 &&
			     lrp_it->inserted(*iter_list) &&
			     fabs(iter_list->getIdentity()
				  -lrp_it->getIdentity())<=para->getIdTolerance())
			    continue;
			  std::list<RangePairSet> lrp2;
			  if(lrp_it->split(*iter_list,lrp2))
			    {
			    found_over=true;
			      for(std::list<RangePairSet>::iterator lrp2_it
				    =lrp2.begin();lrp2_it!=lrp2.end(); lrp2_it++)
				if(lrp2_it->getRangeQ().getLength()
				   >para->getLenFilter())
				  lrp.push_back(*lrp2_it);
			    }
			}
		} //end if
	    }
	    
	  if(found_over)
	    {
	      for(std::list<RangePairSet>::iterator lrp_it=lrp.begin();
		  lrp_it!=lrp.end();lrp_it++)
		if(!lrp_it->empty()
		   && lrp_it->getRangeQ().getLength()>para->getLenFilter())
		  {
		    std::list<RangePairSet>::iterator it
		      =std::lower_bound(iter,rp_list.end(),
					*lrp_it,
					RangePair::lessIdentity);
		    if(it==iter)
		      it++;
		    rp_list.insert(it,*lrp_it);
		  }
	    }
	  else
	    if(iter->getRangeQ().getLength()>para->getLenFilter())
	      insert_path(*iter);

	}

//----------------------------------------------------------------------------
bool BLRMatchMap::isOverlapFound_in_add_split_path(std::list<RangePairSet>::iterator iter, MapPath mapPath, double idTolerance, unsigned lenFilter)
{

    bool found_over=false;
    std::list<RangePairSet> lrp;
    
    lrp.push_back(*iter);
  for(MapPath::iterator m=mapPath.begin(); m!=mapPath.end() ;m++)
    {

      if(m->first.first==iter->getRangeQ().getNumChr())
	{

	  for(std::list<RangePairSet>::iterator lrp_it=lrp.begin();
	      lrp_it!=lrp.end();
	      lrp_it++)
	    for(std::list<RangePairSet>::iterator iter_list=m->second.begin();
	      iter_list!=m->second.end() ;iter_list++)
	      if(lrp_it->overlapQ(*iter_list))
		{
		  if(lrp_it->getLength()>=100 &&
		     lrp_it->inserted(*iter_list) &&
		     fabs(iter_list->getIdentity()
			  -lrp_it->getIdentity())<=idTolerance)
		    continue;
		  std::list<RangePairSet> lrp2;

		  if(lrp_it->split(*iter_list,lrp2))
		    {
		      found_over=true;
		      for(std::list<RangePairSet>::iterator lrp2_it
			    =lrp2.begin();lrp2_it!=lrp2.end(); lrp2_it++)
			if(lrp2_it->getRangeQ().getLength()
			   >lenFilter)
			  lrp.push_back(*lrp2_it);
		    }

		}
	}
    }
	return found_over; 
}

//----------------------------------------------------------------------------
void BLRMatchMap::split_path(void)
// split nested path
{
  std::list<RangePairSet> rp_list;
  for(MapPath::iterator m=path_begin(); m!=path_end();m++)
  {
		while(!m->second.empty())
		{
			RangePairSet rp=m->second.back();
			m->second.pop_back();
			rp_list.push_back(rp);
		}
  }
  map_path.clear();
  // score no longer valid after split of RangePairSet by method
  // split
  rp_list.sort(RangePair::lessIdentity);
  for(std::list<RangePairSet>::iterator i=rp_list.begin();
	  i!=rp_list.end();i++)
	add_split_path(rp_list, i);
}
//----------------------------------------------------------------------------
void BLRMatchMap::insert_path_static(MapPath mapPath, RangePairSet& rangePair)
//insert a RangePairSet in map_path at the right place
{
  std::list<RangePairSet>& al_list
    =mapPath[Key(rangePair.getRangeQ().getNumChr(),
		  rangePair.getRangeS().getNumChr())];

  std::list<RangePairSet>::iterator r=std::lower_bound(al_list.begin(),
						       al_list.end(),rangePair);
  if( rangePair.getRangeQ().getMin() != r->getRangeQ().getMin()
      || rangePair.getRangeQ().getMax() != r->getRangeQ().getMax()
      || rangePair.getRangeS().getMin() != r->getRangeS().getMin()
      || rangePair.getRangeS().getMax() != r->getRangeS().getMax()
      || rangePair.getE_value() != r->getE_value()
      || rangePair.getScore() != r->getScore()
      || rangePair.getIdentity() != r->getIdentity() )
    al_list.insert(r,rangePair);
}
//----------------------------------------------------------------------------
void BLRMatchMap::mapPath(bool joining, bool clean_before, bool clean_after, bool merged, int verbose)
{
  map_path.clear();

  if(joining)
    {
      if(verbose>0)
	std::cout<<"Join parameters: dist_pen="<<para->getDist_pen()
		 <<" gap_pen="<<para->getGap_pen()
		 <<" overlap="<<para->getOverlap()<<std::endl<<std::flush;
      
	
      for(MapAlign::iterator m=map_align.begin(); m!=map_align.end();m++)
		{
		  FragAlign fragAlign(para->getDist_pen(),0,para->getGap_pen(),
					  para->getOverlap());
		  map_path[m->first]=fragAlign.join(m->second);
		  m->second.clear();
		}
      if(verbose>0){
		std::cout<<"After join:\nnb of matches: "<<getNbMatchesInMapPath()<<std::endl;
		std::cout<<"nb of paths: "<<getNbDistinctPaths()<<std::endl;

		std::cout<<"Write joined matches in BED format."<<std::endl;
		SDGString filename = para->getPrefixFileName() + ".joined.bed";
		std::list<RangePairSet> copy_list = copyRpsListFromMapPath();
		SDGString color = "253,63,146";
     		writeBED(filename, copy_list, color, verbose-1);
      }
      if (merged)
      {
		if (verbose>0)
			std::cout<<"Compute score with length."<<std::endl;
		computeScoreWithLength();
	
		// prepare rpsList for merge
		unsigned long id = 1;
		rpsList.clear();
		for( MapPath::iterator m=path_begin(); m!=path_end(); m++ )
		{
			while(!m->second.empty())
			{
				  RangePairSet rp=m->second.front();
				  m->second.pop_front();
				  rp.setId(id);
				  rpsList.push_back(rp);
				  id++;
			 }
		}
		map_path.clear();
	 	// merge
		if (verbose>0)
	       		std::cout<<"Merge on query."<<std::endl;
     	merge(verbose-1);	
      	if (verbose>0){
			std::cout<<"Write merged matches in BED format."<<std::endl;
		SDGString filename = para->getPrefixFileName() + ".merged.bed";
		std::list<RangePairSet> copy_list = copyRpsListFromMapPath();
		SDGString color= "0,255,0";
		writeBED(filename, copy_list, color, verbose-1);
	      	}
      }
      if( ( clean_before | clean_after ) & ( verbose > 0 ) )
    	  std::cout<<"Clean the connections..."<<std::endl<<std::flush;
      if(clean_before)
    	  clean_path(true,verbose-1); // clean when same subject
      if(clean_after)
    	  clean_path(false,verbose-1); // clean with all subjects
      if( ( clean_before | clean_after ) & ( verbose > 0 ) )
      {
    	  std::cout<<"After clean:\nnb of matches: "<<getNbMatchesInMapPath()<<std::endl;
    	  std::cout<<"nb of paths: "<<getNbDistinctPaths()<<std::endl;
    	  std::cout<<"Connections were cleaned when necessary."<<std::endl;
   	
	  std::cout<<"Write cleaned matches in BED format."<<std::endl;
	  SDGString filename = para->getPrefixFileName() + ".cleaned.bed";
	  std::list<RangePairSet> copy_list = copyRpsListFromMapPath();
	  SDGString color = "255,0,0";
     	  writeBED(filename, copy_list, color, verbose-1);
      }
      if(clean_before || clean_after)
		{
		 if(verbose>0)
			std::cout<<"Split the connections..."<<std::endl<<std::flush;
		 split_path();
 		 if(verbose>0){
			std::cout<<"nb of matches: "<<getNbMatchesInMapPath()<<std::endl;
			std::cout<<"nb of paths: "<<getNbDistinctPaths()<<std::endl;
			std::cout<<"Connections were splitted when necessary."<<std::endl;

	  		std::cout<<"Write split matches in BED format."<<std::endl;
	  		SDGString filename = para->getPrefixFileName() + ".split.bed";
	  		std::list<RangePairSet> copy_list = copyRpsListFromMapPath();
			SDGString color = "237,127,16";
     	  		writeBED(filename, copy_list, color, verbose-1);

		  }
		 }

    }
  else // no join
    {
	  int count=0;
	  if(verbose>0)
	  			std::cout<<"No join, considering ";
      for(MapAlign::iterator m=map_align.begin(); m!=map_align.end();m++)
		{
		  std::list<RangePairSet> path;
		  for(std::list<RangePair>::iterator i=m->second.begin();
			  i!=m->second.end();i++)
			{
			  count++;
			  path.push_back(RangePairSet(*i));
			}
		  map_path[m->first]=path;
		  m->second.clear();
		}
      if(verbose>0)
      	  			std::cout<<count<<" matches"<<std::endl;
    }
  map_align.clear();
}
//----------------------------------------------------------------------------
void BLRMatchMap::mapPathJoinOnlyForTest(bool joining, bool clean_before, bool clean_after, int verbose)
{
  map_path.clear();

  if(joining)
    {
      if(verbose>0)
	std::cout<<"Join parameters: dist_pen="<<para->getDist_pen()
		 <<" gap_pen="<<para->getGap_pen()
		 <<" overlap="<<para->getOverlap()<<std::endl<<std::flush;
      for(MapAlign::iterator m=map_align.begin(); m!=map_align.end();m++)
	{
	  FragAlign fragAlign(para->getDist_pen(),0,para->getGap_pen(),
			      para->getOverlap());
	  map_path[m->first]=fragAlign.join(m->second);
	  m->second.clear();
	}
      if(verbose>0){
	std::cout<<"nb of matches: "<<getNbMatchesInMapPath()<<std::endl;
	std::cout<<"nb of paths: "<<getNbDistinctPaths()<<std::endl;}
       }
  return;
}


//----------------------------------------------------------------------------
// TODO is this methos destructive for mapPath attribute ?
std::list<RangePairSet> BLRMatchMap::getRpsListFromMapPath(void) // Build list of RangePairSet
{
	  std::list<RangePairSet> rps_list;
	  for( MapPath::iterator m=path_begin(); m!=path_end(); m++ )
	  {
		  while(!m->second.empty())
		  {
			  RangePairSet rp=m->second.front();
			  m->second.pop_front();
			  rps_list.push_back(rp);
		  }
	  }
	  return rps_list;
};

//----------------------------------------------------------------------------
std::list<RangePairSet> BLRMatchMap::copyRpsListFromMapPath(void)
{
	std::list<RangePairSet> copy_rps_list;
	for( MapPath::iterator m=path_begin(); m!=path_end(); m++ )
	{
		for (std::list<RangePairSet>::iterator it = m->second.begin(); it != m->second.end(); it++){
			copy_rps_list.push_back(*it);
		}
	}
	return copy_rps_list;
	
}
//----------------------------------------------------------------------------
//void BLRMatchMap::join(void* arg)
//{
//	  JoinArg *p=(JoinArg*)arg;
//	  pthread_mutex_t m_lock;
//	  pthread_mutex_init(&m_lock, NULL);
//
//	  FragAlignThreads fragAlign(p->para->getDist_pen(),0,p->para->getGap_pen(),
//	  	  	  			      p->para->getOverlap());
//	  for(std::list<MapAlign::iterator>::iterator i=p->l.begin();i!=p->l.end();i++)
//	  {
//		  std::cout<<"thread "<< pthread_self()<<" working on ("<<(*i)->first.first<<","<<(*i)->first.second<<")"<<std::endl;
//		  std::list<RangePairSet> lpath;
//		  fragAlign.join((*i)->second,lpath);
//		  pthread_mutex_lock(&m_lock);
//		  std::cout<<"thread "<< pthread_self()<<" finished ("<<(*i)->first.first<<","<<(*i)->first.second<<")"<<std::endl;
//		  p->p_map_path->operator []((*i)->first)=lpath;
//		  (*i)->second.clear();
//		  pthread_mutex_unlock(&m_lock);
//
//	  }
//};
//----------------------------------------------------------------------------
//void BLRMatchMap::mapPathWithThreads(bool joining, bool clean_before, bool clean_after, bool merged, int verbose)
//{
//  map_path.clear();
//
//  if(joining)
//    {
//      if(verbose>0)
//    	  std::cout<<"nb threads="<<para->getNbThread()<<std::endl;
//	  ThreadPool tp(para->getNbThread());
//      int ret = tp.initialize_threadpool();
//      if (ret == -1)
//      {
//    	 std::cerr << "Failed to initialize thread pool!" << endl;
//    	 exit(0);
//      }
//      if(verbose>0)
//    	  std::cout<<"Join parameters: dist_pen="<<para->getDist_pen()
//    	  <<" gap_pen="<<para->getGap_pen()
//    	  <<" overlap="<<para->getOverlap()<<std::endl<<std::flush;
//
//      vector<JoinArg*> vja(map_align.size());
//	  //unsigned count=0, nb_tasks=0, nb_count=map_align.size()/para->getNbThread();
//	  unsigned count=0, nb_tasks=0, nb_count=50;
//	  std::list<MapAlign::iterator> lm;
//      for(MapAlign::iterator m=map_align.begin(); m!=map_align.end();m++)
//	  {
//		  lm.push_back(m);
//		  count+=m->second.size();
//		  //count++;
//		  std::cout<<"count="<<count<<"/"<<nb_count<<std::endl;
//		  if(count>=nb_count)
//		  {
//			  JoinArg *ja=new JoinArg();
//			  vja.push_back(ja);
//			  ja->l=lm;
//			  ja->para=para;
//			  ja->p_map_path=&map_path;
//			  Task* t = new Task(&join, (void*) ja);
//			  tp.add_task(t);
//			  nb_tasks++;
//			  std::cout<<nb_tasks<<" task(s) submitted"<<std::endl;
//			  count=0;
//			  lm.clear();
//		  }
//	  }
//	  if(!lm.empty())
//	  {
//		  JoinArg *ja=new JoinArg();
//		  vja.push_back(ja);
//		  ja->l=lm;
//		  ja->para=para;
//		  ja->p_map_path=&map_path;
//		  Task* t = new Task(&join, (void*) ja);
//		  tp.add_task(t);
//		  std::cout<<"Last tasks submitted"<<std::endl;
//	  }
//
//	 //wait thread to finish
//	 sleep(2);
//	 tp.destroy_threadpool();
//	 for(vector<JoinArg*>::iterator i=vja.begin();i!=vja.end();i++)
//		 delete *i;
//
//      if(verbose>0)
//      {
//		std::cout<<"nb of matches: "<<getNbMatchesInMapPath()<<std::endl;
//		std::cout<<"nb of paths: "<<getNbDistinctPaths()<<std::endl;
//      }
//
//
//
//      if (merged)
//      {
//      	if (verbose>0)
//       		std::cout<<"Compute score with length."<<std::endl;
//      	computeScoreWithLength();
//
//      	if (verbose>0)
//       		std::cout<<"Merge on query."<<std::endl;
//     	 merge();
//      }
//      if( ( clean_before | clean_after ) & ( verbose > 0 ) )
//    	  std::cout<<"Clean the connections..."<<std::endl<<std::flush;
//      if(clean_before)
//	clean_path(true,verbose-1); // clean when same subject
//      if(clean_after)
//	clean_path(false,verbose-1); // clean with all subjects
//      if( ( clean_before | clean_after ) & ( verbose > 0 ) )
//      {
//    	  std::cout<<"nb of matches: "<<getNbMatchesInMapPath()<<std::endl;
//    	  std::cout<<"nb of paths: "<<getNbDistinctPaths()<<std::endl;
//    	  std::cout<<"Connections were cleaned when necessary."<<std::endl;
//      }
//
//      if(clean_before || clean_after)
//	{
//	  if(verbose>0)
//	    std::cout<<"Split the connections..."<<std::endl<<std::flush;
//      	  split_path();
//	  if(verbose>0){
//	    std::cout<<"nb of matches: "<<getNbMatchesInMapPath()<<std::endl;
//	    std::cout<<"nb of paths: "<<getNbDistinctPaths()<<std::endl;
//	    std::cout<<"Connections were splitted when necessary."<<std::endl;}
//	}
//    }
//  else
//    {
//      for(MapAlign::iterator m=map_align.begin(); m!=map_align.end();m++)
//	{
//	  std::list<RangePairSet> path;
//	  for(std::list<RangePair>::iterator i=m->second.begin();
//	      i!=m->second.end();i++)
//	    {
//	      path.push_back(RangePairSet(*i));
//	    }
//	  map_path[m->first]=path;
//	  m->second.clear();
//	}
//    }
//  map_align.clear();
//}
//---------------------------------------------------------------------------
void BLRMatchMap::selectQregex(SDGString regex)
{
  regex_t preg;
  regcomp(&preg, regex, REG_EXTENDED | REG_NEWLINE | REG_NOSUB);

  MapAlign::iterator iter_hash,prev;  // iterator to visit hash_align
  iter_hash=begin();
  while (iter_hash!=end())
    {
      prev=iter_hash++;
      if(regexec(&preg,SDGString(num2nameQ[prev->first.first]),1,NULL,0)!=0)
	{
	  map_align.erase(prev);
	}
    }
  regfree(&preg);
}
//---------------------------------------------------------------------------
void BLRMatchMap::select(bool subject, bool clean_before, bool clean_after)
{
  std::ostringstream outfile;
  if(subject)
    {
      RangeMap matchmap;
      MapPath::iterator iter_hash;  // iterator to visit hash_align
      std::list<RangePairSet>::iterator iter_list;  //iterator to visit align_list

      iter_hash=path_begin();
      while (iter_hash!=path_end())
	{
	  iter_list=iter_hash->second.begin();
	  while (iter_list!=iter_hash->second.end())
	    {
	      std::string chrS_name=num2nameS[(iter_list->getRangeS().getNumChr())-1];
		matchmap.add(RangeSeq("",chrS_name,0,0));

	      iter_list++;
	    }
	  iter_hash++;
	}

      if(clean_before || clean_after)
	outfile<<para->getParameterFileName().beforelast(".param")
	       <<".clean_match.subject_selected";
      else
	outfile<<para->getParameterFileName().beforelast(".param")
	       <<".match.subject_selected";
      matchmap.selectSrcSeq(outfile.str().c_str(),subject_db);
    }
  else
    {

      RangeMap matchmap;
      MapPath::iterator iter_hash;  // iterator to visit hash_align
      std::list<RangePairSet>::iterator iter_list;  //iterator to visit align_list
      iter_hash=path_begin();
      while (iter_hash!=path_end())
	{
	  iter_list=iter_hash->second.begin();
	  while (iter_list!=iter_hash->second.end())
	    {
	      std::string chrQ_name=num2nameQ[(iter_list->getRangeQ().getNumChr())-1];
	      matchmap.add(RangeSeq("",chrQ_name,0,0));

	      iter_list++;
	    }
	  iter_hash++;
	}

      if(clean_before || clean_after)
	outfile<<para->getParameterFileName().beforelast(".param")
	       <<".clean_match.query_selected";
      else
	outfile<<para->getParameterFileName().beforelast(".param")
	       <<".match.query_selected";
      matchmap.selectSrcSeq(outfile.str().c_str(),query_db);
    }
  std::cout<<"ok!"<<std::endl;
}
//---------------------------------------------------------------------------
void BLRMatchMap::writeMatch(const SDGString& filename, int verbose)
{
  if(verbose>0)
    std::cout<<"writing 'tab' file..."<<std::flush;
  std::ofstream fout(filename);

  fout<<"query.name"
      <<"\t"<<"query.start"
      <<"\t"<<"query.end"
      <<"\t"<<"query.length"
      <<"\t"<<"query.length.%"
      <<"\t"<<"match.length.%"
      <<"\t"<<"subject.name"
      <<"\t"<<"subject.start"
      <<"\t"<<"subject.end"
      <<"\t"<<"subject.length"
      <<"\t"<<"subject.length.%"
      <<"\t"<<"E.value"
      <<"\t"<<"Score"
      <<"\t"<<"Identity"
      <<"\t"<<"path"
      <<std::endl;

  unsigned path_id=0;
  for(MapPath::iterator iter_hash=path_begin();
      iter_hash!=path_end();iter_hash++)
    {

      std::string query_name=num2nameQ[iter_hash->first.first];
      std::string subject_name;
      if(same_db)
	subject_name=num2nameQ[iter_hash->first.second];
      else
	subject_name=num2nameS[iter_hash->first.second];

      std::string subseqname(subject_name,0,subject_name.find(" "));
      unsigned querylen=query_db.find(query_name).length();
      unsigned subjectlen=subject_db.find(subject_name).length();

      for(std::list<RangePairSet>::iterator iter_list
	    =iter_hash->second.begin();iter_list!=iter_hash->second.end();
	  iter_list++)
	{
	  //iter_list->view();

	  RangeAlignSet rasQ=iter_list->getRangeAlignSetQ();
	  RangeAlignSet rasS=iter_list->getRangeAlignSetS();


	  unsigned rangeQlen=rasQ.getLengthSet();
	  unsigned rangeSlen=rasS.getLengthSet();

	  fout<<query_name
	      <<"\t"<<rasQ.getStart()
	      <<"\t"<<rasQ.getEnd()
	      <<"\t"<<rangeQlen
	      <<"\t"<<(double)(rangeQlen)/querylen
	      <<"\t"<<(double)(rangeQlen)/subjectlen
	      <<"\t"<<subseqname
	      <<"\t"<<rasS.getStart()
	      <<"\t"<<rasS.getEnd()
	      <<"\t"<<rangeSlen
	      <<"\t"<<(double)(rangeSlen)/subjectlen
	      <<"\t"<<iter_list->getE_value()
	      <<"\t"<<iter_list->getScore()
	      <<"\t"<<iter_list->getIdentity()
	      <<"\t"<<++path_id
	      <<std::endl;
	}
    }
  if(verbose>0)
    std::cout<<" done"<<std::endl;
}
//---------------------------------------------------------------------------
void BLRMatchMap::writePath(const SDGString& filename,
		std::list<RangePairSet>& rps_list, int verbose)
{
  if(verbose>0)
    std::cout<<"writing 'path' file..."<<std::flush;

  std::ofstream fout(filename);
  std::ofstream foutRpsAttr(filename + ".attr");
  // TODO DO NOT delete this lines, figure out how to ...
  writeRpsList(rps_list, fout);
  writeRpsListAttribute(rps_list, foutRpsAttr);
  
//  unsigned path_id=0;
//  bool mergedS = false;
//  for(MapPath::iterator iter_hash=path_begin();
//      iter_hash!=path_end();iter_hash++)
//    {
//      std::string query_name=num2nameQ[iter_hash->first.first];
//      std::string subject_name;
//
//      if (iter_hash->first.second == -1)
//	mergedS = true;
//
//      if(same_db)
//	subject_name=num2nameQ[iter_hash->first.second];
//      else
//	subject_name=num2nameS[iter_hash->first.second];
//
//      for(std::list<RangePairSet>::iterator iter_list
//	    =iter_hash->second.begin();iter_list!=iter_hash->second.end();
//	  iter_list++)
//	{
//		unsigned id = ++path_id;
//		if (mergedS)
//			writePathForMergedS(fout, id, query_name, iter_list->getPath());
//		else
//			iter_list->write(fout,id,query_name,subject_name);
//			iter_list->writeRpsAttr(foutRpsAttr,id,query_name,subject_name);
//
//	}
//    }
  if(verbose>0)
    std::cout<<" done"<<std::endl;
}
//---------------------------------------------------------------------------
void BLRMatchMap::writeRpsListAttribute(std::list<RangePairSet>& rps_list, std::ostream& out)
{
  unsigned path_id=0;
  for(std::list<RangePairSet>::iterator iter_list
	    =rps_list.begin();iter_list!=rps_list.end();
	  iter_list++)
    {
      std::string query_name=num2nameQ[iter_list->getNumQuery()];
      std::string subject_name;

      if(same_db)
    	  subject_name=num2nameQ[iter_list->getNumSubject()];
      else
    	  subject_name=num2nameS[iter_list->getNumSubject()];
      unsigned id = ++path_id;
      iter_list->writeRpsAttr(out,id,query_name,subject_name);
    }
}
//---------------------------------------------------------------------------
void BLRMatchMap::writeRpsList(std::list<RangePairSet>& rps_list, std::ostream& out)
{
  unsigned path_id=0;
  for(std::list<RangePairSet>::iterator iter_list
	    =rps_list.begin();iter_list!=rps_list.end();
	  iter_list++)
    {
      std::string query_name=num2nameQ[iter_list->getNumQuery()];
      std::string subject_name;
      if(same_db)
    	  subject_name=num2nameQ[iter_list->getNumSubject()];
      else
    	  subject_name=num2nameS[iter_list->getNumSubject()];
      unsigned id = ++path_id;
      iter_list->write(out,id,query_name,subject_name);
    }
}
//---------------------------------------------------------------------------
void BLRMatchMap::writeBED(const SDGString& filename, const std::list<RangePairSet>& rps_list, const SDGString& color, int verbose)
{
	std::ostringstream bedStream;	
	writeBED(bedStream, rps_list, color, verbose);	
	std::ofstream bedFile(filename);
	bedFile<<bedStream.str();
}
//---------------------------------------------------------------------------
void BLRMatchMap::writeBED(std::ostream& out, const std::list<RangePairSet>& rps_list, const SDGString& color, int verbose)
{
   if(verbose>0)
	   std::cout<<"writing 'bed' file..."<<std::flush;

   // TODO: debug
   std::cout<<" "<<std::endl;
   std::cout<<"writeBED rpsList size "<<rps_list.size()<<std::endl;

   for(std::list<RangePairSet>::const_iterator iter_list
	    =rps_list.begin();iter_list!=rps_list.end();
	  iter_list++)
    {
	std::string query_name=num2nameQ[iter_list->getNumQuery()];
      	std::string subject_name;

      	RangePairSet rps = *iter_list;
      	const std::list<RangePair> rps_path = rps.getPath();
      
      	if (! rps_path.empty()){       
		std::vector<SDGString> vec_subjects;
      		for (std::list<RangePair>::const_iterator it_path = rps_path.begin(); it_path != rps_path.end(); it_path++){
			SDGString tmp_subject_name;
			if (same_db)
	      			tmp_subject_name = num2nameQ[it_path->getRangeS().getNumChr()];
      	     		else
	       			tmp_subject_name = num2nameS[it_path->getRangeS().getNumChr()];
			vec_subjects.push_back(tmp_subject_name);
      		}
	        
		std::stringstream tmp_stream;
		for(size_t i = 0; i < vec_subjects.size(); ++i)
		{
		  if(i != 0)
		    tmp_stream << ",";
		  tmp_stream << vec_subjects[i];
		}
		subject_name = tmp_stream.str();	

      	}else{
      		if(same_db)
    	  		subject_name=num2nameQ[iter_list->getNumSubject()];
      		else
    	  		subject_name=num2nameS[iter_list->getNumSubject()];
    	}
	iter_list->writeBED(out,query_name,subject_name, color);
    }
}
//---------------------------------------------------------------------------
void BLRMatchMap::writePathForMergedS(std::ostream &fout, unsigned path_id, std::string query_name, std::list<RangePair> path)
{
		for(std::list<RangePair>::iterator i=path.begin();i!=path.end();i++)
		{
			
			long numChrS = i->getRangeS().getNumChr();

			std::string nameS = num2nameS[numChrS];		
			fout<<path_id<<"\t"<<query_name
			<<"\t"<<i->getRangeQ().getStart()
			<<"\t"<<i->getRangeQ().getEnd()
			<<"\t"<<nameS
			<<"\t"<<i->getRangeS().getStart()
			<<"\t"<<i->getRangeS().getEnd()
			<<"\t"<<i->getE_value()
			<<"\t"<<i->getScore()
			<<"\t"<<i->getIdentity()
			<<std::endl;
		}
}

//---------------------------------------------------------------------------
RangeMap BLRMatchMap::writeMap(const SDGString& filename, int verbose)
{
  if(verbose>0)
    std::cout<<"writing 'map' file..."<<std::flush;

  RangeMap matchmap;
  unsigned path_id=0;
  for(MapPath::iterator iter_hash=path_begin();
      iter_hash!=path_end();iter_hash++)
    {
      std::string query_name=num2nameQ[iter_hash->first.first];
      std::string subject_name;
      if(same_db)
	subject_name=num2nameQ[iter_hash->first.second];
      else
	subject_name=num2nameS[iter_hash->first.second];

      std::string subseqname(subject_name,0,subject_name.find(" "));



      for(std::list<RangePairSet>::iterator iter_list
	    =iter_hash->second.begin();iter_list!=iter_hash->second.end();
	  iter_list++)
	{
	  RangeAlignSet rasQ=iter_list->getRangeAlignSetQ();
	  RangeAlignSet rasS=iter_list->getRangeAlignSetS();

	  SDGString copyname=subseqname+"."+SDGString(++path_id);
	  
	  if(!rasS.isPlusStrand())
	  	rasQ.reverse();
	  
	matchmap.add(RangeSeq(rasQ,copyname,query_name));
	}
    }

  matchmap.save(filename);
  if(verbose>0)
    std::cout<<" done"<<std::endl;
  return matchmap;
}
//---------------------------------------------------------------------------
void BLRMatchMap::writeSeq(const RangeMap& matchmap,const SDGString& filename, int verbose)
{
  if(verbose>0)
    std::cout<<"writing 'fasta' file..."<<std::endl;
  matchmap.writeSeq(filename,query_db);
  if(verbose>0)
    std::cout<<" done"<<std::endl;
}
//---------------------------------------------------------------------------
void BLRMatchMap::writeMapAlign(std::ostream& out)
{
  for(MapAlign::iterator m=map_align.begin(); m!=map_align.end();m++)
  {
      while (!m->second.empty())
      {
        Key key = m->first;
        SDGString numQuery = SDGString(key.first);
        SDGString numSubject  = SDGString(key.second);
        //std::cout<<" "<<std::endl;
        //std::cout<<"Num query: "<<numQuery<<" Num subject: "<<numSubject<<std::endl;
        RangePair range_pair = m->second.back();
        SDGString numChr = SDGString(range_pair.getRangeQ().getNumChr());
        SDGString nameSeq = range_pair.getRangeQ().getNameSeq();
        //std::cout<<"Query Name Seq: "<<nameSeq<<"Query Num chr: "<<numChr<<std::endl; 
	      m->second.pop_back();
        range_pair.writetxt(out);
      }
  }

}
//---------------------------------------------------------------------------
void BLRMatchMap::contigOverlap(void)
{
  if(para->getQuery()!=para->getBank())
    {
      std::cout<<"MATCHER "<<para->getQuery()<<" vs "<<para->getBank()
	  <<" -> no contig overlap !!"<<std::endl;
      return;
    }


  unsigned count=0;
  RangeMap matchmap;
  MapPath::iterator iter_hash=map_path.begin();
  while (iter_hash!=map_path.end())
    {
      std::list<RangePairSet>::iterator iter_list=iter_hash->second.begin();
      while (iter_list!=iter_hash->second.end())
	{
	  if(
	     (
	      iter_list->getRangeQ().getMin()==1
	      && iter_list->getRangeS().getMax()
	      ==query_db[(iter_list->getRangeS().getNumChr())-1].length()
	      )
	     ||
	     (
	      iter_list->getRangeS().getMin()==1
	      && iter_list->getRangeQ().getMax()
	      ==query_db[(iter_list->getRangeQ().getNumChr())-1].length()
	      )
	     )
	    {
	      std::cout<<*iter_list<<std::endl;

	      SDGString name(++count);
	      std::string chr_name=num2nameQ[(iter_list->getRangeQ().getNumChr())-1];

	      if(iter_list->getRangeQ().isPlusStrand())
		matchmap.add(RangeSeq(name,chr_name,
				      iter_list->getRangeQ().getStart(),
				      iter_list->getRangeQ().getEnd()));
	      else
		matchmap.add(RangeSeq(name,chr_name,
				      iter_list->getRangeQ().getEnd(),
				      iter_list->getRangeQ().getStart()));

	      chr_name=num2nameS[(iter_list->getRangeS().getNumChr())-1];

	      if(iter_list->getRangeS().isPlusStrand())
		matchmap.add(RangeSeq(name,chr_name,
				      iter_list->getRangeS().getStart(),
				      iter_list->getRangeS().getEnd()));
	      else
		matchmap.add(RangeSeq(name,chr_name,
				      iter_list->getRangeS().getEnd(),
				      iter_list->getRangeS().getStart()));
	    }

	  iter_list++;
	}
      iter_hash++;
    }
}
//---------------------------------------------------------------------------
unsigned BLRMatchMap::getNbSseq(void)
{
	if( same_db )
		return num2nameQ.size();
	else
		return num2nameS.size();
}
//---------------------------------------------------------------------------
unsigned BLRMatchMap::getNbMatchesInMapAlign(void)
{
  unsigned nbMatches=0;
  for(MapAlign::iterator m=map_align.begin(); m!=map_align.end(); m++)
    nbMatches += m->second.size();
  return nbMatches;
}
//---------------------------------------------------------------------------
unsigned BLRMatchMap::getNbMatchesInMapPath(void)
{
  unsigned nbMatches=0;
  for(MapPath::iterator m=map_path.begin(); m!=map_path.end(); m++)
    for(std::list<RangePairSet>::iterator i=m->second.begin(); i!=m->second.end(); i++)
      nbMatches += i->getNbRangePairs();
  return nbMatches;
}
//---------------------------------------------------------------------------
unsigned BLRMatchMap::getNbDistinctPaths(void)
{
  unsigned nbPaths=0;
  for(MapPath::iterator m=map_path.begin(); m!=map_path.end(); m++)
    nbPaths += m->second.size();
  return nbPaths;
}
//---------------------------------------------------------------------------
void BLRMatchMap::computeScoreWithLength(std::list<RangePairSet>& rpsList)
{
    for(std::list<RangePairSet>::iterator i=rpsList.begin(); i!=rpsList.end(); i++)
	{
		i->computeScoreWithLength();
	}
}

//---------------------------------------------------------------------------
void BLRMatchMap::computeScoreWithLength()
{
      for (MapPath::iterator m=map_path.begin(); m!=map_path.end(); m++){
	computeScoreWithLength(m->second);
      }
}
//---------------------------------------------------------------------------
void BLRMatchMap::merge(int verbose)
{
	// Assume that mapPath content is set in rpsList
	Graph<unsigned long> graph = clusterizeOverlapingRps(rpsList, verbose-1);	
	if (verbose > 1){
		std::cout<<" "<<std::endl;
		graph.view();
	}

	std::list<RangePairSet> rpsListAfterMerge = mergeOnCluster(rpsList, graph, verbose -1); 
	
	for (std::list<RangePairSet>::iterator it = rpsListAfterMerge.begin(); it != rpsListAfterMerge.end(); it++){
		insert_path(*it);
	}
}
//---------------------------------------------------------------------------
std::list<RangePairSet> BLRMatchMap::mergeOnCluster(std::list<RangePairSet> rpsList, Graph <unsigned long>& graph, int verbose)
{
	std::list<RangePairSet> mergedRpsList;
	
	// index Rps with id
	std::map<unsigned long, RangePairSet> idToRps;
 	for(std::list<RangePairSet>::iterator lrp_it=rpsList.begin(); lrp_it!=rpsList.end();lrp_it++){
		unsigned long id = lrp_it->getId();
		idToRps[id] = *lrp_it;
	}

	std::vector<std::vector<unsigned long> > vec;
	graph.connexComp(vec);	
	
	// merge rps of in each connex comp
        for(std::vector< std::vector<unsigned long> >::iterator it1 =vec.begin();it1!=vec.end();it1++){
		std::vector<unsigned long> currentConnexComp = *it1;
		unsigned size = currentConnexComp.size();
		RangePairSet firstRps = idToRps[currentConnexComp[0]];

		for (unsigned i=1; i<size; i++){
			RangePairSet currentRps = idToRps[currentConnexComp[i]];
			if (verbose > 1){
				std::cout<<" "<<std::endl;
			 	currentRps.view();	
			}
			firstRps.mergeQ(currentRps);
			firstRps.orientSubjects();
		}

		mergedRpsList.push_back(firstRps);
	}	
	return mergedRpsList;
}
//---------------------------------------------------------------------------
Graph <unsigned long> BLRMatchMap::clusterizeOverlapingRps(const std::list<RangePairSet>& rpsList, int verbose)
{
	Graph<unsigned long> graph;
	
  	for(std::list<RangePairSet>::const_iterator lrp_it1=rpsList.begin(); lrp_it1!=rpsList.end();lrp_it1++){
		RangePairSet rps1 = *lrp_it1;
		graph.add_node(rps1.getId());
  		for(std::list<RangePairSet>::const_iterator lrp_it2=rpsList.begin(); lrp_it2!=rpsList.end();lrp_it2++){
			RangePairSet rps2 = *lrp_it2;
			graph.add_node(rps2.getId());
			if (rps1 != rps2 && rps1.first.getNumChr() == rps2.first.getNumChr() && rps1.overlapQ_length(rps2) >= 200){
				graph.add_edge(rps1.getId(),rps2.getId());
			}
		}
	}
	return graph;
}
//---------------------------------------------------------------------------
void BLRMatchMap::insert(RangePairSet& rp)
{
  	unsigned long currentId, previousId=0;
        RangePairSet rps;
	std::list<RangePair> path;

	if (rpsList.size() != 0)
	{
		rps = rpsList.back();
		previousId = rps.getId();
	}

	currentId = rp.getId();
	if (currentId == previousId)
	{	
		rps = rpsList.back();
		rpsList.pop_back();
		path = rps.getPath();
       	} 

	path.push_back(rp);
        rps.setPath(path,para->getDist_pen(),0.0,para->getGap_pen());
	rpsList.push_back(rps);

}
//----------------------------------------------------------------------------
void BLRMatchMap::mapPathJoinAndComputeScoreWithLengthOnly(bool joining, bool clean_before, bool clean_after, int verbose)
{
  map_path.clear();

  if(joining)
    {
      if(verbose>0)
	std::cout<<"Join parameters: dist_pen="<<para->getDist_pen()
		 <<" gap_pen="<<para->getGap_pen()
		 <<" overlap="<<para->getOverlap()<<std::endl<<std::flush;
      for(MapAlign::iterator m=map_align.begin(); m!=map_align.end();m++)
	{
	  FragAlign fragAlign(para->getDist_pen(),0,para->getGap_pen(),
			      para->getOverlap());
	  map_path[m->first]=fragAlign.join(m->second);
	  m->second.clear();
	}
      if(verbose>0){
	std::cout<<"nb of matches: "<<getNbMatchesInMapPath()<<std::endl;
	std::cout<<"nb of paths: "<<getNbDistinctPaths()<<std::endl;}
      
      if( ( clean_before | clean_after ) & ( verbose > 0 ) )
    	  std::cout<<"Clean the connections..."<<std::endl<<std::flush;
      if(clean_before)
	clean_path(true,verbose-1); // clean when same subject
      if(clean_after)
	clean_path(false,verbose-1); // clean with all subjects
      if( ( clean_before | clean_after ) & ( verbose > 0 ) )
      {
    	  std::cout<<"nb of matches: "<<getNbMatchesInMapPath()<<std::endl;
    	  std::cout<<"nb of paths: "<<getNbDistinctPaths()<<std::endl;
    	  std::cout<<"Connections were cleaned when necessary."<<std::endl;
      }
      if(verbose>0){
	std::cout<<"Recompute score with length... "<<getNbMatchesInMapPath()<<std::endl;
      }
     
      for (MapPath::iterator m=map_path.begin(); m!=map_path.end(); m++){
	computeScoreWithLength(m->second);
      }	 
    }
  return;
}
// TODO remove
void BLRMatchMap::viewMapPath(){
  for(BLRMatchMap::MapPath::iterator m=path_begin(); m!=path_end();m++)
  {
      while (!m->second.empty())
      {
        BLRMatchMap::Key key = m->first;
        SDGString keyFirst = SDGString(key.first);
        SDGString keySecond  = SDGString(key.second);
        std::cout<<" "<<std::endl;
        std::cout<<"("<<keyFirst<<","<<keySecond<<")"<<std::endl;
        std::cout<<" "<<std::endl;
        RangePairSet range_pair = m->second.back();
        SDGString numChr = SDGString(range_pair.getRangeQ().getNumChr());
        SDGString nameSeq = range_pair.getRangeQ().getNameSeq();
	      m->second.pop_back();
        range_pair.view();
      }
  }
}
