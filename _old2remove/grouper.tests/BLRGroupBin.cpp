/**
 * \file BLRGroupBin.cpp
 */

#include <cmath>
#include <map>

#include "BLRGroupBin.h"


void BLRGroup::addMember( GROUPLIST::iterator& gi, Member& mem )
{
	mem.id = count_memb;
	vec_memb.push_back( mem );
	gi->insert( gi->begin(), mem.id );
	memIdx->insert_member_idx( mem, gi );
	count_memb++;
}


void BLRGroup::addGroup( Member& memQ, Member& memS )
{
	std::list<unsigned> g;
	GROUPLIST::iterator group_it=group_list.insert( group_list.begin(), g );
	addMember( group_it, memQ );
	addMember( group_it, memS );
}


void BLRGroup::mergeGroup( GROUPLIST::iterator& iQ, GROUPLIST::iterator& iS )
{
	for( std::list<unsigned>::iterator i=iS->begin(); i!=iS->end(); i++ )
		memIdx->set_groupit_member_idx( *i, iQ );
	iQ->splice( iQ->begin(), *iS, iS->begin(), iS->end() );
	group_list.erase( iS );
}


/*---------------------------------------------------------------------
  - Function merge two Member
  - Arguments: a-> reference on current alignment
  -            b-> reference on current  Member
  - c-> reference on current complement alignment a.S<-a.Q || a.Q<-a.S.
  ----------------------------------------------------------------------*/
void BLRGroup::mergeMember( GROUPLIST::iterator& ga, unsigned ma,
		GROUPLIST::iterator& gb, unsigned mb,
		GROUPLIST::iterator& gc )
{
  // Case when c is assigned but not a
  // roam_group must assign a group at current Member,
  // function Compare must manage this:
	if( ga==group_list.end() )
	{
		std::cout<<"###VERIFIER COMPARE QUAND Si==group_list.end()###";
		return;
	}
	// merge if a and b is in different group
	if( ga != gb )
	{
		for( std::list<unsigned>::iterator i=gb->begin(); i!=gb->end(); i++ )
			memIdx->set_groupit_member_idx( *i, ga );
		ga->splice( ga->begin(), *gb, gb->begin(), gb->end() );
		if( gc == gb )
			gc = ga; // To change gc before erase gb
		group_list.erase( gb );
		gb = ga;
	}
	// Normally case where current Member (a) or/and this complement (c) found group before
	// merge if a and b is in the same group
	if( ga==gb && ma!=mb )
	{
		vec_memb[ma].merge( vec_memb[mb] );
		vec_memb[ma].idlist.splice( vec_memb[ma].idlist.end(), vec_memb[mb].idlist );
		if( vec_memb[ma].getMin() >= memIdx->get_start_member_idx( ma )
				&& vec_memb[ma].getMax() <= memIdx->get_end_member_idx( ma ) )
			return;
		memIdx->adjust_member_idx( ma, vec_memb[ma].getMin(), vec_memb[ma].getMax() );
		memIdx->erase_member_idx( mb );
		ga->remove( mb );
		mb = ma;
	}
}


void BLRGroup::mergeAlign( RangeAlignSet& align, long long id, unsigned member_id )
{
	vec_memb[ member_id ].idlist.push_back( id );
	bool isMemberIncluded = vec_memb[ member_id ].getIncluded();
	bool isAlignIncluded = align.getIncluded();
	vec_memb[ member_id ].merge( align );
	if( ! isMemberIncluded && isAlignIncluded && vec_memb[ member_id ].hasSameRanges( align ) )
		vec_memb[ member_id ].setIncluded( false );
	if( align.getMin() >= memIdx->get_start_member_idx( member_id )
			&& align.getMax() <= memIdx->get_end_member_idx( member_id ) )
		return;
	memIdx->adjust_member_idx( member_id, vec_memb[member_id].getMin(), vec_memb[member_id].getMax() );
}


void BLRGroup::Inverse( GROUPLIST::iterator& g )
{
	for( std::list<unsigned>::iterator member_iter=g->begin();
	member_iter!=g->end();
	member_iter++ )
		vec_memb[ *member_iter ].reverse();
}


bool BLRGroup::compareQ( RangeAlignSet& alignQ, RangeAlignSet& alignS,
		long long id,
		unsigned mi,
		GROUPLIST::iterator& Qi, GROUPLIST::iterator& Si,
		GROUPLIST::iterator& Gi,
		int verbose )
{
	id=id*(-1);

	Member& mb = vec_memb[mi];

	if( alignQ.getNumChr() != mb.getNumChr() )
		return false;

	// fit member to query if they are overlapping, and modify strand if necessary
	if( mb.overlap( alignQ ) )
	{
		if( alignQ.isContained( mb ) )
			alignQ.setIncluded( true );

		unsigned lenQ = mb.overlap_length( alignQ );
		unsigned lenM = std::max( mb.getLengthSet(), alignQ.getLengthSet() );
		double coverQ = (double)(lenQ) / lenM;
		unsigned diffAlignQ = alignQ.getLengthSet() - lenQ;
		unsigned diffMb = mb.getLengthSet() - lenQ;
		if(verbose>0)
		{
			if( cover_limit <= 1 )
				std::cout<<"coverage: "<<std::floor(100*coverQ)<<"% ("<<lenQ<<"/"<<lenM<<")"<<std::endl;
			else
				std::cout<<"coverage: diffAlignQ="<<diffAlignQ<<" diffMb="<<diffMb<<""<<std::endl;
		}
		if( cover_limit <= 1 && coverQ < cover_limit )
			return false;
		if( cover_limit > 1 && ( diffAlignQ > cover_limit || diffMb > cover_limit ) )
			return false;

		if( mb.isPlusStrand() != alignQ.isPlusStrand() )
		{
			if( Qi==group_list.end() && Si==group_list.end() )
			{
				// alignment did not find a group before
				// could change strand of subject and  query.
				alignQ.reverse();
				alignS.reverse();
			}
			else if( Si!=Gi && Qi!=Gi )
			{
				// query or subject find another group before !!
				// to push it into this group, strand of other group's Member
				// must be reversed
				Inverse( Gi );
			}
		}
		mergeAlign( alignQ, id, mi );
		return true;
	}
	return false;
}


bool BLRGroup::compareS( RangeAlignSet& alignQ, RangeAlignSet& alignS,
		long long id,
		unsigned mi,
		GROUPLIST::iterator& Qi, GROUPLIST::iterator& Si,
		GROUPLIST::iterator& Gi,
		int verbose )
{
	Member& mb = vec_memb[mi];

	if( alignS.getNumChr() != mb.getNumChr() )
		return false;

	// fit member to subject if there are overlapping, and modify strand if necessary
	if( mb.overlap( alignS ) )
	{
		if( alignS.isContained( mb ) )
			alignS.setIncluded( true );

		unsigned lenS = mb.overlap_length( alignS );
		unsigned lenM = std::max( mb.getLengthSet(), alignS.getLengthSet() );
		double coverS = (double)(lenS) / lenM;
		unsigned diffAlignS = alignS.getLengthSet() - lenS;
		unsigned diffMb = mb.getLengthSet() - lenS;
		if(verbose>0)
		{
			if( cover_limit <= 1 )
				std::cout<<"coverage: "<<std::floor(100*coverS)<<"% ("<<lenS<<"/"<<lenM<<")"<<std::endl;
			else
				std::cout<<"coverage: diffAlignS="<<diffAlignS<<" diffMb="<<diffMb<<""<<std::endl;
		}
		if( cover_limit <= 1 && coverS < cover_limit )
			return false;
		if( cover_limit > 1 && ( diffAlignS > cover_limit || diffMb > cover_limit ) )
			return false;

		if( mb.isPlusStrand() != alignS.isPlusStrand() )
		{
			if( Qi==group_list.end() && Si==group_list.end() )
			{
				// alignment did not find a group before
				// could change strand of subject and  query.
				alignS.reverse();
				alignQ.reverse();
			}
			else if( Si!=Gi && Qi!=Gi )
			{
				// query or subject find another group before !!
				// to push it into this group, strand of other group's Member
				// must be reversed
				Inverse( Gi );
			}
		}
		mergeAlign( alignS, id, mi );
		return true;
	}
	return false;
}


void BLRGroup::roam_group( RangePairSet& align, int verbose )
{
	long long id = align.getId();

	RangeAlignSet alignQ( align.getRangeAlignSetQ() );
	RangeAlignSet alignS( align.getRangeAlignSetS() );

	GROUPLIST::iterator group_iter = group_list.end();
	GROUPLIST::iterator group_iterQ = group_list.end();
	GROUPLIST::iterator group_iterS = group_list.end();
	unsigned member_idQ=0, member_idS=0 ;

	if(verbose>0)
		std::cout<<"query ("<<alignQ.getRangeSetSize()<<"): "
		<<match_map.numQ2name( alignQ.getNumChr() )<<" "
		<< alignQ.getStart()<<" "<<alignQ.getEnd()
		<<std::endl<<std::flush;
	if(verbose>1 && alignQ.getRangeSetSize()>1)
		alignQ.view();
	std::vector<unsigned> idQ = memIdx->search_member( alignQ.getNumChr(),
			alignQ.getMin(),
			alignQ.getMax() );
	if(verbose>0)
		std::cout<<"nb overlapping members: "<<idQ.size()<<std::endl;
	for( std::vector<unsigned>::iterator i=idQ.begin();
	i!=idQ.end(); i++ )
	{
		unsigned member_id=*i;
		group_iter = memIdx->get_groupit_member_idx( member_id );
		if( compareQ( alignQ, alignS, id, member_id, group_iterQ, group_iterS,
				group_iter, verbose-1 ) )
		{
			// the alignment query overlaps with a member
			if( group_iterQ == group_list.end() )
			{
				// first group found
				if(verbose>1)
					std::cout<<"first group found"<<std::endl;
				group_iterQ = group_iter;
				member_idQ = member_id;
			}
			else
			{
				// query found a member before
				if(verbose>1)
					std::cout<<"merge member"<<std::endl;
				mergeMember( group_iterQ, member_idQ, group_iter,
						member_id, group_iterS );
			}
		}
	}
	if(verbose>1)
		std::cout<<"included: "<<std::boolalpha<<alignQ.getIncluded()<<std::noboolalpha<<std::endl<<std::flush;

	if(verbose>0)
	{
		std::string subjectName;
		if( grouper_parameter->getBank() == grouper_parameter->getQuery() )
			subjectName = match_map.numQ2name( alignS.getNumChr() );
		else
			subjectName = match_map.numS2name( alignS.getNumChr() );
		std::cout<<"subject ("<<alignS.getRangeSetSize()<<"): "
		<<subjectName<<" "<< alignS.getStart()<<" "<<alignS.getEnd()
		<<std::endl<<std::flush;
	}
	if(verbose>1 && alignS.getRangeSetSize()>1)
		alignS.view();
	std::vector<unsigned> idS=memIdx->search_member(alignS.getNumChr(),
			alignS.getMin(),
			alignS.getMax());
	if(verbose>0)
		std::cout<<"nb overlapping members: "<<idS.size()<<std::endl;
	for( std::vector<unsigned>::iterator i=idS.begin();
	i!=idS.end(); i++ )
	{
		unsigned member_id=*i;
		group_iter = memIdx->get_groupit_member_idx( member_id );
		if( compareS( alignQ, alignS, id, member_id, group_iterQ, group_iterS,
				group_iter, verbose-1 ) )
		{
			// the alignment subject overlaps with a member
			if( group_iterS == group_list.end() )
			{
				// first group found
				group_iterS = group_iter;
				member_idS = member_id;
			}
			else
			{
				// subject found a member before
				if(verbose>1)
					std::cout<<"merge member"<<std::endl;
				mergeMember( group_iterS, member_idS, group_iter,
						member_id, group_iterQ );
			}
		}
	}
	if(verbose>1)
		std::cout<<"included: "<<std::boolalpha<<alignS.getIncluded()<<std::noboolalpha<<std::endl<<std::flush;

	build_group( alignQ, alignS, id, group_iterQ, group_iterS, verbose );
}

//void BLRGroup::searchOverlap(void *args)
//{
//	struct searchOverlapArgs *p = (searchOverlapArgs*)args;
//	std::vector<unsigned> vect;
//	pthread_mutex_t m_lock;
//	pthread_mutex_init(&m_lock, NULL);
//
//	vect = p->memIdx->search_member( p->ras.getNumChr(),
//			p->ras.getMin(),
//			p->ras.getMax() );
//
//	pthread_mutex_lock(&m_lock);
//	*(p->p_vid)=vect;
//	pthread_mutex_unlock(&m_lock);
//};

//void BLRGroup::roam_group_threads( RangePairSet& align, int verbose )
//{
//	long long id = align.getId();
//
//	RangeAlignSet alignQ( align.getRangeAlignSetQ() );
//	RangeAlignSet alignS( align.getRangeAlignSetS() );
//
//	GROUPLIST::iterator group_iter = group_list.end();
//	GROUPLIST::iterator group_iterQ = group_list.end();
//	GROUPLIST::iterator group_iterS = group_list.end();
//	unsigned member_idQ=0, member_idS=0 ;
//
//	if(verbose>0)
//		std::cout<<"query ("<<alignQ.getRangeSetSize()<<"): "
//		<<match_map.numQ2name( alignQ.getNumChr() )<<" "
//		<< alignQ.getStart()<<" "<<alignQ.getEnd()
//		<<std::endl<<std::flush;
//	if(verbose>1 && alignQ.getRangeSetSize()>1)
//		alignQ.view();
//
//
//	//pthread version
//	std::vector<unsigned> idQ,idS;
//	searchOverlapArgs* argsQ=new searchOverlapArgs(memIdx,alignQ,idQ);
//    searchOverlapArgs* argsS=new searchOverlapArgs(memIdx,alignS,idS);
//
//	if(grouper_parameter->getNbThread()==1)
//	{
//		searchOverlap(argsQ);
//		searchOverlap(argsS);
//	}
//	else
//	{
//		Task* tQ = new Task(&searchOverlap, (void*) argsQ);
//		tp.add_task(tQ);
//		Task* tS = new Task(&searchOverlap, (void*) argsS);
//		tp.add_task(tS);
//	}
//	delete argsQ;
//	delete argsS;
//	//end pthread version
//
//	if(verbose>0)
//		std::cout<<"nb overlapping members: "<<idQ.size()<<std::endl;
//	for( std::vector<unsigned>::iterator i=idQ.begin();
//	i!=idQ.end(); i++ )
//	{
//		unsigned member_id=*i;
//		group_iter = memIdx->get_groupit_member_idx( member_id );
//		if( compareQ( alignQ, alignS, id, member_id, group_iterQ, group_iterS,
//				group_iter, verbose-1 ) )
//		{
//			// the alignment query overlaps with a member
//			if( group_iterQ == group_list.end() )
//			{
//				// first group found
//				if(verbose>1)
//					std::cout<<"first group found"<<std::endl;
//				group_iterQ = group_iter;
//				member_idQ = member_id;
//			}
//			else
//			{
//				// query found a member before
//				if(verbose>1)
//					std::cout<<"merge member"<<std::endl;
//				mergeMember( group_iterQ, member_idQ, group_iter,
//						member_id, group_iterS );
//			}
//		}
//	}
//	if(verbose>1)
//		std::cout<<"included: "<<std::boolalpha<<alignQ.getIncluded()<<std::noboolalpha<<std::endl<<std::flush;
//
//	if(verbose>0)
//	{
//		std::string subjectName;
//		if( grouper_parameter->getBank() == grouper_parameter->getQuery() )
//			subjectName = match_map.numQ2name( alignS.getNumChr() );
//		else
//			subjectName = match_map.numS2name( alignS.getNumChr() );
//		std::cout<<"subject ("<<alignS.getRangeSetSize()<<"): "
//		<<subjectName<<" "<< alignS.getStart()<<" "<<alignS.getEnd()
//		<<std::endl<<std::flush;
//	}
//	if(verbose>1 && alignS.getRangeSetSize()>1)
//		alignS.view();
//
//
//	if(verbose>0)
//		std::cout<<"nb overlapping members: "<<idS.size()<<std::endl;
//	for( std::vector<unsigned>::iterator i=idS.begin();
//	i!=idS.end(); i++ )
//	{
//		unsigned member_id=*i;
//		group_iter = memIdx->get_groupit_member_idx( member_id );
//		if( compareS( alignQ, alignS, id, member_id, group_iterQ, group_iterS,
//				group_iter, verbose-1 ) )
//		{
//			// the alignment subject overlaps with a member
//			if( group_iterS == group_list.end() )
//			{
//				// first group found
//				group_iterS = group_iter;
//				member_idS = member_id;
//			}
//			else
//			{
//				// subject found a member before
//				if(verbose>1)
//					std::cout<<"merge member"<<std::endl;
//				mergeMember( group_iterS, member_idS, group_iter,
//						member_id, group_iterQ );
//			}
//		}
//	}
//	if(verbose>1)
//		std::cout<<"included: "<<std::boolalpha<<alignS.getIncluded()<<std::noboolalpha<<std::endl<<std::flush;
//
//	build_group( alignQ, alignS, id, group_iterQ, group_iterS, verbose );
//}
void BLRGroup::build_group( RangeAlignSet& alignQ, RangeAlignSet& alignS,
		long long id,
		GROUPLIST::iterator& iQ, GROUPLIST::iterator& iS,
		int verbose )
{
	Member mQ( alignQ, id*(-1) );
	Member mS( alignS, id );

	if( iQ!=group_list.end() && iS==group_list.end() )
	{
		if(verbose>0)
			std::cout<<"add member subject in group query"<<std::endl<<std::flush;
		addMember( iQ, mS );
	}
	else if( iQ==group_list.end() && iS!=group_list.end() )
	{
		if(verbose>0)
			std::cout<<"add member query in group subject"<<std::endl<<std::flush;
		addMember( iS, mQ );
	}
	else if( iQ!=group_list.end() && iS!=group_list.end() && iQ!=iS )
    {
		if(verbose>0)
			std::cout<<"merge groups"<<std::endl<<std::flush;
		mergeGroup( iQ, iS );
    }
	else if( iQ==group_list.end() && iS==group_list.end() )
    {
		if(verbose>0)
			std::cout<<"add group"<<std::endl<<std::flush;
		addGroup( mQ, mS );
    }
	else if( iQ!=group_list.end() && iS!=group_list.end() && iQ==iS )
		if(verbose>0)
			std::cout<<"nothing"<<std::endl<<std::flush;
}

std::list<RangePairSet> BLRGroup::getRpsListAfterLoad( int verbose )
{
      std::list<RangePairSet> rp_list;
     
      if (!grouper_parameter->getLoad_path()){

      	if(verbose>0)
    	  	std::cout<<"Load the matches..."<<std::endl<<std::flush;

      	match_map.load();

      	if(verbose>0)
    	  	std::cout<<"Matches were loaded."<<std::endl;

      	if(grouper_parameter->getJoin_frag())
      	{
    	  	if(verbose>0)
    		  	std::cout<<"Connect the fragments..."<<std::endl<<std::flush;
   	 	  	match_map.mapPath(true, false, false, false, 2);
    	  	if(verbose>0)
    		  	std::cout<<"Fragments were connected."<<std::endl;
      	}
      	else
    	  match_map.mapPath(false, false, false, false, 0);

      	if(verbose>0)
    	  	std::cout<<"Roam groups (coverage "<<grouper_parameter->getCoverage()<<")..."<<std::endl<<std::flush;

      	cover_limit=grouper_parameter->getCoverage();

      	for( BLRMatchMap::MapPath::iterator m=match_map.path_begin();
      	m!=match_map.path_end(); m++ )
      	{
    	  	while(!m->second.empty())
    	  	{
    		  	RangePairSet rp=m->second.back();
    		  	m->second.pop_back();
    		  	rp_list.push_back(rp);
    	  	}
      	}
      }else{
     	if(verbose>0)
    	  	std::cout<<"Load the matches from path file..."<<std::endl<<std::flush;
	match_map.loadPath();
	rp_list = match_map.getRpsList();
      }

      memIdx=new BLRMemIdx(match_map.getNbQseq());
      match_map.clear();
      
      rp_list.sort(RangePair::greaterLengthQ);
      return rp_list;
}


void BLRGroup::group( int verbose )
{
	if( verbose > 0 )
		std::cout<<"Group building..."<<std::endl<<std::flush;

	unsigned long nbChains=0;  // a chain is a set of connected matches

    std::list<RangePairSet> rp_list;
     
    if (!grouper_parameter->getLoad_path())
    	{
			if(verbose>0)
				std::cout<<"Load the matches..."<<std::endl<<std::flush;
			match_map.load( verbose );
			if(verbose>0)
				std::cout<<"Matches were loaded."<<std::endl;

			if( grouper_parameter->getJoin_frag() )
				{
					if(verbose>0)
						std::cout<<"Connect the fragments..."<<std::endl<<std::flush;
					match_map.mapPath(true, false, false, false, verbose-1);
					//match_map.mapPathWithThreads(true, false, false, false, verbose-1);
			if (verbose>0)
				match_map.writePath(grouper_parameter->getPrefixFileName() +  ".gpath");
			if(verbose>0)
	    		  	std::cout<<"Fragments were connected."<<std::endl;
	      		}
	      	else
	    	  match_map.mapPath(false, false, false, false, 0);
	    	  
			unsigned long id = 1;
		 	for( BLRMatchMap::MapPath::iterator m=match_map.path_begin();
		      	m!=match_map.path_end(); m++ )
		      	{
		    	  	while(!m->second.empty())
		    	  	{
		    		  	RangePairSet rp=m->second.front();
		    		  	m->second.pop_front();
					rp.setId(id);
		    		  	rp_list.push_back(rp);
					id++;
		    	  	}
		      	}
    	}
    else
    	{
     		if(verbose>0)
    	  		std::cout<<"Load the matches from path file..."<<std::endl<<std::flush;
			//match_map.loadPath();
			match_map.loadPath(grouper_parameter->getPath_filename(),1);
			rp_list = match_map.getRpsList();
      	}
      	
	cover_limit = grouper_parameter->getCoverage();
	if(verbose>0)
	{
		std::cout<<"Roam groups (coverage limit: ";
		if( cover_limit <= 1 )
			std::cout<<std::floor(100*cover_limit)<<"%)..."<<std::endl<<std::flush;
		else
			std::cout<<cover_limit<<" bp)..."<<std::endl<<std::flush;
	}
	
	memIdx = new BLRMemIdx( match_map.getNbQseq() );
	match_map.clear();

	rp_list.sort( RangePair::greaterScoreIdLenCoord );
	
	//for debuging
	if (verbose > 0)
		{
      	if (!grouper_parameter->getLoad_path())
      		writePath(grouper_parameter->getPrefixFileName() + ".rpsListNoLoadPath",rp_list,1);
      	else
      		writePath(grouper_parameter->getPrefixFileName() + ".rpsListLoadPath",rp_list,1);
      } 	
	std::list<RangePairSet>::iterator r=rp_list.begin();

	std::string chainfile = grouper_parameter->getPrefixFileName() + ".chains.path";
	std::ofstream chainstream( chainfile.c_str() );
	std::string queryName, subjectName;

	while( r!=rp_list.end() )
	{
		nbChains++;
		if(verbose>1)
			std::cout<<"chain "<<nbChains<<" (score="<<r->getScore()<<")"<<std::endl;
		r->setId( nbChains );

		queryName = match_map.numQ2name( r->getNumQuery() );
		if( grouper_parameter->getBank() == grouper_parameter->getQuery() )
			subjectName = match_map.numQ2name( r->getNumSubject() );
		else
			subjectName = match_map.numS2name( r->getNumSubject() );
		r->write( chainstream, nbChains, queryName, subjectName );

		//roam_group_threads( *r, verbose-1 );
		roam_group( *r, verbose-1 );

		rp_list.pop_front();
		r = rp_list.begin();

		if(verbose>2)
			show_group();
	}

	chainstream.close();

	if(verbose>0)
	{
		std::cout<<"number of chains: "<<nbChains<<std::endl;
		std::cout<<"number of members: "<<getNbMembers()<<std::endl;
		std::cout<<"number of groups: "<<group_list.size()<<std::endl;
		std::cout<<"Roaming done."<<std::endl;
	}

	mergeMembersInGroups( verbose );

	if(verbose>0)
	{
		std::cout<<"number of chains: "<<nbChains<<std::endl;
		std::cout<<"number of members: "<<getNbMembers()<<std::endl;
		std::cout<<"number of groups: "<<group_list.size()<<std::endl;
	}

	if( nbChains==0 )
	{
		std::cout<<"Error when building groups: no alignment!"<<std::endl;
		exit( EXIT_FAILURE );
	}

	if( verbose > 0 )
		std::cout<<"Groups were built."<<std::endl<<std::flush;
}


void BLRGroup::mergeMembersInGroups( int verbose )
{
	if( verbose > 0 )
	std::cout<<"Merge members in groups..."<<std::endl;

	GROUPLIST::iterator it_group = group_list.begin();
	unsigned countGroups = 1;
	while( it_group != group_list.end() )
	{
		unsigned nbMembers = it_group->size();
		if( nbMembers < 2 )
		{
			++ countGroups;
			++ it_group;
			continue;
		}
		if( verbose > 1 )
			std::cout<<"group "<<countGroups<<": "<<nbMembers<<" members"<<std::endl;

		std::list<unsigned>::iterator it_idPrevMb = it_group->begin();
		Member* pt_prevMb;
		unsigned nbMerges = 0;
		std::list<unsigned>::iterator it_beforeLastMb = it_group->end();
		it_beforeLastMb--;
		while( it_idPrevMb != it_beforeLastMb )
		{
			pt_prevMb = &vec_memb[ *it_idPrevMb ];
			if( pt_prevMb->isEmpty() )
			{
				++ it_idPrevMb;
				continue;
			}
			if( verbose > 2 ) {
				std::cout<<"prevMb: "; pt_prevMb->view(); }
			std::list<unsigned>::iterator it_idCurrMb = it_idPrevMb;
			++ it_idCurrMb;
			while( it_idCurrMb != it_group->end() )
			{
				Member* pt_currMb = &vec_memb[ *it_idCurrMb ];
				if( verbose > 2 && ! pt_currMb->isEmpty() ) {
					std::cout<<"currMb: "; pt_currMb->view(); }
				unsigned long minLength = std::min( pt_prevMb->getLength(), pt_currMb->getLength() );
				if( ! pt_currMb->isEmpty() && pt_prevMb->getNumChr() == pt_currMb->getNumChr() && *pt_prevMb == *pt_currMb )
//				if( ! pt_currMb->isEmpty() && pt_prevMb->overlap_length( *pt_currMb ) >= 0.25 * minLength )
				{
					if( verbose > 2 )
						std::cout<<"merge"<<std::endl;
					++ nbMerges;
					bool isIncludedPrevMb = pt_prevMb->getIncluded();
					bool isIncludedCurrMb = pt_currMb->getIncluded();
					unsigned lengthPrevMb = pt_prevMb->getLengthSet();
					unsigned lengthCurrMb = pt_currMb->getLengthSet();
					pt_prevMb->merge( *pt_currMb );
					if( ! isIncludedPrevMb && isIncludedCurrMb && lengthPrevMb <= lengthCurrMb )
						pt_prevMb->setIncluded( true );
					pt_currMb->reset();
				}
				else
				{
					if( verbose > 2 && ! pt_currMb->isEmpty() )
						std::cout<<"overlap "<<pt_prevMb->overlap_length( *pt_currMb )<<" < "<<0.25 * minLength<<std::endl;
				}
				++ it_idCurrMb;
			}
			++ it_idPrevMb;
		}

		if( verbose > 1 && nbMerges > 0 )
		{
			std::cout<<nbMerges<<" merge(s)"<<std::endl;
			unsigned nbMembersAfter = 0;
			for( std::list<unsigned>::iterator member_iter=it_group->begin(); member_iter!=it_group->end(); member_iter++ )
				if( vec_memb[ *member_iter ].getLengthSet() > 0 )
					nbMembersAfter++;
			std::cout<<"nbMembersAfter="<<nbMembersAfter<<std::endl;
		}
		++ countGroups;
		++ it_group;
	}

	if( verbose > 0 )
		std::cout<<"Members in groups were merged."<<std::endl;
}

void BLRGroup::writePath(const SDGString& filename,std::list<RangePairSet> l, int verbose)
{
  if(verbose>0)
    std::cout<<"writing 'path' file..."<<std::flush;

  std::ofstream fout(filename);
  std::ofstream foutRpsAttr(filename + ".attr");  
 
  for(std::list<RangePairSet>::iterator iter_list
	    =l.begin();iter_list!=l.end();
	  iter_list++)
	{
		unsigned long id = iter_list->getId();
		iter_list->write(fout,id,iter_list->first.getNameSeq(),iter_list->second.getNameSeq());
		iter_list->writeRpsAttr(foutRpsAttr,id,iter_list->first.getNameSeq(),iter_list->second.getNameSeq());
	}
    
  if(verbose>0)
    std::cout<<" done"<<std::endl;
}



void BLRGroup::show_group( std::ostream& out )
{
	int count_group=1;
	int count_mb=1;
	SDGString chr_name;

	GROUPLIST::iterator group_iter=group_list.begin();

	bool samedb=false;
	if( grouper_parameter->getBank() == grouper_parameter->getQuery() )
		samedb = true;

	while( group_iter != group_list.end() )
	{
		std::list<unsigned>::iterator member_iter=group_iter->begin();
		out<<"GROUP NUMBER "<<count_group<<" CLUSTER "<<gr2clust[count_group];
		out<<" contains "<<getNbMembers( group_iter )<<" members:"<<std::endl;
		out<<"\tmember \tseq_num\tstart \tend  \tfrag \tlength \tinclude\ttotal \talign_list"<<std::endl;
		while( member_iter != group_iter->end() )
		{
			Member memb = vec_memb[ *member_iter ];
			if( memb.getLengthSet() > 0 )
			{
				if( memb.idlist.front() > 0 ) //subject
					if( samedb )
						chr_name = match_map.numQ2name( memb.getNumChr() );
					else
						chr_name = match_map.numS2name( memb.getNumChr() );
				else  // query
					chr_name = match_map.numQ2name( memb.getNumChr() );
				out<<count_group<<"\t"<<count_mb<<"\t";
				out<<chr_name<<"\t";
				out<<memb.getStart()<<"\t";
				out<<memb.getEnd()<<"\t";
				out<<memb.getRangeSet().size()<<"\t";
				out<<memb.getLengthSet()<<"\t";
				out<<std::boolalpha<<memb.getIncluded()<<"\t";
				out<<memb.idlist.size()<<"\t";
				for( std::list<long long>::iterator i=memb.idlist.begin();
				i!=memb.idlist.end(); i++ )
					out<<*i<<" ";
				out<<std::endl;
				count_mb++;
			}
			member_iter++;
		}
		count_group++;
		group_iter++;
	}

	unsigned count=0;
	for( std::vector< std::vector<unsigned> >::iterator l=clust.begin(); l!=clust.end(); l++ )
	{
		out<<"cluster #"<<++count<<" size="<<l->size()<<" :";
		for( std::vector<unsigned>::iterator i=l->begin(); i!=l->end(); i++ )
		{
			out<<" "<<*i;
		}
		out<<std::endl;
	}
	out<<"Members found="<<--count_mb<<std::endl;
	out<<"Groups found="<<--count_group<<std::endl;
	out<<"Clusters found="<<clust.size()<<std::endl;
}


void BLRGroup::save( int verbose )
{
	if( verbose > 0 )
		std::cout<<"Write the results..."<<std::endl<<std::flush;

	unsigned long count_gp;            // Number of groups
	unsigned long count_mb;            // Number of group members
	unsigned long count_all;           // Total number of members
	GROUPLIST::iterator group_iter;    // counter of loop to visit group_list
	std::list<unsigned>::iterator member_iter; // counter of loop to visit members of group

	RangeMap membmap;

	count_gp=1;
	count_all=0;

	bool samedb = false;
	if( grouper_parameter->getBank() == grouper_parameter->getQuery() )
		samedb = true;

	group_iter = group_list.begin();
	while( group_iter != group_list.end() )
	{
		if(verbose>1)
			std::cout<<"Working on group number "<<count_gp<<":"<<std::endl;
		member_iter = group_iter->begin();
		count_mb = 0;

		while( member_iter != group_iter->end() )
		{
			Member memb = vec_memb[ *member_iter ];
			if( ! memb.isEmpty() )
			{
				count_mb++;
				count_all++;
				SDGString subseqname;
				std::string chr_name;
				if(verbose>1)
				{
					std::cout<<"\tsubseq from "<<memb.getNumChr()<<std::flush;
					std::cout<<" "<<memb.getStart()<<".."<<memb.getEnd()
					<<" length member:"<<memb.getLengthSet()<<std::endl;
				}
				if( memb.idlist.front() > 0 ) //subject
				{
					if( samedb )
						chr_name = match_map.numQ2name( memb.getNumChr() );
					else
						chr_name = match_map.numS2name( memb.getNumChr() );
					subseqname = "MbS"+SDGString(count_all)
					+"Gr"+SDGString(count_gp)
					+"Cl"+SDGString(gr2clust[count_gp]);
				}
				else  // query
				{
					chr_name = match_map.numQ2name( memb.getNumChr() );
					subseqname = "MbQ"+SDGString(count_all)
					+"Gr"+SDGString(count_gp)
					+"Cl"+SDGString(gr2clust[count_gp]);
				}
				if(verbose>1)
				{
					std::cout<<"\tseqname:"<<chr_name<<std::endl;
					std::cout<<"\tsubseqname:"<<subseqname<<std::endl;
				}
				membmap.add(RangeSeq(memb,subseqname,chr_name));
			}
			member_iter++;
		}
		count_gp++;
		group_iter++;
	}

	std::ostringstream filename;
	filename<<grouper_parameter->getPrefixFileName()<<".group.c"
	<<grouper_parameter->getCoverage();

//	if(verbose>0)
//		std::cout<<"Writing 'map' file..."<<std::flush;
//	SDGString mapfile = filename.str() + ".map";
//	membmap.save( mapfile );
//	if(verbose>0)
//		std::cout<<" done"<<std::endl;

	if(verbose>0)
		std::cout<<"Writing 'set' file..."<<std::flush;
	SDGString setfile = filename.str() + ".set";
	membmap.saveSet( setfile );
	if(verbose>0)
		std::cout<<" done"<<std::endl;

	if(verbose>0)
		std::cout<<"Writing 'fasta' file..."<<std::flush;
	SDGString seqfile = filename.str() + ".fa";
	membmap.writeSeq( seqfile, grouper_parameter->getBank() );
	if( grouper_parameter->getBank() != grouper_parameter->getQuery() )
		membmap.writeSeq( seqfile, grouper_parameter->getQuery() );
	if(verbose>0)
		std::cout<<" done"<<std::endl;

	if(verbose>0)
		std::cout<<"Writing 'param' file..."<<std::flush;
	SDGString parafile = filename.str() + ".param";
	std::ofstream parastream( parafile );
	grouper_parameter->view( parastream );
	parastream.close();
	if(verbose>0)
		std::cout<<" done"<<std::endl;

	if(verbose>0)
		std::cout<<"Writing 'txt' file..."<<std::flush;
	SDGString txtfile = filename.str() + ".txt";
	std::ofstream txtstream( txtfile );
	show_group( txtstream );
	txtstream.close();
	if(verbose>0)
		std::cout<<" done"<<std::endl;

	if( verbose > 0 )
		std::cout<<"Results were written."<<std::endl<<std::flush;
}


void BLRGroup::cluster( int verbose )
{
	if( verbose > 0 )
		std::cout<<"Cluster building (length filter="
		<<grouper_parameter->getGraphfilter()<<")..."
		<<std::endl<<std::flush;

	Graph<unsigned> graph;

	std::vector<unsigned> group_size;
	unsigned nb_group=0;
	for( GROUPLIST::iterator i=group_list.begin(); i!=group_list.end();
	i++, nb_group++ )
	{
		group_size.push_back( i->size() );
		nb_group++;
	};

	if(verbose>0)
		std::cout<<"Building graph..."<<std::flush;
	unsigned long count_gp1=0;
	unsigned long count_gp2=0;
	unsigned long count_mb1=0;

	GROUPLIST::iterator group_iter1 = group_list.begin();
	while( group_iter1 != group_list.end() )
	{
		std::map<unsigned long ,unsigned > count;
		std::list<unsigned>::iterator member_iter1=group_iter1->begin();
		count_gp1++;

		while( member_iter1 != group_iter1->end() )
		{
			Member memb1=vec_memb[*member_iter1];
			if(!memb1.isEmpty())
			{
				count_mb1++;
				std::vector<unsigned> member_id
				=memIdx->search_member(memb1.getNumChr(),memb1.getMin(), memb1.getMax());
				for( std::vector<unsigned>::iterator i=member_id.begin();
				  i!=member_id.end(); i++ )
				{
					Member memb2=vec_memb[*i];
					int cover=0;
					if( memb1.isIncluded( memb2 ) )
					{
						cover=(unsigned)
						memb1.overlap_length( memb2 );
					}
					else if( memb1.isContained( memb2 ) )
					{
						cover=(unsigned)
						memb1.overlap_length(memb2);
					}
					else if(// member1 and member2 overlap and member1 first
							memb1.overlap(memb2)
							&& memb1>memb2)
					{
						cover=(unsigned)
						memb1.overlap_length(memb2);
					}
					else if(// member1 and member2 overlap and member 2 first
							memb1.overlap(memb2)
							&& (memb1<memb2
									|| memb1==memb2))
					{
						cover=(unsigned)
						memb1.overlap_length(memb2);
					}
					if( cover>=grouper_parameter->getGraphfilter())
					{
						count_gp2=0;
						for(GROUPLIST::iterator gi=group_list.begin();
						gi!=group_list.end(); gi++)
						{
							count_gp2++;
							if(memIdx->get_groupit_member_idx(memb2.id)==gi)
								break;
						}
						count[count_gp2]++;
					}
				}
			}
			member_iter1++;
		}
		graph.add_node(count_gp1);
		for(std::map<unsigned long,unsigned>::iterator i=count.begin();i!=count.end();i++)
		{
			graph.add_edge(count_gp1,i->first);
		}
		group_iter1++;
	}
	if(verbose>0)
		std::cout<<" done (size="<<graph.size()<<")"<<std::endl;

	if(verbose>0)
		std::cout<<"Search connected components..."<<std::flush;
	graph.connexComp( clust );
	if(verbose>0)
		std::cout<<" done"<<std::endl;

	std::ostringstream clus_filename;
	clus_filename<<grouper_parameter->getPrefixFileName()<<".group.c"
	<<grouper_parameter->getCoverage()<<".cluster.dot";
	graph.toDot( clus_filename.str() );

	unsigned count=0;
	for( std::vector< std::vector<unsigned> >::iterator l=clust.begin(); l!=clust.end(); l++ )
	{
		count++;
		for( std::vector<unsigned>::iterator i=l->begin(); i!=l->end(); i++ )
			gr2clust[*i]=count;
	}
	if( verbose > 0 )
		std::cout<<"Clusters were built."<<std::endl<<std::flush;
}


/*-----------------------------------------------------------------------------
Remove groups in which less than 3 members are not included in others.
(Keep groups in which at least 3 members are not included)
-----------------------------------------------------------------------------*/
void BLRGroup::include_filter( int verbose )
{
	if( verbose > 0 )
	{
		std::cout<<"Filtering groups with included members..."<<std::endl<<std::flush;
		std::cout<<"before: "<<group_list.size()<<" groups"<<std::endl;
	}
	std::list<GROUPLIST::iterator> lGrpItToErase;
	unsigned countGroups = 1;

	for( GROUPLIST::iterator it_group_list=group_list.begin();
	it_group_list!=group_list.end();
	it_group_list++ )
	{
		if( verbose > 1 )
			std::cout<<"group "<<countGroups<<std::endl;
		unsigned nbUniqMembers = 0;
		for( std::list<unsigned>::iterator it_member=it_group_list->begin();
		it_member!=it_group_list->end();
		it_member++ )
		{
			Member mb = vec_memb[ *it_member ];
			if( mb.isEmpty() )
				continue;
			if( verbose > 2 )
			{
				mb.view();
				std::cout<<"included: "<<std::boolalpha<<mb.getIncluded()<<std::endl;
			}
			if( ! mb.getIncluded() && mb.getStart() != 0 && mb.getEnd() != 0 )
				++ nbUniqMembers;
		}
		if( verbose > 1 )
			std::cout<<nbUniqMembers<<" unique member(s)"<<std::endl;
		if( nbUniqMembers < grouper_parameter->getIncludeFilter() )
		{
			lGrpItToErase.push_back( it_group_list );
			if( verbose > 1 )
				std::cout<<"erase"<<std::endl;
		}
		++ countGroups;
	}

	for( std::list<GROUPLIST::iterator>::iterator i=lGrpItToErase.begin();
	i!=lGrpItToErase.end();
	i++ )
	{
		for( std::list<unsigned>::iterator m=(*i)->begin(); m!=(*i)->end(); m++ )
			memIdx->erase_member_idx(*m);
		group_list.erase(*i);
	}

	if( verbose > 0 )
	{
		std::cout<<"after: "<<group_list.size()<<" groups"<<std::endl;
		std::cout<<"Groups with included members were filtered."<<std::endl;
	}
}

void BLRGroup::include_filter_old( int verbose )
{
	if(verbose>0)
		std::cout<<"Including filter..."<<std::flush;
	GROUPLIST::iterator group_iter1=group_list.begin();
	std::list<GROUPLIST::iterator> gr2erase;

	while(group_iter1!=group_list.end())
	{
		std::list<unsigned>::iterator member_iter1=group_iter1->begin();
		// all over member of a group
		while(member_iter1!=group_iter1->end())
		{
			Member memb1=vec_memb[*member_iter1];
			std::vector<unsigned> member_id
			=memIdx->search_member(memb1.getNumChr(),memb1.getMin(), memb1.getMax());
			for(std::vector<unsigned>::iterator i=member_id.begin();
			i!=member_id.end(); i++)
			{
				if(*i!=*member_iter1)
				{
					Member memb2=vec_memb[*i];
					if(//member2 included in member1
							memb1.isIncluded(memb2))
					{
						if(std::find(gr2erase.begin(),gr2erase.end(),
								memIdx->get_groupit_member_idx(*i))
						==gr2erase.end())
							gr2erase.push_back(memIdx->get_groupit_member_idx(*i));
					}
				}
			}
			member_iter1++;
		}
		group_iter1++;
	}

	for(std::list<GROUPLIST::iterator>::iterator i=gr2erase.begin();
	i!=gr2erase.end();i++)
	{
		for(std::list<unsigned>::iterator m=(*i)->begin();
		m!=(*i)->end();m++)
			memIdx->erase_member_idx(*m);
		group_list.erase(*i);
	}
	if(verbose>0)
		std::cout<<" done"<<std::endl;
}


void BLRGroup::group_size_filter( int verbose )
{
	if(verbose>0)
	{
		std::cout<<"Filtering groups on size..."<<std::endl<<std::flush;
		std::cout<<"before: "<<group_list.size()<<" groups"<<std::endl;
	}
	GROUPLIST::iterator group_iter=group_list.begin();
	while( group_iter != group_list.end() )
	{
		if( group_iter->size() < grouper_parameter->getSizefilter() )
		{
			for( std::list<unsigned>::iterator m=(group_iter)->begin();
			m!=(group_iter)->end(); m++ )
				memIdx->erase_member_idx( *m );
			group_iter = group_list.erase( group_iter );
		}
		else
			group_iter++;
	}
	if(verbose>0)
	{
		std::cout<<"after: "<<group_list.size()<<" groups"<<std::endl;
		std::cout<<"Groups were filtered on size."<<std::endl;
	}
}


unsigned BLRGroup::getNbMembers( void )
{
	unsigned nbMembers = 0;
	for( GROUPLIST::iterator group_iter=group_list.begin();
			group_iter!=group_list.end(); group_iter++ )
		nbMembers += getNbMembers( group_iter );
	return nbMembers;
}


unsigned BLRGroup::getNbMembers( GROUPLIST::iterator group_iter )
{
	unsigned nbMembers = 0;
	for( std::list<unsigned>::iterator member_iter=group_iter->begin();
			member_iter!=group_iter->end(); ++ member_iter )
		if( vec_memb[ *member_iter ].getLengthSet() > 0 )
			++ nbMembers;
	return nbMembers;
}
