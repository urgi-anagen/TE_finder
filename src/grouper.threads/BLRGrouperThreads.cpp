/**
 * \file BLRGroupThreads.cpp
 */

#include <cmath>
#include <map>

#include "BLRGrouperThreads.h"

//------------------------------------------------------------------------------------------------------------
bool BLRGrouper::compareQ( RangeAlignSet& alignQ, RangeAlignSet& alignS,
		long long id,
		unsigned mi,
		GROUPLIST::iterator& Qi, GROUPLIST::iterator& Si,
		GROUPLIST::iterator& Gi, BLRGroup& gr,
		int verbose )
{
	id=id*(-1);

	Member& mb = gr.getRefMember(mi);

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
			if( Qi==gr.end() && Si==gr.end() )
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
				gr.Inverse( Gi );
			}
		}
		gr.mergeAlign( alignQ, id, mi );
		return true;
	}
	return false;
}
//------------------------------------------------------------------------------------------------------------
bool BLRGrouper::compareS( RangeAlignSet& alignQ, RangeAlignSet& alignS,
		long long id,
		unsigned mi,
		GROUPLIST::iterator& Qi, GROUPLIST::iterator& Si,
		GROUPLIST::iterator& Gi, BLRGroup& gr,
		int verbose )
{
	Member& mb = gr.getRefMember(mi);

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
			if( Qi==gr.end() && Si==gr.end() )
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
				gr.Inverse( Gi );
			}
		}
		gr.mergeAlign( alignS, id, mi );
		return true;
	}
	return false;
}
//------------------------------------------------------------------------------------------------------------
void BLRGrouper::roam_group( RangePairSet& align, BLRGroup& gr, int verbose )
{
	long long id = align.getId();

	RangeAlignSet alignQ( align.getRangeAlignSetQ() );
	RangeAlignSet alignS( align.getRangeAlignSetS() );

	GROUPLIST::iterator group_iter = gr.end();
	GROUPLIST::iterator group_iterQ = gr.end();
	GROUPLIST::iterator group_iterS = gr.end();
	unsigned member_idQ=0, member_idS=0 ;


	if(verbose>0)
	{
		std::cout<<"query ("<<alignQ.getRangeSetSize()<<"): "
		<< alignQ.getStart()<<" "<<alignQ.getEnd()
		<<std::endl<<std::flush;

		std::cout<<"subject ("<<alignS.getRangeSetSize()<<"): "
		<< alignS.getStart()<<" "<<alignS.getEnd()
		<<std::endl<<std::flush;
	}
	if(verbose>1 && alignQ.getRangeSetSize()>1)
		alignQ.view();
	if(verbose>1 && alignS.getRangeSetSize()>1)
		alignS.view();

	//std::vector<unsigned> idQ,idS;

	/* debug thread
	std::cout<<"threads launched"<<std::endl;
	std::cout<<"memIdx ori:"<<memIdx<<std::endl;
	*/

    /* unthreaded implementation
	threadedIdxSearch(std::ref(memIdx),alignQ.getNumChr(),alignQ.getMin(),alignQ.getMax(),std::ref(idQ),verbose);
	threadedIdxSearch(std::ref(memIdx),alignQ.getNumChr(),alignQ.getMin(),alignQ.getMax(),std::ref(idQ),verbose);
     */

	/* thread implementation
	std::thread t1(threadedIdxSearch,std::ref(memIdx),alignQ.getNumChr(),alignQ.getMin(),alignQ.getMax(),std::ref(idQ),verbose);
    std::thread t2(threadedIdxSearch,std::ref(memIdx),alignS.getNumChr(),alignS.getMin(),alignS.getMax(),std::ref(idS),verbose);
    */
	/* threadpool implementation
	BLRMemIdx* local_memIdx=memIdx; // for pool.run lambda function)
	auto fQ = pool.run([&]() -> std::vector<unsigned>
			{ return threadedIdxSearch(local_memIdx,alignQ.getNumChr(),alignQ.getMin(),alignQ.getMax(),verbose); });
	auto fS = pool.run([&] () -> std::vector<unsigned>
			{ return threadedIdxSearch(local_memIdx,alignS.getNumChr(),alignS.getMin(),alignS.getMax(),verbose); });
	 */


    //Join the thread with the main thread
	//std::cout<<"waiting threads "<<std::endl;
	/* thread implementation
    t1.join();
	t2.join();
	*/
	/* threadpool implementation
	pool.wait();
	idQ=fQ.get();
	idS=fS.get();
	*/
	
	/* debug thread
	std::cout<<"threads joined"<<std::endl;
	std::cout<<"nb Q overlapping members: "<<idQ.size()<<std::endl;
	std::cout<<"nb S overlapping members: "<<idS.size()<<std::endl;
	*/

	// search on query
	std::vector<unsigned> idQ = gr.search_member( alignQ.getNumChr(),
			alignQ.getMin(),
			alignQ.getMax() );
	if(verbose>0)
		std::cout<<"nb overlapping members: "<<idQ.size()<<std::endl;
	// end search on query

	// search on subject
	std::vector<unsigned> idS=gr.search_member(alignS.getNumChr(),
				alignS.getMin(),
				alignS.getMax());
		if(verbose>0)
			std::cout<<"nb overlapping members: "<<idS.size()<<std::endl;
	// end search on subject

	for( std::vector<unsigned>::iterator i=idQ.begin();
	i!=idQ.end(); i++ )
	{
		unsigned member_id=*i;
		if (gr.get_groupit_member_idx( member_id, group_iter))
			if( compareQ( alignQ, alignS, id, member_id, group_iterQ, group_iterS,
					group_iter, gr, verbose-1 ) )
			{
				// the alignment query overlaps with a member
				if( group_iterQ == gr.end() )
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
					gr.mergeMember( group_iterQ, member_idQ, group_iter,
							member_id, group_iterS );
				}
			}
	}
	if(verbose>1)
		std::cout<<"included: "<<std::boolalpha<<alignQ.getIncluded()<<std::noboolalpha<<std::endl<<std::flush;

	for( std::vector<unsigned>::iterator i=idS.begin();
	i!=idS.end(); i++ )
	{
		unsigned member_id=*i;
		if (gr.get_groupit_member_idx( member_id, group_iter))
			if( compareS( alignQ, alignS, id, member_id, group_iterQ, group_iterS,
					group_iter, gr, verbose-1 ) )
			{
				// the alignment subject overlaps with a member
				if( group_iterS == gr.end() )
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
					gr.mergeMember( group_iterS, member_idS, group_iter,
							member_id, group_iterQ );
				}
			}
	}
	if(verbose>1)
		std::cout<<"included: "<<std::boolalpha<<alignS.getIncluded()<<std::noboolalpha<<std::endl<<std::flush;

	gr.build_group( alignQ, alignS, id, group_iterQ, group_iterS, verbose );
}
//------------------------------------------------------------------------------------------------------------
void BLRGrouper::group( std::list<RangePairSet>& rp_list,  BLRGroup& gr, int verbose )
{
	if( verbose > 0 )
		std::cout<<"Group building..."<<std::endl<<std::flush;

	unsigned long nbChains=0;  // a chain is a set of connected matches

	if(verbose>0)
	{
		std::cout<<"Roam groups (coverage limit: ";
		if( cover_limit <= 1 )
			std::cout<<std::floor(100*cover_limit)<<"%)..."<<std::endl<<std::flush;
		else
			std::cout<<cover_limit<<" bp)..."<<std::endl;
		std::cout<<"matches="<<rp_list.size()<<std::endl;
	}

	rp_list.sort( RangePair::greaterScoreIdLenCoord );
	
	std::list<RangePairSet>::iterator r=rp_list.begin();

	while( r!=rp_list.end() )
	{
		nbChains++;
		if(verbose>1)
			std::cout<<"chain "<<nbChains<<" (score="<<r->getScore()<<")"<<std::endl;
		r->setId( nbChains );

		roam_group( *r, gr, verbose-1 );

		rp_list.pop_front();
		r = rp_list.begin();

		if(verbose>2)
			gr.show_group();
	}

	if(verbose>0)
	{
		std::cout<<"number of chains: "<<nbChains<<std::endl;
		std::cout<<"number of members: "<<gr.getNbMembers()<<std::endl;
		std::cout<<"number of groups: "<<gr.getNbGroups()<<std::endl;
		std::cout<<"Roaming done."<<std::endl;
	}

	gr.mergeMembersInGroups( verbose );

	if(verbose>0)
	{
		std::cout<<"number of chains: "<<nbChains<<std::endl;
		std::cout<<"number of members: "<<gr.getNbMembers()<<std::endl;
		std::cout<<"number of groups: "<<gr.getNbGroups()<<std::endl;
	}

	if( nbChains==0 )
	{
		std::cout<<"Error when building groups: no alignment!"<<std::endl;
		exit( EXIT_FAILURE );
	}

	if( verbose > 0 )
		std::cout<<"Groups were built."<<std::endl<<std::flush;
}
//------------------------------------------------------------------------------------------------------------
void BLRGrouper::mergeGroupsLists(  BLRGroup* gr1_ptr, BLRGroup* gr2_ptr, int verbose )
{
	if( verbose > 0 )
		std::cout<<"Merge 2 group lists..."<<std::endl<<std::flush;

	std::list<GROUPLIST::iterator> found_gr_it;

	GROUPLIST::iterator group_iter1 = gr1_ptr->begin();
	while( group_iter1 != gr1_ptr->end() )
	{
		std::list<unsigned>::iterator member_iter1=group_iter1->begin();
		while( member_iter1 != group_iter1->end() )
		{
			Member memb1=gr1_ptr->getRefMember(*member_iter1);
			if(!memb1.isEmpty())
			{
				std::vector<unsigned> member_id
				=gr2_ptr->search_member(memb1.getNumChr(),memb1.getMin(), memb1.getMax());
				for( std::vector<unsigned>::iterator i=member_id.begin();
				  i!=member_id.end(); i++ )
				{
					Member memb2=gr2_ptr->getRefMember(*i);
					if(memb2.isEmpty()) continue;
					double cover=0;
                    if(memb1.overlap(memb2)) {
                        cover = (unsigned)
                                memb1.overlap_length(memb2);
                        if (cover/memb1.getLengthSet() >= cover_limit) {
                            GROUPLIST::iterator gr_it;
                            if (gr2_ptr->get_groupit_member_idx(*i, gr_it)) {
                                if (std::find(found_gr_it.begin(), found_gr_it.end(), gr_it) == found_gr_it.end()
                                    && gr_it != gr2_ptr->end()) {
                                    gr1_ptr->mergeWithExtGroup(group_iter1, gr_it, *gr2_ptr, member_iter1, *i);
                                    found_gr_it.push_back(gr_it);
                                }
                            }
                        }
                    }
				}
			}
			member_iter1++;
		}
		group_iter1++;
	}
	//add remaining groups
	if( verbose > 0 )
		std::cout<<"adding remaining groups"<<std::endl<<std::flush;
	gr1_ptr->addExtGroup(*gr2_ptr);
	if( verbose > 0 )
		std::cout<<"Group lists merged"<<std::endl<<std::flush;
}
//------------------------------------------------------------------------------------------------------------
void BLRGrouper::cluster( BLRGroup& gr, int verbose )
{
	if( verbose > 0 )
		std::cout<<"Cluster building (length filter="
		<<grouper_parameter.getGraphfilter()<<")..."
		<<std::endl<<std::flush;

	Graph<unsigned> graph;

	std::vector<unsigned> group_size;
	unsigned nb_group=0;
	for( GROUPLIST::iterator i=gr.begin(); i!=gr.end();
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

	GROUPLIST::iterator group_iter1 = gr.begin();
	while( group_iter1 != gr.end() )
	{
		std::map<unsigned long ,unsigned > count;
		std::list<unsigned>::iterator member_iter1=group_iter1->begin();
		count_gp1++;

		while( member_iter1 != group_iter1->end() )
		{
			Member memb1=gr.getRefMember(*member_iter1);
			if(!memb1.isEmpty())
			{
				count_mb1++;
				std::vector<unsigned> member_id
				=gr.search_member(memb1.getNumChr(),memb1.getMin(), memb1.getMax());
				for( std::vector<unsigned>::iterator i=member_id.begin();
				  i!=member_id.end(); i++ )
				{
					Member memb2=gr.getRefMember(*i);
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
					if( cover>=grouper_parameter.getGraphfilter())
					{
						count_gp2=0;
						for(GROUPLIST::iterator gi=gr.begin();
						gi!=gr.end(); gi++)
						{
							count_gp2++;
							GROUPLIST::iterator gm;
							if(gr.get_groupit_member_idx(memb2.id,gm) && gm==gi)
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
	graph.connexComp( gr.clust );
	if(verbose>0)
		std::cout<<" done"<<std::endl;

	std::ostringstream clus_filename;
	clus_filename<<grouper_parameter.getPrefixFileName()<<".group.c"
	<<grouper_parameter.getCoverage()<<".cluster.dot";
	graph.toDot( clus_filename.str() );

	unsigned count=0;
	for( std::vector< std::vector<unsigned> >::iterator l=gr.clust.begin(); l!=gr.clust.end(); l++ )
	{
		count++;
		for( std::vector<unsigned>::iterator i=l->begin(); i!=l->end(); i++ )
			gr.gr2clust[*i]=count;
	}
	if( verbose > 0 )
		std::cout<<"Clusters were built."<<std::endl<<std::flush;
}
/*-----------------------------------------------------------------------------
Remove groups in which less than 3 members are not included in others.
(Keep groups in which at least 3 members are not included)
-----------------------------------------------------------------------------*/
void BLRGrouper::include_filter( BLRGroup& gr, int verbose )
{
	if( verbose > 0 )
	{
		std::cout<<"Filtering groups with included members..."<<std::endl<<std::flush;
		std::cout<<"before: "<<gr.getNbGroups()<<" groups"<<std::endl;
	}
	std::list<GROUPLIST::iterator> lGrpItToErase;
	unsigned countGroups = 1;

	for( GROUPLIST::iterator it_group_list=gr.begin();
	it_group_list!=gr.end();
	it_group_list++ )
	{
		if( verbose > 1 )
			std::cout<<"group "<<countGroups<<std::endl;
		unsigned nbUniqMembers = 0;
		for( std::list<unsigned>::iterator it_member=it_group_list->begin();
		it_member!=it_group_list->end();
		it_member++ )
		{
			Member mb = gr.getRefMember( *it_member );
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
		if( nbUniqMembers < grouper_parameter.getIncludeFilter() )
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
			gr.erase_member_idx(*m);
		gr.erase_group(*i);
	}

	if( verbose > 0 )
	{
		std::cout<<"after: "<<gr.getNbGroups()<<" groups"<<std::endl;
		std::cout<<"Groups with included members were filtered."<<std::endl;
	}
}
//------------------------------------------------------------------------------------------------------------
void BLRGrouper::group_size_filter( BLRGroup& gr, int verbose )
{
	if(verbose>0)
	{
		std::cout<<"Filtering groups on size..."<<std::endl<<std::flush;
		std::cout<<"before: "<<gr.getNbGroups()<<" groups"<<std::endl;
	}
	GROUPLIST::iterator group_iter=gr.begin();
	while( group_iter != gr.end() )
	{
		if( group_iter->size() < grouper_parameter.getSizefilter() )
		{
			for( std::list<unsigned>::iterator m=(group_iter)->begin();
			m!=(group_iter)->end(); m++ )
				gr.erase_member_idx(*m);
			group_iter = gr.erase_group( group_iter );
		}
		else
			group_iter++;
	}
	if(verbose>0)
	{
		std::cout<<"after: "<<gr.getNbGroups()<<" groups"<<std::endl;
		std::cout<<"Groups were filtered on size."<<std::endl;
	}
}
