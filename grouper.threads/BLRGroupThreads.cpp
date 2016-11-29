/**
 * \file BLRGroupThreads.cpp
 */

#include <cmath>
#include <map>


#include "BLRGroupThreads.h"

//------------------------------------------------------------------------------------------------------------
std::list<unsigned>::iterator BLRGroup::addMember( GROUPLIST::iterator gi, Member& mem )
{
	mem.id = count_memb++;
	vec_memb.push_back( mem );
	memIdx->insert_member_idx( mem, gi );
	return 	gi->insert( gi->begin(), mem.id );
}
//------------------------------------------------------------------------------------------------------------
std::list<unsigned>::iterator BLRGroup::eraseMember( GROUPLIST::iterator gi, std::list<unsigned>::iterator i )
{
	if(vec_memb[*i].isEmpty()) return gi->end();
	vec_memb[*i].reset();
	memIdx->erase_member_idx( *i );
	return gi->erase(i);
}
//------------------------------------------------------------------------------------------------------------
GROUPLIST::iterator BLRGroup::addGroup( Member& memQ, Member& memS )
{
	std::list<unsigned> g;
	GROUPLIST::iterator group_it=group_list.insert( group_list.begin(), g );
	addMember( group_it, memQ );
	addMember( group_it, memS );
	return group_it;
}
//------------------------------------------------------------------------------------------------------------
void BLRGroup::mergeWithExtGroup( GROUPLIST::iterator iQ, GROUPLIST::iterator& iS, BLRGroup& ext_gr,
		std::list<unsigned>::iterator m_it1, unsigned m2 )
{

	vec_memb[*m_it1].merge(ext_gr.vec_memb[m2]);
	vec_memb[*m_it1].idlist.splice( vec_memb[*m_it1].idlist.end(), ext_gr.vec_memb[m2].idlist );
	if( vec_memb[*m_it1].getMin() < memIdx->get_start_member_idx( *m_it1 )
			|| vec_memb[*m_it1].getMax() > memIdx->get_end_member_idx( *m_it1 ) )
		memIdx->adjust_member_idx( *m_it1, vec_memb[*m_it1].getMin(), vec_memb[*m_it1].getMax() );
	ext_gr.eraseMember(iS,std::find(iS->begin(),iS->end(),m2));

	std::list<unsigned>::iterator i=iS->begin();
	while(i!=iS->end())
	{
		addMember(iQ,ext_gr.vec_memb[*i]);
		i=ext_gr.eraseMember(iS,i);
	}
	iS=ext_gr.group_list.erase( iS );
}

//------------------------------------------------------------------------------------------------------------
void BLRGroup::addExtGroup( BLRGroup& ext_gr )
{
	GROUPLIST::iterator i=ext_gr.group_list.begin();
	while(i!=ext_gr.group_list.end())
	{
		std::list<unsigned> lm;
		group_list.push_back(lm);
		GROUPLIST::iterator lm_it=group_list.end();
		lm_it--;
		std::list<unsigned>::iterator m=i->begin();
		while(m!=i->end())
		{
			addMember(lm_it,ext_gr.vec_memb[*m]);
			m=ext_gr.eraseMember(i,m);
		}
		i=ext_gr.group_list.erase( i );
	}
}
//------------------------------------------------------------------------------------------------------------
void BLRGroup::mergeGroup( GROUPLIST::iterator& iQ, GROUPLIST::iterator& iS )
{
	for( std::list<unsigned>::iterator i=iS->begin(); i!=iS->end(); i++ )
		memIdx->set_groupit_member_idx( *i, iQ );
	iQ->splice( iQ->begin(), *iS, iS->begin(), iS->end() );
	iS=group_list.erase( iS );
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
		if( vec_memb[ma].getMin() < memIdx->get_start_member_idx( ma )
				|| vec_memb[ma].getMax() > memIdx->get_end_member_idx( ma ) )
			memIdx->adjust_member_idx( ma, vec_memb[ma].getMin(), vec_memb[ma].getMax() );
		memIdx->erase_member_idx( mb );
		ga->remove( mb );
	}
}
//------------------------------------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------------------------------------
void BLRGroup::Inverse( GROUPLIST::iterator& g )
{
	for( std::list<unsigned>::iterator member_iter=g->begin();
	member_iter!=g->end();
	member_iter++ )
		vec_memb[ *member_iter ].reverse();
}
//------------------------------------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------------------------------------
void BLRGroup::mergeMembersInGroups( int verbose )
{
	if( verbose > 0 )
	std::cout<<"Merge members in groups..."<<std::endl;

	GROUPLIST::iterator it_group = group_list.begin();
	while( it_group != group_list.end() )
	{
		if( it_group->size() < 2 )
		{
			it_group++;
			continue;
		}

		std::list<unsigned>::iterator it_idPrevMb = it_group->begin();
		std::list<unsigned>::iterator it_beforeLastMb = it_group->end(); it_beforeLastMb--;
		while(true)
		{
			while(it_idPrevMb != it_beforeLastMb && vec_memb[ *it_idPrevMb ].isEmpty()) it_idPrevMb++;
			if( it_idPrevMb != it_beforeLastMb ) break;
			Member* pt_prevMb = &vec_memb[ *it_idPrevMb ];
			std::list<unsigned>::iterator it_idCurrMb = it_idPrevMb;
			++ it_idCurrMb;
			while( it_idCurrMb != it_group->end() )
			{
				Member* pt_currMb = &vec_memb[ *it_idCurrMb ];
				unsigned long minLength = std::min( pt_prevMb->getLength(), pt_currMb->getLength() );
				while( it_idCurrMb != it_group->end() && ! pt_currMb->isEmpty() && pt_prevMb != pt_currMb
						&& pt_prevMb->overlap_length( *pt_currMb ) >= 0.95 * minLength )
				{
					bool isIncludedPrevMb = pt_prevMb->getIncluded();
					bool isIncludedCurrMb = pt_currMb->getIncluded();
					unsigned lengthPrevMb = pt_prevMb->getLengthSet();
					unsigned lengthCurrMb = pt_currMb->getLengthSet();
					pt_prevMb->merge( *pt_currMb );
					pt_prevMb->idlist.splice( pt_prevMb->idlist.end(), pt_currMb->idlist );
					if( ! isIncludedPrevMb && isIncludedCurrMb && lengthPrevMb <= lengthCurrMb )
						pt_prevMb->setIncluded( true );
					it_idCurrMb=it_group->erase(it_idCurrMb);
					pt_currMb = &vec_memb[ *it_idCurrMb ];
				}
				if(it_idCurrMb != it_group->end() ) break;
				++ it_idCurrMb;
			}
			++ it_idPrevMb;
		}
		it_group++;
	}

	if( verbose > 0 )
		std::cout<<"Members in groups were merged."<<std::endl;
}
//------------------------------------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------------------------------------
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
						chr_name = match_map->numQ2name( memb.getNumChr() );
					else
						chr_name = match_map->numS2name( memb.getNumChr() );
				else  // query
					chr_name = match_map->numQ2name( memb.getNumChr() );
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
//------------------------------------------------------------------------------------------------------------
void BLRGroup::save( void )
{
	unsigned verbose=grouper_parameter->getVerbose();

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
						chr_name = match_map->numQ2name( memb.getNumChr() );
					else
						chr_name = match_map->numS2name( memb.getNumChr() );
					subseqname = "MbS"+SDGString(count_all)
					+"Gr"+SDGString(count_gp)
					+"Cl"+SDGString(gr2clust[count_gp]);
				}
				else  // query
				{
					chr_name = match_map->numQ2name( memb.getNumChr() );
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
	membmap.writeSeq( seqfile, match_map->getRefQueryDB() );
	if( !samedb && grouper_parameter->getQuery()!="<not set>")
		membmap.writeSeq( seqfile, match_map->getRefSubjectDB() );
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
//------------------------------------------------------------------------------------------------------------
unsigned BLRGroup::getNbMembers( void )
{
	unsigned nbMembers = 0;
	for( GROUPLIST::iterator group_iter=group_list.begin();
			group_iter!=group_list.end(); group_iter++ )
		nbMembers += getNbMembers( group_iter );
	return nbMembers;
}
//------------------------------------------------------------------------------------------------------------
unsigned BLRGroup::getNbMembers( GROUPLIST::iterator group_iter )
{
	unsigned nbMembers = 0;
	for( std::list<unsigned>::iterator member_iter=group_iter->begin();
			member_iter!=group_iter->end(); ++ member_iter )
		if( vec_memb[ *member_iter ].getLengthSet() > 0 )
			++ nbMembers;
	return nbMembers;
}
