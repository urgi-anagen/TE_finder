/**
 * \file BLRGroupThreads.h
 * \brief Header file for the class BLRGroup
 */

#ifndef BLRGROUPTHREADS_H
#define BLRGROUPTHREADS_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include <cmath>
#include <thread>

#include "Reference.h"
#include "SDGString.h"
#include "BlastMatch.h"
#include "BLRBioSeqDB.h"
#include "BLRMatchMap.h"
#include "BLRGrouperParameter.h"
#include "RangeSeq.h"
#include "RangeMap.h"
#include "BLRMemIdxBin.h"
//#include "threadpool/threadpool.h"

/**
 * \class BLRGroup
 * \brief Class implementing GROUPER data structure
 *
 */

class BLRGrouper;

class BLRGroup
{
	friend class BLRGrouper;

	private:
		BLRGrouperParameter grouper_parameter;
		BLRMatchMap match_map;

		GROUPLIST group_list;   // list of lists, each sub-list containing the identifier (unsigned integer) of the members in the group
		VMEMB vec_memb;         // vector containing all the members
		BLRMemIdx* memIdx;
		unsigned count_memb;

		BLRBioSeqDB bank_db;
		BLRBioSeqDB query_db;

		std::vector< std::vector<unsigned> > clust;
		std::map<unsigned,unsigned> gr2clust;

		//threadpool::ThreadPool pool;


	public:

		BLRGroup(BLRGrouperParameter gp, BLRMatchMap mm,unsigned nb_seq):
			grouper_parameter(gp), match_map(mm), group_list(0), vec_memb(0), count_memb(1)
			{
				Member null_mb;
				vec_memb.push_back(null_mb);
				memIdx=new BLRMemIdx(nb_seq);
			};

//		BLRGroup(const BLRGroup& g ) // not sure it works properly! to be tested?
//			{
//			grouper_parameter=g.grouper_parameter;
//			match_map=g.match_map;
//
//			group_list=g.group_list;
//			vec_memb=g.vec_memb;
//			memIdx=new BLRMemIdx(*(g.memIdx)); //doesn't work!
//			count_memb=g.count_memb;
//
//			bank_db=g.bank_db;
//			query_db=g.query_db;
//			clust=g.clust;
//			gr2clust=g.gr2clust;
//			};

		~BLRGroup( void )
		{
			delete memIdx;
		};

		void view_a_group( GROUPLIST::iterator& g )
		{
			std::cout<<"\ngroup";
			for(std::list<unsigned>::iterator member_iter=g->begin();
			member_iter!=g->end();member_iter++)
				std::cout<<":"<<*member_iter<<std::flush;
		}

		/** \fn void addMember( GROUPLIST::iterator&, Member& );
		 * \brief Add a new member in MEMBLIST
		 * \param gi: reference on current group
		 * \param mem: reference on member object
		 */
		std::list<unsigned>::iterator addMember( GROUPLIST::iterator, Member& );
		std::list<unsigned>::iterator eraseMember( GROUPLIST::iterator gi, std::list<unsigned>::iterator );

		/** \fn void addGroup( Member&, Member& );
		 * \brief Add a new group in group_list
		 * \param memQ: member created by query match
		 * \param memS: member created by subject match
		 */
		GROUPLIST::iterator addGroup( Member&, Member& );

		/** \fn void mergeGroup( GROUPLIST::iterator&, GROUPLIST::iterator& );
		 * \brief Merge two groups
		 * \param iQ: current group find by query match
		 * \param iS: current group find by subject match
		 */
		void mergeGroup( GROUPLIST::iterator&, GROUPLIST::iterator& );
		void mergeWithExtGroup( GROUPLIST::iterator iQ, GROUPLIST::iterator& iS, BLRGroup& ext_gr,
				std::list<unsigned>::iterator m_it1, unsigned m2 );
		void addExtGroup( BLRGroup& ext_gr );

		/** \fn
		 * \brief Merge two members
		 * \param ma reference on current alignment
		 * \param ga reference on Grouplist's current alignment
		 * \param ma reference on Listmember's current alignment
		 * \param gb reference on Grouplist's current member
		 * \param mb reference on Listmember's current member
		 * \param c reference on current Complement alignement
		 */
		void mergeMember( GROUPLIST::iterator& ga, unsigned ma,
				GROUPLIST::iterator& gb, unsigned mb,
				GROUPLIST::iterator& gc );

		/** \fn void mergeAlign( RangeAlignSet& align, long long id, unsigned member_id );
		 * \brief Merge a RangeAlignSet with an existing Member
		 * \param align RangeAlignSet instance to merge
		 * \param member_id member id
		 */
		void mergeAlign( RangeAlignSet& align, long long id, unsigned member_id );

		/** \fn void Inverse( GROUPLIST::iterator& );
		 * \brief Change strand of all members in MEMBLIST
		 * \param g group to inverse
		 */
		void Inverse( GROUPLIST::iterator& );

		/** \fn
		 * \brief Build group_list according to match
		 * \param ai current alignment
		 * \param iQ group found by Query
		 * \param iS group found by Subject
		 */
		void build_group( RangeAlignSet& alignQ, RangeAlignSet& alignS,
				long long id,
				GROUPLIST::iterator&, GROUPLIST::iterator&,
				int verbose=0 );

		void writePath(const SDGString& filename,std::list<RangePairSet> l, int verbose);

		Member& getRefMember(unsigned memb_num) { return vec_memb[memb_num];};

		void erase_member_idx(unsigned memb_num) {memIdx->erase_member_idx(memb_num);};
		GROUPLIST::iterator erase_group(GROUPLIST::iterator it) {return group_list.erase(it);};

		GROUPLIST::iterator begin(void) {return group_list.begin();};
		GROUPLIST::iterator end(void) {return group_list.end();};

		std::vector<unsigned> search_member(unsigned seqnum, unsigned start, unsigned end)
				{return memIdx->search_member(seqnum, start, end);};
		bool get_groupit_member_idx( unsigned member_id, GROUPLIST::iterator& gi)
		{
				return memIdx->get_groupit_member_idx( member_id, gi );
		};

		void show_group(  std::ostream& out=std::cout );

		void save( void );

		unsigned getNbGroups( void ) {return group_list.size();};

		/** \fn unsigned getNbMembers( void );
		 * \brief Get the number of members in group_list
		 */
		unsigned getNbMembers( void );

		/** \fn unsigned getNbMembers( GROUPLIST::iterator group_iter );
		 * \brief Get the number of members in group corresponding to group_iter
		 */
		unsigned getNbMembers( GROUPLIST::iterator group_iter );

		/** \fn void mergeMembersInGroups( int verbose=0 );
		 * \brief Merge members within each group
		 * \param verbose verbosity level
		 */
		void mergeMembersInGroups( int verbose=0 );
		
		std::list<unsigned>::iterator group_begin( GROUPLIST::iterator group_iter )
		{return group_iter->begin();};

		std::list<unsigned>::iterator group_end( GROUPLIST::iterator group_iter )
		{return group_iter->end();};

		GROUPLIST getGroupList(void) {return group_list;};

//		std::list<RangePairSet> getRpsListAfterLoad( int verbose ); not used?

		bool operator==(const BLRGroup& g) const
		{
			if(grouper_parameter!=g.grouper_parameter) return false;
			if(match_map!=g.match_map) return false;
			if(group_list!=g.group_list) return false;
			if(vec_memb!=g.vec_memb) return false;
			if(count_memb!=g.count_memb) return false;
			if(bank_db!=g.bank_db) return false;
			if(query_db!=g.query_db) return false;
			if(clust!=g.clust) return false;
			if(gr2clust!=g.gr2clust) return false;
			return true;
		};
};


#endif
