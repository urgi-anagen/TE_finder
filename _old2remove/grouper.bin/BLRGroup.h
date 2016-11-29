/**
 * \file BLRGroup.h
 * \brief Header file for the class BLRGroup
 */

#ifndef BLRGROUP_H
#define BLRGROUP_H

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include <cmath>

#include "Reference.h"
#include "SDGString.h"
#include "BlastMatch.h"
#include "BLRGrouperParameter.h"
#include "BLRBioSeqDB.h"
#include "BLRMatchMap.h"
#include "RangeSeq.h"
#include "RangeMap.h"
#include "Graph.h"
#include "BLRMemIdx.h"


/**
 * \class BLRGroup
 * \brief Class implementing GROUPER
 *
 * The GROUPER program takes as input matches in the "align" format and sequences in the "fasta" format.
 * It can join matches into chains, cluster chains based on coverage, and filter groups.
 */
class BLRGroup
{

	private:

		BLRGrouperParameter *grouper_parameter;
		double cover_limit;
		GROUPLIST group_list;   // list of lists, each sub-list containing the identifier (unsigned integer) of the members in the group
		VMEMB vec_memb;         // vector containing all the members
		BLRMemIdx* memIdx;

		void view_group( GROUPLIST::iterator& g )
		{
			std::cout<<"\ngroup";
			for(std::list<unsigned>::iterator member_iter=g->begin();
			member_iter!=g->end();member_iter++)
				std::cout<<":"<<*member_iter<<std::flush;
		}

		BLRMatchMap match_map;
		unsigned count_memb;

		BLRBioSeqDB bank_db;
		BLRBioSeqDB query_db;

		std::vector< std::vector<unsigned> > clust;
		std::map<unsigned,unsigned> gr2clust;

		//ThreadPool tp;

		/** \fn void addMember( GROUPLIST::iterator&, Member& );
		 * \brief Add a new member in MEMBLIST
		 * \param gi: reference on current group
		 * \param mem: reference on member object
		 */
		void addMember( GROUPLIST::iterator&, Member& );

		/** \fn void addGroup( Member&, Member& );
		 * \brief Add a new group in group_list
		 * \param memQ: member created by query match
		 * \param memS: member created by subject match
		 */
		void addGroup( Member&, Member& );

		/** \fn void mergeGroup( GROUPLIST::iterator&, GROUPLIST::iterator& );
		 * \brief Merge two groups
		 * \param iQ: current group find by query match
		 * \param iS: current group find by subject match
		 */
		void mergeGroup( GROUPLIST::iterator&, GROUPLIST::iterator& );

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
		 * \brief Compare query match to current member
		 * \param mi current member
		 * \param Qi group assigned to Q
		 * \param Si group assigned to S
		 * \param Gi current group
		 * \param verbose verbosity level
		 * \return boolean true if Qi and mi are same
		 */
		bool compareQ( RangeAlignSet& alignQ, RangeAlignSet& alignS,
				long long id,
				unsigned mi,
				GROUPLIST::iterator& Qi, GROUPLIST::iterator& Si,
				GROUPLIST::iterator& Gi,
				int verbose=0 );
		/** \fn
		 * \brief Compare subject match to current member
		 * \param mi current member
		 * \param Qi group assigned to Q
		 * \param Si group assigned to S
		 * \param Gi current group
		 * \param verbose verbosity level
		 * \return boolean true if Si and mi are same
		 */
		bool compareS( RangeAlignSet& alignQ, RangeAlignSet& alignS,
				long long id,
				unsigned mi,
				GROUPLIST::iterator& Qi, GROUPLIST::iterator& Si,
				GROUPLIST::iterator& Gi,
				int verbose=0 );

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

		struct searchOverlapArgs //to pass arguments to pthread
		{
			BLRMemIdx* memIdx;
			RangeAlignSet ras;
			std::vector<unsigned>* p_vid;
			searchOverlapArgs(BLRMemIdx* m, RangeAlignSet& r, std::vector<unsigned>& v){
				memIdx=m;
				ras=r;
				p_vid=&v;
			}
		};

		static void searchOverlap(void *args);

  void writePath(const SDGString& filename,std::list<RangePairSet> l, int verbose);

	public:

		BLRGroup( BLRGrouperParameter *para ):
			grouper_parameter(para),cover_limit(0),
			group_list(0), vec_memb(0),memIdx(0),
			match_map(grouper_parameter), count_memb(0)
			{

			};

		~BLRGroup( void )
		{
			delete memIdx;
		};

		/** \fn void roam_group( RangePairSet& align, int verbose=0 );
		 * \brief Visit all MEMBLIST of group_list with this alignment
		 * \param align
		 * \param verbose verbosity level
		 */
		void roam_group( RangePairSet& align, int verbose=0 );

		/** \fn void group( int verbose=0 );
		 * \brief Load the alignments, connect them and define groups
		 * \param verbose verbosity level
		 */
		void group( int verbose=0 );

		/** \fn void include_filter( int verbose=0 );
		 * \brief Keep groups in which at least X members are not included
		 * \param verbose verbosity level
		 */
		void include_filter( int verbose=0 );
		void include_filter_old( int verbose=0 );

		/** \fn void group_size_filter( int verbose=0 );
		 * \brief Filter groups with less than a given number of members
		 * \param verbose verbosity level
		 */
		void group_size_filter( int verbose=0 );

		/** \fn void show_group( std::ostream& out=std::cout );
		 * \brief Show all MEMBLIST of group_list
		 * \param out stream for the output
		 */
		void show_group( std::ostream& out=std::cout );

		/** \fn void cluster( int verbose=0 );
		 * \brief Cluster similar groups together
		 * \param verbose verbosity level
		 */
		void cluster( int verbose=0 );

		/** \fn void save( int verbose=0 );
		 * \brief Save all sequence members in fasta format
		 * \param verbose verbosity level
		 */
		void save( int verbose=0 );

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
		
		std::list<RangePairSet> getRpsListAfterLoad( int verbose );

};

#endif
