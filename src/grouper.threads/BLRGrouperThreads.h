/**
 * \file BLRGrouperThreads.h
 * \brief Header file for the class BLRGrouper
 */

#ifndef BLRGROUPERTHREADS_H
#define BLRGROUPERTHREADS_H

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
#include "BLRGrouperParameter.h"
#include "BLRMatchMap.h"
#include "RangeSeq.h"
#include "RangeMap.h"
#include "Graph.h"
#include "BLRGroupThreads.h"


/**
 * \class BLRGrouper
 * \brief Class implementing GROUPER algorithms
 *
 * The GROUPER program takes as input matches in the "align" format (or "path") and sequences in the "fasta" format.
 * It can join matches into chains, cluster chains based on coverage, and filter groups.
 */
class BLRGrouper
{

	private:

		BLRGrouperParameter grouper_parameter;
        double cover_limit ;

		bool compareQ( RangeAlignSet& alignQ, RangeAlignSet& alignS,
				long long id,
				unsigned mi,
				GROUPLIST::iterator& Qi, GROUPLIST::iterator& Si,
				GROUPLIST::iterator& Gi, BLRGroup& gr,
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
				GROUPLIST::iterator& Gi, BLRGroup& gr,
				int verbose=0 );

		/** \fn
		 * \brief Build group_list according to match
		 * \param ai current alignment
		 * \param iQ group found by Query
		 * \param iS group found by Subject
		 */

	public:

		BLRGrouper( const BLRGrouperParameter& para ):
			grouper_parameter(para)
			{cover_limit = grouper_parameter.getCoverage();};

		/** \fn void roam_group( RangePairSet& align, int verbose=0 );
		 * \brief Visit all MEMBLIST of group_list with this alignment
		 * \param align
		 * \param verbose verbosity level
		 */
		void roam_group( RangePairSet& align, BLRGroup& gr, int verbose=0 );

		/** \fn void group( int verbose=0 );
		 * \brief Load the alignments, connect them and define groups
		 * \param verbose verbosity level
		 */
		void group( std::list<RangePairSet>& rp_list, BLRGroup& gr, int verbose=0 );

		void mergeGroupsLists( BLRGroup* gr1_ptr, BLRGroup* gr2_ptr, int verbose );

		/** \fn void include_filter( int verbose=0 );
		 * \brief Keep groups in which at least X members are not included
		 * \param verbose verbosity level
		 */
		void include_filter( BLRGroup& gr, int verbose=0 );
		//void include_filter_old( int verbose=0 );

		/** \fn void group_size_filter( int verbose=0 );
		 * \brief Filter groups with less than a given number of members
		 * \param verbose verbosity level
		 */
		void group_size_filter(BLRGroup& gr, int verbose=0 );

		/** \fn void show_group( std::ostream& out=std::cout );
		 * \brief Show all MEMBLIST of group_list
		 * \param out stream for the output
		 */
		void cluster( BLRGroup& gr, int verbose=0 );

		/** \fn void save( int verbose=0 );
		 * \brief Save all sequence members in fasta format
		 * \param verbose verbosity level
		 */
};


#endif
