/**
 *
 * BLRGroup.h
 *
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
#include "../../grouper.threads/BLRGrouperParameter.h"
#include "BLRBioSeqDB.h"
#include "BLRMatchMap.h"
#include "RangeSeq.h"
#include "RangeMap.h"
#include "Graph.h"
#include "BLRMemIdx.h"

class BLRGroup
{

 private:

  //! Member is sequence of blast match

  BLRGrouperParameter *grouper_parameter;
  double cover_limit;
  GROUPLIST group_list;
  VMEMB vec_memb;
  BLRMemIdx* memIdx;

  BLRMatchMap match_map;
  unsigned count_memb;

  BLRBioSeqDB bank_db;
  BLRBioSeqDB query_db;

  std::vector< std::vector<unsigned> > clust;
  std::map<unsigned,unsigned> gr2clust;

  /*! Function to add new member in MEMBLIST of LMEMBLIST iterator
    \param gi: reference on current group
    \param mem: reference on member object*/
  void addMember(GROUPLIST::iterator&,Member&);

  /*! Function to add new group in group_list
    \param memQ: member created by query match
    \param memS: member created by subject match*/
  void addGroup(Member&,Member&);

  /*! Function to merge two group
    \param iQ: current group find by query match
    \param iS: current group find by subject match*/
  void mergeGroup(GROUPLIST::iterator&,GROUPLIST::iterator&);

  /*! Function to merge two member
  -  a: reference on current alignment \param ga: on Grouplist \param ma: on Listmember
  -  b: reference on current memBer     \param gb: on Grouplist \param mb: on Listmember
  -  c: reference on current Complement alignement              \param gb: on Grouplist */
  void mergeMember(GROUPLIST::iterator& ga, unsigned ma,
		   GROUPLIST::iterator& gb, unsigned mb,
		   GROUPLIST::iterator& gc);

  /*--------------------------------------------------------
    - Function merge a RangeAlignSet to an existing Member
    - Arguments: align-> RangeAlignSet instance to merge
    -            member_id-> member id
    --------------------------------------------------------*/
  void mergeAlign( RangeAlignSet& align,long long id,unsigned member_id);

  /*! Function to change strand of all member in MEMBLIST
    \param g: group to inverse*/
  void Inverse(GROUPLIST::iterator&);

  /*! Function to compare query match to current member
    \param mi: current member
    \param Qi: group assigned to Q
    \param Si: group assigned to S
    \param Gi: current group
    \return booleen true if Qi and mi are same */
  bool compareQ(RangeAlignSet& alignQ,RangeAlignSet& alignS,
			long long id,
			unsigned mi,
			GROUPLIST::iterator& Qi,GROUPLIST::iterator& Si,
		GROUPLIST::iterator& Gi);
  /*! Function to compare subject match to current member
    \param mi: current member
    \param Qi: group assigned to Q
    \param Si: group assigned to S
    \param Gi: current group
    \return booleen true if Qi and mi are same */
  bool compareS(RangeAlignSet& alignQ,RangeAlignSet& alignS,
			long long id,
			unsigned mi,
			GROUPLIST::iterator& Qi,GROUPLIST::iterator& Si,
		GROUPLIST::iterator& Gi);

  /*! Function to build group_list according to match
    \param ai: current alignment
    \param iQ: group found by Query
    \param iS: group found by Subject*/
  void build_group(RangeAlignSet& alignQ,RangeAlignSet& alignS,
		   long long id,
		   GROUPLIST::iterator&,GROUPLIST::iterator&);

  void writePath(const SDGString& filename,std::list<RangePairSet> l, int verbose);

  public:

  BLRGroup(BLRGrouperParameter *para):
    grouper_parameter(para),cover_limit(0),
    group_list(0), vec_memb(0),memIdx(0),
    match_map(grouper_parameter), count_memb(0)
    {};
    ~BLRGroup(void)
      {
	delete memIdx;
      };


  //! Function to visit all MEMBLIST of group_list with this alignment
  void roam_group( RangePairSet& align, int verbose=0 );

  void group( int verbose=0 );

  void include_filter( int verbose=0 );

  void group_size_filter( int verbose=0 );

  //! Function to show all MEMBLIST of group_list
  void show_group(std::ostream& out=std::cout);

  /*! Function to save all sequence member in fasta format
    \param para: need to translate sequence and write in good file*/
  void save( int verbose=0 );

  void cluster( int verbose=0 );

  unsigned getNbMembers(void);

  GROUPLIST getGroupList(void){return group_list;};
  
  void view_group(GROUPLIST::iterator& g)
  {
    std::cout<<"\ngroup";
    for(std::list<unsigned>::iterator member_iter=g->begin();
	member_iter!=g->end();member_iter++)
      std::cout<<":"<<*member_iter<<std::flush;
  }

  VMEMB getVecMemb(void){return vec_memb;};
  
  std::list<RangePairSet> getRpsListAfterLoad( int verbose);	


}; // BLRGroup

#endif

