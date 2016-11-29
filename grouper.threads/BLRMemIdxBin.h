/*
 * \file BLRMemdIdxBin.h
 */

#ifndef BLRMEMIDXBIN_H
#define BLRMEMIDXBIN_H

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include <cmath>


#include "BLRBinSearch.h"


#include "RangeAlignSet.h"


struct Member : public RangeAlignSet
{
  // identifier of the chains composing the Member
  std::list<long long> idlist;

  // identifier of the Member, given in BLRGroup::addMember()
  unsigned id;

  Member(void):
    RangeAlignSet(),idlist(),id(0)
  {};

  Member(const Member& m):
    RangeAlignSet(m),idlist(m.idlist),id(m.id)
  {};

  Member(const RangeAlignSet& rs,long long l)
    :RangeAlignSet(rs),idlist(),id(0)
  {idlist.push_front(l);};

  ~Member(void){};

  void reset( void )
  {
	RangeAlignSet::reset();
	id = 0;
  	idlist.clear();
  };
  bool isEmpty(void)
  {
	  if(RangeAlignSet::isEmpty() || id==0) return true;
	  return false;
  }

  void view( void )
  {
	std::cout<<std::endl<<id<<"->"<<std::flush;
	std::cout<<" id list:";
	for(std::list<long long>::iterator i=idlist.begin(); i!=idlist.end();i++)
	{
		std::cout<<*i<<" ";
	}
	std::cout<<" -- ";
	RangeAlignSet::view();
  };
};


typedef std::vector<Member> VMEMB;

/*
 * \note each sub-list corresponds to a given group
 * \note it contains the indices of VMEMB of the Members in this group
 */
typedef std::list< std::list<unsigned> > GROUPLIST;


class BLRMemIdx
{
	struct Member_idx
	{
		unsigned id_member;
		GROUPLIST::iterator it_group;
		unsigned seq;
		unsigned start;
		unsigned end;
		Member_idx(const Member& m,GROUPLIST::iterator& gi)
		  :id_member(m.id),it_group(gi)
		{
		  start=m.getMin();
		  end=m.getMax();
		  seq=m.getNumChr();
		};
		Member_idx(unsigned numseq,unsigned s, unsigned e)
		  :id_member(0),it_group(), seq(numseq),start(s), end(e)
		{};

		friend bool operator<(const Member_idx m1, const Member_idx m2)
		{
			if(m1.start<m2.start)
				return true;
			else
				if(m1.start==m2.start && m1.end<m2.end)
					return true;
			return false;
		};
	};

  std::vector< Member_idx* > vec_idx;

  typedef BLRBinSearch SearchRange;

  class SearchSeqRange
  {
  	std::vector<SearchRange*> vsrch;

  public:
  	SearchSeqRange(void) : vsrch(0) {};

  	virtual ~SearchSeqRange(void)
  	{
  		for( std::vector<SearchRange*>::iterator i=vsrch.begin();
  				i!=vsrch.end();i++)
  			delete *i;
  	};

  	SearchSeqRange(unsigned nseq): vsrch(0)
		{
  		resize(nseq);
  		for( std::vector<SearchRange*>::iterator i=vsrch.begin();
  				i!=vsrch.end();i++)
  			*i=new SearchRange();
		};

  	void resize(unsigned nseq, int verbose=0)
  	{
  		if(verbose>0)
  			std::cout<<nseq<<" seq,"<<nseq;
  		vsrch.resize(nseq);
  		if(verbose>0)
  			std::cout<<" vsrch resized to "<<vsrch.size()<<std::endl;
  	};

  	std::vector<unsigned> search(unsigned seqnum, unsigned start, unsigned end)
			{
  		return vsrch[seqnum-1]->search(start,end);
			};

  	void insert(unsigned seqnum, unsigned start, unsigned end, unsigned id)
  	{
  		vsrch[seqnum-1]->insert(start,end,id);
  	};

  	void erase(unsigned seqnum, unsigned start, unsigned end, unsigned id)
  	{
  		vsrch[seqnum-1]->erase(start,end,id);
  	};

  	unsigned long size( void )
  	{
  		return vsrch.size();
  	}
  }; // SearchSeqRange

  SearchSeqRange* srange;


public:

  /*
   * \note a BLRMemIdx instance is initialized in BLRGroup::group() with the number of queries
   * \note it has two attributes: "srange" and "vec_idx"
   * \note it is filled in BLRGroup::addMember() when a new Member is created
   * \note srange is a ptr_SearchSeqRange which is itself a vector of ptr_SearchSeqRange
   * \note
   */
  BLRMemIdx(unsigned nseq) : vec_idx(0)
  {
  	srange=new SearchSeqRange(nseq);
  }

  ~BLRMemIdx(void)
  {
  	delete srange;
  	for(std::vector<Member_idx*>::iterator i=vec_idx.begin();i!=vec_idx.end();i++)
  		if(*i!=NULL)
  			delete *i;
  }

  friend std::ostream& operator<<(std::ostream&, Member&);

  friend std::ostream& operator<<(std::ostream&, const BLRMemIdx::Member_idx&);

  void insert_member_idx( const Member& m, GROUPLIST::iterator& gi );

  void erase_member_idx( unsigned id );

  void adjust_member_idx(unsigned id, unsigned start, unsigned end);

  std::vector<unsigned> search_member(unsigned seqnum, unsigned start, unsigned end)
  {
  	return srange->search(seqnum, start, end);
  }

  unsigned get_start_member_idx(unsigned id)
  {
	  if(vec_idx[id]==NULL) return 0;
	  return vec_idx[id]->start;
  }

  unsigned get_end_member_idx(unsigned id)
  {
	  if(vec_idx[id]==NULL) return 0;
	  return vec_idx[id]->end;
  }

  unsigned get_seq_member_idx(unsigned id)
  {
	  if(vec_idx[id]==NULL) return 0;
	  return vec_idx[id]->seq;
  }

  unsigned get_id_member_idx(unsigned id)
  {
	  if(vec_idx[id]==NULL) return 0;
	  return vec_idx[id]->id_member;
  }

  bool get_groupit_member_idx(unsigned id, GROUPLIST::iterator& it_group)
  {
	  if(vec_idx[id]==NULL) return false;
	  it_group=vec_idx[id]->it_group;
	  return true;
  }

  void set_groupit_member_idx(unsigned id, GROUPLIST::iterator& i)
  {
	  if(vec_idx[id]!=NULL)
		  vec_idx[id]->it_group=i;
  }

}; // BLRMemIdxBin

#endif
