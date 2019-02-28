/*
 * \file BLRMemdIdx.h
 */

#ifndef BLRMEMIDX_H
#define BLRMEMIDX_H

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include <cmath>
#include <mutex>

//#include "BLRExhaustSearch.h"
#include "BLRBinSearch.h"
//#include "BLRVecBinSearch.h"

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

  Member(const RangeAlignSet& rs)
    :RangeAlignSet(rs),idlist(),id(0)
  {};

  Member(const RangeAlignSet* rs)
    :RangeAlignSet(*rs),idlist(),id(0)
  {};

  virtual ~Member(void){};

  Member& operator=(const Member& m)
  {
    *this=RangeAlignSet::operator=(m);
    id=m.id;
    idlist=m.idlist;
    return *this;
  };

  virtual void* clone(void) const
  {
    return (void*) new Member(this);
  };

  void reset( void )
  {
  	id = 0;
  	idlist.clear();
  	RangeAlignSet::reset();
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

  //typedef BLRExhaustSearch SearchRange;
  typedef BLRBinSearch SearchRange;
  //typedef BLRVecBinSearch SearchRange;

  class SearchSeqRange
  {
  	std::vector<SearchRange*> vsrch;

  public:
  
    friend std::vector<unsigned> threadedIdxSearch(BLRMemIdx* memIdx, unsigned seqnum, unsigned start, unsigned end, int verbose);
 
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
  
  friend std::vector<unsigned> threadedIdxSearch(BLRMemIdx* memIdx, unsigned seqnum, unsigned start, unsigned end, int verbose);
  
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
  	return vec_idx[id]->start;
  }

  unsigned get_end_member_idx(unsigned id)
  {
  	return vec_idx[id]->end;
  }

  unsigned get_seq_member_idx(unsigned id)
  {
  	return vec_idx[id]->seq;
  }

  unsigned get_id_member_idx(unsigned id)
  {
  	return vec_idx[id]->id_member;
  }

  GROUPLIST::iterator get_groupit_member_idx(unsigned id)
  {
  	return vec_idx[id]->it_group;
  }

  void set_groupit_member_idx(unsigned id, GROUPLIST::iterator& i)
  {
  	vec_idx[id]->it_group=i;
  }

}; // BLRMemIdx

#endif
