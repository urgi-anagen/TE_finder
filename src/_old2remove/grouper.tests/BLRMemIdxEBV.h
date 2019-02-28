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

#include "BLRExhaustSearch.h"
#include "BLRBinSearch.h"
#include "BLRVecBinSearch.h"

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

  typedef BLRExhaustSearch SearchRangeE;
  typedef BLRBinSearch SearchRangeB;
  typedef BLRVecBinSearch SearchRangeV;

  class SearchSeqRangeE
  {
  	std::vector<SearchRangeE*> vsrch;

  public:
  	SearchSeqRangeE(void) : vsrch(0) {};
  	virtual ~SearchSeqRangeE(void)
  	{
  		for( std::vector<SearchRangeE*>::iterator i=vsrch.begin();
  				i!=vsrch.end();i++)
  			delete *i;
  	};

  	SearchSeqRangeE(unsigned nseq): vsrch(0)
		{
  		resize(nseq);
  		for( std::vector<SearchRangeE*>::iterator i=vsrch.begin();
  				i!=vsrch.end();i++)
  			*i=new SearchRangeE();
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
  }; // SearchSeqRangeE

  class SearchSeqRangeB
  {
  	std::vector<SearchRangeB*> vsrch;

  public:
  	SearchSeqRangeB(void) : vsrch(0) {};
  	virtual ~SearchSeqRangeB(void)
  	{
  		for( std::vector<SearchRangeB*>::iterator i=vsrch.begin();
  				i!=vsrch.end();i++)
  			delete *i;
  	};

  	SearchSeqRangeB(unsigned nseq): vsrch(0)
		{
  		resize(nseq);
  		for( std::vector<SearchRangeB*>::iterator i=vsrch.begin();
  				i!=vsrch.end();i++)
  			*i=new SearchRangeB();
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
  }; // SearchSeqRangeB

  class SearchSeqRangeV
  {
  	std::vector<SearchRangeV*> vsrch;

  public:
  	SearchSeqRangeV(void) : vsrch(0) {};
  	virtual ~SearchSeqRangeV(void)
  	{
  		for( std::vector<SearchRangeV*>::iterator i=vsrch.begin();
  				i!=vsrch.end();i++)
  			delete *i;
  	};

  	SearchSeqRangeV(unsigned nseq): vsrch(0)
		{
  		resize(nseq);
  		for( std::vector<SearchRangeV*>::iterator i=vsrch.begin();
  				i!=vsrch.end();i++)
  			*i=new SearchRangeV();
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
  }; // SearchSeqRangeV

  SearchSeqRangeE* srangeE;
  SearchSeqRangeB* srangeB;
  SearchSeqRangeV* srangeV;
  
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
  	srangeE=new SearchSeqRangeE(nseq);
  	srangeB=new SearchSeqRangeB(nseq);
  	srangeV=new SearchSeqRangeV(nseq);
  }

  ~BLRMemIdx(void)
  {
  	delete srangeE;
  	for(std::vector<Member_idx*>::iterator i=vec_idx.begin();i!=vec_idx.end();i++)
  		if(*i!=NULL)
  			delete *i;
  	delete srangeB;
  	for(std::vector<Member_idx*>::iterator i=vec_idx.begin();i!=vec_idx.end();i++)
  		if(*i!=NULL)
  			delete *i;
  	delete srangeV;
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
	std::vector<unsigned> vE=srangeE->search(seqnum, start, end);
	std::vector<unsigned> vB=srangeB->search(seqnum, start, end);
	std::vector<unsigned> vV=srangeV->search(seqnum, start, end);

	std::sort (vE.begin(),vE.end());
	std::sort (vB.begin(),vB.end());
	std::sort (vV.begin(),vV.end());
	std::vector<unsigned>::iterator ivB=vB.begin();
	for(std::vector<unsigned>::iterator ivE=vE.begin(); ivE!=vE.end(); ivE++)
	{
		if(*ivB!=*ivE)
		{
			std::cout<<"search differs with BinSearch:"<<*ivB<<" expected "<<*ivE<<std::endl;
			std::cout<<"ExSearch:"<<vec_idx[*ivE]->seq<<" "<<vec_idx[*ivE]->start<<" "<<vec_idx[*ivE]->end<<" "<<vec_idx[*ivE]->id_member<<std::endl;
			std::cout<<"BinSearch:"<<vec_idx[*ivB]->seq<<" "<<vec_idx[*ivB]->start<<" "<<vec_idx[*ivB]->end<<" "<<vec_idx[*ivB]->id_member<<std::endl;
			exit(EXIT_FAILURE);
		}
		ivB++;
	}

	std::vector<unsigned>::iterator ivV=vV.begin();
	for(std::vector<unsigned>::iterator ivE=vE.begin(); ivE!=vE.end(); ivE++)
	{
		if(*ivV!=*ivE)
		{
			std::cout<<"search differs with VecBinSearch:"<<*ivV<<" expected "<<*ivE<<std::endl;
			std::cout<<"ExSearch:"<<vec_idx[*ivE]->seq<<" "<<vec_idx[*ivE]->start<<" "<<vec_idx[*ivE]->end<<" "<<vec_idx[*ivE]->id_member<<std::endl;
			std::cout<<"VecSearch:"<<vec_idx[*ivV]->seq<<" "<<vec_idx[*ivV]->start<<" "<<vec_idx[*ivV]->end<<" "<<vec_idx[*ivV]->id_member<<std::endl;
			exit(EXIT_FAILURE);
		}
		ivV++;
	}
  	return vE;
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
