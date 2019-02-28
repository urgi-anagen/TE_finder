/**
 *
 * BLRMemIdx.h
 *
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
#include "RangeAlignSet.h"


struct Member : public RangeAlignSet
{
  std::list<long long> idlist;      //!< id of match which adjust membre

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

  Member& operator=(const Member& m)
  {
    *this=RangeAlignSet::operator=(m);
    id=m.id;
    idlist=m.idlist;
    return *this;
  };

  virtual ~Member(void){};
  virtual void* clone(void) const
  {
    return (void*) new Member(this);
  };
};


typedef std::vector<Member> VMEMB;
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

  typedef BLRExhaustSearch SearchRange;

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

      SearchSeqRange(unsigned ntree): vsrch(0)
	{
	  resize(ntree);
	  for( std::vector<SearchRange*>::iterator i=vsrch.begin();
	       i!=vsrch.end();i++)
	    *i=new SearchRange();
	};

      void resize(unsigned ntree, int verbose=0)
      {
    	  if(verbose>0)
    		  std::cout<<"seq #:"<<ntree;
    	  vsrch.resize(ntree);
    	  if(verbose>0)
    		  std::cout<<" resized to:"<<vsrch.size()<<std::endl;
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
  };

  SearchSeqRange* srange;


  std::vector< Member_idx* > vec_idx;

 public:

  BLRMemIdx(unsigned nseq) : vec_idx(0)
    {
      srange=new SearchSeqRange(nseq);
    }
  ~BLRMemIdx(void)
    {
    	delete srange;
	for(std::vector<Member_idx*>::iterator i=vec_idx.begin();
	    i!=vec_idx.end();i++)
	  if(*i!=NULL)
	    delete *i;
    }

  std::vector<unsigned> search_member(unsigned seqnum, unsigned start, unsigned end)
    {
      return srange->search(seqnum,start,end);
    }

  void insert_member_idx(const Member& m,GROUPLIST::iterator& gi)
    {

      srange->insert(m.getNumChr(),m.getMin(),m.getMax(),m.id);

      if(vec_idx.size()<=m.id)
	vec_idx.resize(m.id+100);
      vec_idx[m.id]=new Member_idx(m,gi);
    }

  void erase_member_idx(unsigned id)
    {
      if(vec_idx[id]!=NULL)
	{
	  srange->erase(vec_idx[id]->seq,
		       vec_idx[id]->start,
		       vec_idx[id]->end,
		       vec_idx[id]->id_member);
	  delete vec_idx[id];
	  vec_idx[id]=NULL;
	}
    }

  void adjust_member_idx(unsigned id, unsigned start, unsigned end)
    {
      if(vec_idx[id]!=NULL)
	{
	  srange->erase(vec_idx[id]->seq,
		       vec_idx[id]->start,
		       vec_idx[id]->end,
		       vec_idx[id]->id_member);
	  vec_idx[id]->start=start;
	  vec_idx[id]->end=end;
	  srange->insert(vec_idx[id]->seq,
		       vec_idx[id]->start,
		       vec_idx[id]->end,
		       vec_idx[id]->id_member);
	}
    }

  unsigned get_start_member_idx(unsigned id)
  {return vec_idx[id]->start;}
  unsigned get_end_member_idx(unsigned id)
  {return vec_idx[id]->end;}
  unsigned get_seq_member_idx(unsigned id)
  {return vec_idx[id]->seq;}
   GROUPLIST::iterator get_groupit_member_idx(unsigned id)
  {return vec_idx[id]->it_group;}
   void set_groupit_member_idx(unsigned id,GROUPLIST::iterator& i)
  {vec_idx[id]->it_group=i;}

  //! overload operateur <<
  friend std::ostream& operator<<(std::ostream&,Member&);

  friend std::ostream& operator<<(std::ostream&,const BLRMemIdx::Member_idx&);

}; // BLRMemIdx

#endif
