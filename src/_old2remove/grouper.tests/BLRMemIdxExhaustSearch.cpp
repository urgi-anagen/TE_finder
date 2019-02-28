/*
 * \file BLRMemIdx.cpp
 */

#include "BLRMemIdxExhaustSearch.h"

std::ostream& operator<<(std::ostream& out, Member& mem)
{
  out<<mem.id<<"\t"; 
  out<<mem.getNumChr()<<"\t";
  out<<mem.getStart()<<"\t";
  out<<mem.getEnd()<<"\t";
  out<<mem.getRangeSetSize()<<"\t";
  out<<mem.getLengthSet()<<"\t";
  out<<mem.idlist.size()<<"\t";
  for(std::list<long long>::iterator i=mem.idlist.begin();
      i!=mem.idlist.end();i++)
    out<<*i<<" "; 
  return out;
};

std::ostream& operator<<(std::ostream& out, const BLRMemIdx::Member_idx& mem)
{
  out<<mem.id_member<<"\t"; 
  out<<mem.seq<<"\t";
  out<<mem.start<<"\t";
  out<<mem.end;
  return out;
};

void BLRMemIdx::insert_member_idx(const Member& m, GROUPLIST::iterator& gi)
{
	srange->insert( m.getNumChr(), m.getMin(), m.getMax(), m.id );
	if( vec_idx.size() <= m.id )
		vec_idx.resize( m.id+100 );
	vec_idx[ m.id ] = new Member_idx( m, gi );
}

void BLRMemIdx::erase_member_idx( unsigned id )
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

void BLRMemIdx::adjust_member_idx(unsigned id, unsigned start, unsigned end)
{
	if(vec_idx[id]!=NULL)
	{
		srange->erase( vec_idx[id]->seq,
				vec_idx[id]->start,
				vec_idx[id]->end,
				vec_idx[id]->id_member );
		vec_idx[id]->start = start;
		vec_idx[id]->end = end;
		srange->insert( vec_idx[id]->seq,
				vec_idx[id]->start,
				vec_idx[id]->end,
				vec_idx[id]->id_member );
	}
}
