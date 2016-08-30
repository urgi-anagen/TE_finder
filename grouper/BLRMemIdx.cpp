/**
 *
 * BLRMemIdx.cpp
 *
 **/
#include "BLRMemIdx.h"

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
