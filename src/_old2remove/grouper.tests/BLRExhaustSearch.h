/**
 *
 * BLRExhaustSearch.h
 *
 */

#ifndef BLREXHAUSTSEARCH_H
#define BLREXHAUSTSEARCH_H

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include <cmath>

#include "RangeAlignSet.h"


class BLRExhaustSearch
{
  struct coord
  {
    unsigned start;
    unsigned end;
    unsigned id;

    coord(unsigned s, unsigned e, unsigned i) : start(s), end(e), id(i) {};

    friend bool operator==(const coord& c1, const coord& c2)
      {
	if(c1.start==c2.start && c1.end==c2.end && c1.id==c2.id)
	  return true;
	return false;
      };
    
  };
  std::list<coord> lcoord;

 public:
  BLRExhaustSearch(void)
    {
    };

  virtual ~BLRExhaustSearch(void)
    {
    };

  std::vector<unsigned> search(unsigned start, unsigned end)
    {
      std::vector<unsigned> vm;
      for(std::list<coord>::iterator i=lcoord.begin();
	  i!=lcoord.end();i++)
	if((i->start>=start && i->start<=end)
	   ||(i->end>=start && i->end<=end)
	   ||(i->start<=start && i->end>=end))
	  vm.push_back(i->id);
      return vm;
    };

    void insert(unsigned start, unsigned end, unsigned id)
    {
      lcoord.insert(lcoord.begin(),coord(start,end,id));
    };

    void erase(unsigned start, unsigned end, unsigned id)
    {
      lcoord.remove(coord(start,end,id));
    };
};


#endif



