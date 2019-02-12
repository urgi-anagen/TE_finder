/**
 *
 * BLRBinSearch.h
 *
 */

#ifndef BLRBINSEARCH_H
#define BLRBINSEARCH_H

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include <cmath>
#include <mutex>

#include "RangeAlignSet.h"


class BLRBinSearch
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
  
  friend void threadedSearch(BLRBinSearch* bs, unsigned start, unsigned end, std::vector<unsigned>& vm);
  
  std::map<unsigned, std::list<coord> > mcoord;

  const static unsigned min_lvl=3;
  const static unsigned max_lvl=6;
  const static unsigned long max_coord=10000000000;

  unsigned long nb_bin_per_lvl;
  unsigned bin_size [8];

  unsigned getIdx(unsigned val, unsigned bin_lvl) const
  {
	  return unsigned((bin_lvl-min_lvl+1)*nb_bin_per_lvl)+getBin(val,bin_lvl);
  }
  unsigned getBin(unsigned val, unsigned bin_lvl) const // return the bin number on a given level where val reside
  {
	  return unsigned(val/bin_size[bin_lvl]);
  }

 public:
  BLRBinSearch(void): bin_size{1,10,100,1000,10000,100000,1000000,10000000}
    {
	  nb_bin_per_lvl=(max_coord/pow((float)10,(int)(max_lvl-min_lvl)));
    };

  virtual ~BLRBinSearch(void)
    {
    };

  std::vector<unsigned> search(unsigned start, unsigned end)
    {
       std::vector<unsigned> vm;
       for( unsigned bin_lvl=min_lvl; bin_lvl<=max_lvl; bin_lvl++)
       {
          	unsigned start_bin=getBin(start,bin_lvl);
            unsigned end_bin=getBin(end,bin_lvl);
            for(unsigned bin=start_bin; bin<=end_bin; bin++)
            {
            	unsigned idx=unsigned((bin_lvl-min_lvl+1)*nb_bin_per_lvl)+bin;
            	std::list<coord>& lcoord=mcoord[idx];
            	for(std::list<coord>::iterator i=lcoord.begin();i!=lcoord.end();i++)
            		if((i->start>=start && i->start<=end)
            				||(i->end>=start && i->end<=end)
            				||(i->start<=start && i->end>=end))
            			vm.push_back(i->id);
    	   }
       }
      return vm;
    };

    void insert(unsigned start, unsigned end, unsigned id)
    {
    	unsigned idx=0;
        for( unsigned bin_lvl=min_lvl; bin_lvl<max_lvl; bin_lvl++)
         {
			if(getBin(start,bin_lvl)==getBin(end,bin_lvl))
        	{
        		idx=getIdx(start,bin_lvl);
        		break;
        	}
        }
       if(idx==0)
       {
    	   idx=getIdx(start,max_lvl);
       }
       std::list<coord>& lcoord=mcoord[idx];
       lcoord.insert(lcoord.begin(),coord(start,end,id));
    };

    void erase(unsigned start, unsigned end, unsigned id)
    {
       	unsigned idx=0;
        for( unsigned bin_lvl=min_lvl; bin_lvl<max_lvl; bin_lvl++)
          {
			if(getBin(start,bin_lvl)==getBin(end,bin_lvl))
			{
				idx=getIdx(start,bin_lvl);
				break;
			}
		}
	   if(idx==0)
	   {
		   idx=getIdx(start,max_lvl);
	   }
	   std::list<coord>& lcoord=mcoord[idx];
	   lcoord.remove(coord(start,end,id));
    };
};

//threaded search
void threadedSearch(BLRBinSearch* bs, unsigned start, unsigned end, std::vector<unsigned>& vm); 
#endif



