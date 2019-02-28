/**
 *
 * BLRVecBinSearch.h
 *
 */

#ifndef BLRVECBINSEARCH_H
#define BLRVECBINSEARCH_H

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <cmath>

#include "RangeAlignSet.h"



class BLRVecBinSearch
{
	const static unsigned nb_lvl=3;
	const static unsigned long max_coord=10000000; //10Mb

	unsigned bin_size [nb_lvl];

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

  class BLRBinVec
  {
	  std::vector< std::vector< std::list<coord> > > lvl;
	  unsigned bsize [nb_lvl];



  public:
	  BLRBinVec(void): lvl(nb_lvl),bsize{1000,10000,100000}
	  {
		  for(unsigned bin_lvl=0;bin_lvl<nb_lvl;bin_lvl++)
		  {
				 lvl[bin_lvl].resize((unsigned long)floor((unsigned long)((max_coord/bsize[bin_lvl])+0.5)));
		  }
	  };

	  virtual ~BLRBinVec(void){};

	  void insertCoord(unsigned bin_lvl, unsigned bin, const coord& c)
	  {
			  lvl[bin_lvl][bin].push_back(c);
	  };

	  void removeCoord(unsigned bin_lvl, unsigned bin, const coord& c)
	  {
			 lvl[bin_lvl][bin].remove(c);
	  };

	  const std::list<coord>& getCoordList(unsigned bin_lvl, unsigned bin)
	  {
		  return lvl[bin_lvl][bin];
	  };

  };


  BLRBinVec mcoord;

  unsigned getBin(unsigned val, unsigned bin_lvl)
  {
	  return unsigned(floor(val/bin_size[bin_lvl]));
  }

 public:
  BLRVecBinSearch(void): bin_size{1000,10000,100000}
    { };

  virtual ~BLRVecBinSearch(void)
    { };

  std::vector<unsigned> search(unsigned start, unsigned end)
    {
       std::vector<unsigned> vm;
       for( unsigned bin_lvl=0; bin_lvl<nb_lvl; bin_lvl++)
       {
          	unsigned start_bin=getBin(start,bin_lvl);
            unsigned end_bin=getBin(end,bin_lvl);
            for(unsigned bin=start_bin; bin<=end_bin; bin++)
            {
    	       const std::list<coord>& lcoord=mcoord.getCoordList(bin_lvl,bin);
    	       for(std::list<coord>::const_iterator i=lcoord.begin();i!=lcoord.end();i++)
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
        for( unsigned bin_lvl=0; bin_lvl<nb_lvl; bin_lvl++)
         {
			if(getBin(start,bin_lvl)==getBin(end,bin_lvl))
        	{
        		mcoord.insertCoord(bin_lvl,getBin(start,bin_lvl),coord(start,end,id));
        		return;
        	}
        }
        mcoord.insertCoord(nb_lvl-1,getBin(start,nb_lvl-1),coord(start,end,id));
    };

    void erase(unsigned start, unsigned end, unsigned id)
    {
        for( unsigned bin_lvl=0; bin_lvl<nb_lvl; bin_lvl++)
          {
			if(getBin(start,bin_lvl)==getBin(end,bin_lvl))
			{
				mcoord.removeCoord(bin_lvl,getBin(start,bin_lvl),coord(start,end,id));
				return;
			}
		}
        mcoord.removeCoord(nb_lvl-1,getBin(start,nb_lvl-1),coord(start,end,id));
    };
};
#endif



