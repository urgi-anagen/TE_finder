/*
 *
 * Diag.h
 *
 */

#ifndef DIAG_H
#define DIAG_H

#include <list>
#include <utility>
#include "Range.h"

class Diag : public std::list<Range>
{

  unsigned nwmatch;

 public:
  Diag(void):nwmatch(0){};

  void add(const Range& r)
    {
      if(back().getEnd()==r.getEnd()-1)
	{
	  back().setEnd(r.getEnd());
	}
      else
	{
	  if(back().getStart()==r.getStart()+1)
	    {
	      back().setStart(r.getStart());
	    }
	  else
	    {      	
	      push_back(r);
	    }
	}
      nwmatch++;
    }

  unsigned getLength(void){return  nwmatch;};
};

#endif
