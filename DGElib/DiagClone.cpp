#include "DiagClone.h"

void DiagClone::add(const Range& r)
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
