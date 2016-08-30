/*
 *
 * Diag.h
 *
 */

#ifndef DIAGCLONE_H
#define DIAGCLONE_H

#include <list>
#include <utility>
#include "Range.h"

class DiagClone : public std::list<Range> 
	
{
  private:
  unsigned nwmatch;

  //std::list<Range> ranges;
  

  public:
 /*
  typedef std::list<Range>::iterator iterator;
  iterator begin() { return ranges.begin(); }
  iterator end() { return ranges.end(); }
*/

  DiagClone(void):nwmatch(0){};
  void add(const Range& r);
  unsigned getLength(void){return  nwmatch;};
  /*
  Range& back(void);
  void push_back(const Range& r);
*/
};

#endif
