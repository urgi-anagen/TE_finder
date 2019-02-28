#include <Range.h>
#include <iostream>


/*******************************************************
 *
 *    Constructeur de la classe
 *
 *******************************************************/

Range::Range(const Range *o)
{
  start=o->start;
  end=o->end;
}

Range::Range(const Range  &o)
{
  start=o.start;
  end=o.end;
}

Range::Range(ulong s, ulong e) 
  : start(s),end(e)
{
}


/**************************************************************
 * 
 * Destructeur
 *
 **************************************************************/

Range::~Range(void)
{}

void *Range::clone(void) const
{
  return (void*) new Range(this);
};



