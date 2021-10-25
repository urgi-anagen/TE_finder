#include <RangeAlign.h>
#include <iostream>


/*******************************************************
 *
 *    Constructeur de la classe
 *
 *******************************************************/

RangeAlign::RangeAlign(const RangeAlign *o) : Range(o)
{
  num_chr=o->num_chr;
  name_seq=o->name_seq;
}

RangeAlign::RangeAlign(const RangeAlign &o) : Range(o)
{
  num_chr=o.num_chr;
  name_seq=o.name_seq;
}

RangeAlign::RangeAlign(long c, ulong s, ulong e) 
  : Range(s,e),num_chr(c),name_seq("")
{
}
RangeAlign::RangeAlign(std::string name, long c, ulong s, ulong e) 
  : Range(s,e),num_chr(c),name_seq(name)
{
}

/**************************************************************
 * 
 * Destructeur
 *
 **************************************************************/

RangeAlign::~RangeAlign(void)
{}

void *RangeAlign::clone(void) const
{
  return (void*) new RangeAlign(this);
};



