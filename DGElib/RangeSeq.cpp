#include <RangeSeq.h>
#include <iostream>


/*******************************************************
 *
 *    Constructeur de la classe
 *
 *******************************************************/

RangeSeq::RangeSeq(const RangeSeq *o) : Range(o)
{
  chr=o->chr;
  name=o->name;
  range_set=o->range_set;
}

RangeSeq::RangeSeq(const RangeSeq &o) : Range(o)
{
  chr=o.chr;
  name=o.name;
  range_set=o.range_set;
};

RangeSeq::RangeSeq(const std::string& n, const std::string& c, ulong s, ulong e) 
  : Range(s,e),name(n),chr(c) {}

RangeSeq::RangeSeq(const RangeAlignSet& r,
	     const std::string& n, const std::string& c)
{set(r,n,c);};

/**************************************************************
 * 
 * Destructeur
 *
 **************************************************************/

RangeSeq::~RangeSeq()
{
}

void *RangeSeq::clone() const
{
  return (void*) new RangeSeq(this);
};



