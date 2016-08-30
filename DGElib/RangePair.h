/**
 * \file RangePair.h
 * \brief Header file for the class RangePair
 */

#ifndef RANGEPAIR_H
#define RANGEPAIR_H

#include <iostream>
#include <SDGString.h>
#include <utility>
#include <cmath>
#include "BlastMatch.h"
#include "RangeAlign.h"
#include "SDGString.h"

/**
 * \class RangePair
 * \brief Implement a match as two fragments with E-value, score and identity
 *
 * This class has a pair of RangeAlign objects as well as other attributes.
 */
class RangePair: public std::pair<RangeAlign,RangeAlign>
{
  struct Less
  {
    Less(void){};
    int operator () (const RangePair& a, const RangePair& b) const
    {
      if(a.first<b.first) return true;
      if(a.first==b.first && a.second<b.second) return true;
      return false;
    };
  };

  struct Greater
  {
    Greater(void){};
    int operator () (const RangePair& a, const RangePair& b) const
    {
      if(a.first>b.first) return true;
      if(a.first==b.first && a.second>b.second) return true;
      return false;
    };
  };

  struct GreaterScore
  {
    GreaterScore(void){};
    int operator () (const RangePair& a, const RangePair& b) const
    {
      if(a.score>b.score)
	return true;
      if(a.score==b.score && a.length>b.length)
	return true;
      return false;
    };
  };

  struct GreaterScoreIdLenCoord
  {
	GreaterScoreIdLenCoord(void){};
    int operator () (const RangePair& a, const RangePair& b) const
    {
      if(a.score>b.score) return true;
      if (a.score<b.score) return false;
      if(a.identity>b.identity) return true;
      if(a.identity<b.identity) return false;
      if(a.length>b.length) return true;
      if(a.length<b.length) return false;
      if(a.first>b.first) return true;
      if(a.first<b.first) return false;
      return false;
    };
  };

  struct LessIdentity
  {
    LessIdentity (void){};
    int operator () (const RangePair& a, const RangePair& b) const
    {
      if(a.identity<b.identity)
	return true;
      if(a.identity==b.identity && a.length>b.length)
	return true;
      return false;
    };
  };

  struct GreaterLengthQ
  {
    GreaterLengthQ(void){};
    int operator () (const RangePair& a, const RangePair& b) const
    {
      if(a.first.getLength()>b.first.getLength())
	return true;
      else if (a.first.getLength()==b.first.getLength()
	    && (a.second.getLength()>b.second.getLength()) )
	return true;
      return false;
    };
  };

  struct GreaterLengthIdent
  {
    GreaterLengthIdent(void){};
    int operator () (const RangePair& a, const RangePair& b) const
    {
      if(a.length>b.length)
	return true;
      else
	{
	  if (a.length==b.length
	      && (a.length>b.length) )
	    return true;
	  else
	    if ((a.length==b.length)
		&& (a.length==b.length)
		&&(a.identity>b.identity))
	      return true;
	}
      return false;
    };
  };

  struct StrictLess
  {
    StrictLess(void){};
    int operator () (const RangePair& a, const RangePair& b) const
    {
     if(a.first.getMax()<b.first.getMin()) return true;
      return false;
    };
  };

 private:
  inline bool isQueryStartHasMoved(unsigned queryStart)
	{
	return queryStart != first.getStart();
	}

  inline bool isQueryEndHasMoved(unsigned queryEnd)
	{
	return queryEnd !=first.getEnd();
	}

 protected:
  double e_value;
  double identity;
  unsigned length;
  long score;
  unsigned long id;

 public:
  static const Less less;
  static const Greater greater;
  static const GreaterScore greaterScore;
  static const GreaterScoreIdLenCoord greaterScoreIdLenCoord;
  static const LessIdentity lessIdentity;
  static const GreaterLengthQ greaterLengthQ;
  static const GreaterLengthIdent greaterLengthIdent;
  static const StrictLess strictLess;

  RangePair(void): e_value(0), identity(0), length(0), score(0),id(0){};
  RangePair(const RangeAlign& r1,const RangeAlign& r2)
    :std::pair<RangeAlign,RangeAlign>(r1,r2),e_value(0), identity(0), length(0), score(0),id(0)
    {};
  RangePair( BlastMatch al );
  RangePair( SDGString line );


  friend bool operator==( const RangePair &rp1, const RangePair &rp2 );

  void clear(void)
    {
      first.set();second.set();
      e_value=-1.0;identity=0.0;
      length=score=id=0;
    };

  bool empty(void) const
    {return first.empty() || second.empty();};

  void setStrand(void)
    {
      unsigned minq=first.getMin();
      unsigned maxq=first.getMax();
      unsigned mins=second.getMin();
      unsigned maxs=second.getMax();
      if(
	 (first.isPlusStrand() && second.isPlusStrand())
	 ||(!first.isPlusStrand() && !second.isPlusStrand()))
	{
	  first.setStart(minq);
	  first.setEnd(maxq);
	  second.setStart(mins);
	  second.setEnd(maxs);
	}
      else
	{
	  first.setStart(minq);
	  first.setEnd(maxq);
	  second.setStart(maxs);
	  second.setEnd(mins);
	}
    }

  bool isPlusStrand(void) const
    { return second.isPlusStrand();};

  BlastMatch toAlign(void);

  RangeAlign& getRangeQ(void){return first;};
  RangeAlign& getRangeS(void){return second;}
  const RangeAlign& getRangeQ(void) const {return first;};
  const RangeAlign& getRangeS(void) const {return second;}
  void setRangeQ(const RangeAlign& r){first=r;};
  void setRangeS(const RangeAlign& r){second=r;};

  void set( SDGString );

  unsigned long getId(void){return id;};
  unsigned long getId(void) const {return id;};
  void setId(unsigned long i){id=i;};

  double getE_value(void){return e_value;};
  double getE_value(void) const {return e_value;};
  void setE_value(double p){e_value=p;};

  long getScore(void){return score;};
  long getScore(void) const {return score;};
  void setScore(long sc){score=sc;};

  double getIdentity(void){return identity;};
  double getIdentity(void) const {return identity;};
  void setIdentity(double i){identity=i;};

  unsigned getLength(void){return length;};
  unsigned getLength(void) const {return length;};
  void setLength(unsigned l){length=l;};

  void view(void)
    {
      std::cout<<(*this)<<std::endl;
    };
  void view(void) const
    {
      std::cout<<(*this)<<std::endl;
    };
  void viewWithLabel(void)
  {
     std::cout<<"rangeQ "<<std::endl;
     first.view(); 
     std::cout<<"rangeS "<<std::endl;
     second.view(); 
     std::cout<<"id "<<id<<" e_value "<<e_value<<" identity "<<identity<<" length "<<length<<" score "<<score<<std::endl;
  };
  void writetxt(std::ostream& out);

  void readReputer(std::istream& in);
  void readtxt(std::istream& in);
  void readlst(std::istream& in);

  void merge(RangePair& r);
  void merge(const RangePair& r);

  bool overlap(const RangePair& r) const
    {
      return first.overlap(r.first) && second.overlap(r.second);
    };

  bool overlap(const RangePair& r)
    {
      return first.overlap(r.first) && second.overlap(r.second);
    };

  bool overlap1(const RangePair& r) const
    {
      if(first.overlap(r.first)
	 || first.overlap(r.second)
	 || second.overlap(r.first)
	 || second.overlap(r.second))
	return true;
      else
	return false;
    };

  bool overlap1(const RangePair& r)
    {
      if(first.overlap(r.first)
	 || first.overlap(r.second)
	 || second.overlap(r.first)
	 || second.overlap(r.second))
	return true;
      else
	return false;
    };

  bool overlapQ(const RangePair& r) const
    {
      return first.overlap(r.first);
    };

  RangePair diffQ(const RangePair& r);

  friend std::ostream& operator<<(std::ostream& out, const RangePair& r)
    {
	  out<<r.first<<"\t";
	  out<<r.second<<"\t";
	  out<<"score="<<r.score;
	  out<<" e_value="<<r.e_value;
	  out<<" identity="<<r.identity;
	  out<<" length="<<r.length;
	  return out;
    };

  friend std::istream& operator>>(std::istream& in, RangePair& r)
    {
	  in>>r.first>>r.second>>r.score>>r.e_value>>r.identity>>r.length;
	  return in;
    };

  friend bool operator<(const RangePair& r1,const RangePair& r2)
    {
      if(r1.first<r2.first) return true;
      else
	if(r1.first==r2.first && r1.second<r2.second) return true;
	else return false;
    };
  friend bool operator>(const RangePair& r1,const RangePair& r2)
    {
      if(r1.first>r2.first) return true;
      else
	if(r1.first==r2.first && r1.second>r2.second) return true;
	else return false;
    };

  void invertQuerySubject( void )
  {
  	RangeAlign tmp = getRangeQ();
  	setRangeQ( getRangeS() );
  	setRangeS( tmp );
  	if( ! first.isPlusStrand() )
  	{
  		first.reverse();
  		second.reverse();
  	}
  }

  void reComputeSubjectCoords(RangePair& new_r, unsigned qs, unsigned qe);

  inline double computeSizeRatioOnQuery(unsigned queryStartInit, unsigned queryEndInit)
  {
    unsigned size_init = std::abs((long int)(queryEndInit-queryStartInit)) + 1;
    double ratio_self=double(first.getLength())/size_init;
    return ratio_self;
  }

  void orientSubjects();
};

#endif
