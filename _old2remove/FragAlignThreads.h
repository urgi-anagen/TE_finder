/**
 *
 * FragAlignThreads.h
 *
 */

#ifndef FRAGALIGNTHREADS_H
#define FRAGALIGNTHREADS_H

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>
#include <limits>
#include <pthread.h>

#include "ThreadPool.h"
#include "RangePair.h"
#include "RangePairSet.h"

class FragAlignThreads
{
  struct bound
    {
      unsigned long end_point;
      unsigned long y_high;
      unsigned long y_low;
      bool left;
      std::list<RangePair>::iterator iter_range_pair;
      unsigned rect;
      bound():end_point(0),y_high(0),y_low(0),left(true),rect(0){};

      struct Greater 
      {
	Greater(void){};
	int operator () (const bound& a, const bound& b) const
	{
	  if(b.end_point<a.end_point) return true;
	  if(b.end_point==a.end_point && b.y_high<a.y_high) return true;
	  if(b.end_point==a.end_point && b.y_high==a.y_high 
	     && b.y_low<a.y_low) return true;
	  return false;
	};
      };

      struct Less 
      {
	Less(void){};
	int operator () (const bound& a, const bound& b) const
	{
	  if(b.end_point>a.end_point) return true;
	  if(b.end_point==a.end_point && b.y_low>a.y_low) return true;
	  if(b.end_point==a.end_point && b.y_low==a.y_low 
	     && b.y_high>a.y_high) return true;
	  return false;
	};
      };

      static const Greater greater;
      static const Less less;

  };

  struct triple
    {
      unsigned long lower;
      long long path_val;
      std::list<RangePair>::iterator iter_range_pair;
      unsigned rect;
      triple():lower(0),path_val(0),rect(0){};

      struct Less 
      {
	Less(void){};
	int operator () (const triple& t1, const triple& t2) const
	{return t1.lower<t2.lower;};
      };
      struct LessI 
      {
	LessI(void){};
	int operator () (const triple& t, const long unsigned int& i) const
	{return t.lower<i;};
      };

      struct Greater 
      {
	Greater(void){};
	int operator () (const triple& t1, const triple& t2) const
	{return t1.lower>t2.lower;};
      };
      struct Equal
      {
	Equal(void){};
	int operator () (const triple& t1, const triple& t2) const
	{return t1.lower==t2.lower;};
      };

      static const Greater greater;
      static const Less less;
      static const LessI lessI;
      static const Equal equal;

  };



  double gapo_pen;
  double gape_pen;
  double mism_pen;
  unsigned over;
  unsigned connect_dist_limit;

  std::list<bound> I; // list of horizontal rectangle coordinates (left and right corners)
  std::vector<long long> V; // vector of the best chain scores that finish in a given rectangle  curr_path.clear();

  //long long align(std::list<bound> *pI, std::vector<long long> *pV, std::list<std::list<RangePair>::iterator> *pPath, unsigned nb_frag );
  void align(double go,double ge, double m, unsigned c, std::list<bound>& i,
		  std::vector<long long>& v, unsigned nf,
		  std::list<std::list<RangePair>::iterator> &path_result);
 public:

  FragAlignThreads(double mism,double gapo_p,double gape_p, unsigned o,
	    unsigned l=20000)
    {
      gapo_pen=gapo_p;
      gape_pen=gape_p;
      mism_pen=mism;
      connect_dist_limit=l;
      over=o;
    };

 	 ~FragAlignThreads(void){}

  void alignDirectDirect(std::list<RangePair>& l, std::list<std::list<RangePair>::iterator>& path_result);
  void alignDirectCompl(std::list<RangePair>& l, std::list<std::list<RangePair>::iterator>& path_result);
  void alignComplDirect(std::list<RangePair>& l, std::list<std::list<RangePair>::iterator>& path_result);
  void alignComplCompl(std::list<RangePair>& l, std::list<std::list<RangePair>::iterator>& path_result);
  void join(std::list<RangePair>& l,std::list<RangePairSet>& finished_list);
};
#endif





