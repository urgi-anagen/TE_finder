/**
 *
 * FragAlign.h
 *
 */

#ifndef FRAGALIGN_H
#define FRAGALIGN_H

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>
#include <limits>

#include "RangePair.h"
#include "RangePairSet.h"

class FragAlign
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


  std::list<bound> I; // list of horizontal rectangle coordinates (left and right corners)
  unsigned nb_frag;
  std::vector<long long> V; // vector of the best chain scores that finish in a given rectangle
  std::list<std::list<RangePair>::iterator> curr_path, 
    path_dd, path_cd, path_dc, path_cc; 
  unsigned k;

  double gapo_pen;
  double gape_pen;
  double mism_pen;
  unsigned over;
  unsigned connect_dist_limit;

  void align(void);

 public:

    FragAlign(double mism, double gapo_p, double gape_p, unsigned o,
              unsigned l = 20000) {
        gapo_pen = gapo_p;
        gape_pen = gape_p;
        mism_pen = mism;
        connect_dist_limit = l;
        over = o;
    };

    void alignDirectDirect(std::list<RangePair> &l);

    void alignDirectCompl(std::list<RangePair> &l);

    void alignComplDirect(std::list<RangePair> &l);

    void alignComplCompl(std::list<RangePair> &l);

    std::list<RangePairSet> join(std::list<RangePair> &l);

};
#endif





