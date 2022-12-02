/**
 *
 * FragJoin.h
 *
 */

#ifndef FRAGJOIN_H
#define FRAGJOIN_H

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

class FragJoin
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

  void align(unsigned nb_frag, std::list<bound>& I, std::vector<long long>& V, std::list<std::list<RangePair>::iterator>& join_path);

 public:

    FragJoin(double mism, double gapo_p, double gape_p, unsigned o,
              unsigned l = 20000) {
        gapo_pen = gapo_p;
        gape_pen = gape_p;
        mism_pen = mism;
        connect_dist_limit = l;
        over = o;
    };


    void splitFromStrand(std::list<RangePair> &list_in,
                         std::list<RangePair> &list_out_dd,
                         std::list<RangePair> &list_out_id,
                         std::list<RangePair> &list_out_di,
                         std::list<RangePair> &list_out_ii);

    void align_dd(std::list<RangePair> &l_in, std::list<RangePairSet> &l_out);

    void align_dc(std::list<RangePair> &l_in, std::list<RangePairSet> &l_out);

    void align_cd(std::list<RangePair> &l_in, std::list<RangePairSet> &l_out);

    void align_cc(std::list<RangePair> &l_in, std::list<RangePairSet> &l_out);

    void align_all(std::list<RangePair> &l, std::list<RangePairSet> &l_out);

    void join_path2rp_list(std::list<std::list<RangePair>::iterator>& join_path, std::list<RangePair> &l_in, std::list<RangePairSet> &l_out);
};
#endif





