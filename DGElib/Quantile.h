#ifndef QUANTILE_H
#define QUANTILE_H 1

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "CountStat.h"

class Quantile : public CountStat
{
  std::vector<double> values_list;
  
  void order(void){ sort(values_list.begin(),values_list.end());};

 public:
  Quantile(void):CountStat(){};

  void init() {CountStat::init();values_list.clear();};
  void reset() {init();};

  void operator+(double val){add(val);};

  void add(double d)
    {
      CountStat::add(d);
      values_list.push_back(d);
    };
  
  double quantile(double q)
    {
      order();
      int pos=(int)floor(q*size());
      return values_list[pos];
    };

  int locate(double q)
    {
      order();
      int count=0;
      std::vector<double>::iterator i;
      for(i=values_list.begin();i!=values_list.end() && *i<q;i++)
	count++;
      if(i==values_list.end()) count=size()-1;
      return count;
    };
};
#endif
