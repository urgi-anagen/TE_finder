#ifndef COUNTSTAT_H
#define COUNTSTAT_H 1

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

//-----------------------------------------------------------------------------
class CountStat
{
      double som,som2;
      double _min,_max;
      int n;

public:
      CountStat() {init();};
      void init() {_min=_max=som=som2=0.0; n=0;};
      void reset() {init();};

      friend int operator==(const CountStat& s1, const CountStat& s2)
      {
         if(s1.som!=s2.som) return 0;
         if(s1.som2!=s2.som2) return 0;
         if(s1._min!=s2._min) return 0;
         if(s1._max!=s2._max) return 0;
         if(s1.n!=s2.n) return 0;
	 return 1;
      }

      friend int operator!=(const CountStat& s1, const CountStat& s2)
      {
         if(s1==s2) return 0;
      	 return 1;
      }
      void operator+=(double val){add(val);};
      void add(double val)
	{
	  if(n==0) _min=_max=val;
	  else
	    {
	      if(_min>val) _min=val;
	      if(_max<val) _max=val;
	    }
	  som+=val;
	  som2+=(val*val);
	  n++;
	};

      double sum() const {return som;}

      double mean() const 
      {if(n<1) return -1; else return som/n;};

      double var() const 
      {if(n<2) return 0; else return (som2-((som*som)/n))/(n-1);};

      double sd() const {return sqrt(var());};
      double cv() const {if(som==0) return 0; else return sd()/mean();};
      int size() const {return n;};
      double getMin(void) const {return _min;};
      double getMax(void) const {return _max;};

      void save(std::ofstream& file)
	{
	  file<<_min<<"\n";
	  file<<_max<<"\n";
	  file<<som<<"\n";
	  file<<som2<<"\n";
	  file<<n<<"\n";
	};

      void load(std::ifstream& file)
	{
	  char buff[81];
	  
	  file.getline(buff,80);
	  _min=atof(buff);
	  
	  file.getline(buff,80);
	  _max=atof(buff);
	  
	  file.getline(buff,80);
	  som=atof(buff);
	  
	  file.getline(buff,80);
	  som2=atof(buff);
	  
	  file.getline(buff,80);
	  n=atoi(buff);
	};

      friend std::ostream& operator<<(std::ostream& out, const CountStat& s)
      {
	out<<"mean="<<s.mean();
	out<<" var="<<s.var();
	out<<" sd="<<s.sd();
	out<<" cv="<<s.cv();
	out<<" min="<<s.getMin();
	out<<" max="<<s.getMax();
	out<<" n="<<s.size();
	return out;
      };

};
#endif

