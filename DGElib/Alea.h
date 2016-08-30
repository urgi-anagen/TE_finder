/*
 * Copyright (c) Laboratoire de Dynamique du Genome et Evolution - Institut 
 * Jacques Monod, 1997. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by Hadi Quesneville at Laboratoire de Dynamique 
 * du Genome et Evolution, Institut Jacques Monod, 2 place Jussieu, 
 * 75251 Paris Cedex 05, France.
 * 
 * Contact: Hadi Quesneville
 * Laboratoire de Dynamique du Genome et Evolution,
 * Institut Jacques Monod,
 * 2, place Jussieu, 75251 Paris Cedex 05, France
 * e-mail: hq@ccr.jussieu.fr
 *
 *
 * Laboratoire de Dynamique du Genome et Evolution disclaims all warranties
 * with regard to this software.
 */
//----------------------------------------------------------------------------
// RANDOM.H - contains random number generator and related utilities.
//----------------------------------------------------------------------------
#ifndef random_h
#define random_h 1

#include <limits.h>
//#include <ACG.h>
//#include <Random.h>
#include <vector>

#include <SDGError.h>
#include <SDGString.h>
#include "CountStat.h"

class Alea
{
  //  enum {mACG,mSMG} methode;
  
  // Substractive Methode Generator
  double oldrand[55];  /* Array of 55 random numbers */
  int jrand;           /* current random number */
  double rndx2;        /* used with random normal deviate */
  int rndcalcflag;     /* used with random normal deviate */
  unsigned long count;
  unsigned seed;
  
  void warmupRandom(double random_seed);
  double randomNormalDeviate(void);
  void initRandomNormalDeviate(void);
  
  //  Additive Congruential Generator     
  // ACG randGenerator;

public:
  //    Alea(void): randGenerator(){methode=mACG;init();};
  //	Alea(unsigned seed):randGenerator(){init(seed);};
        Alea(void) {init();};
	Alea(unsigned seed){init(seed);};
        void init(void);
        void initParam(SDGString);
	void showParam(std::ostream& out);
        void initByTime(void);
	void init(unsigned seed);

	unsigned long getCount(){return count;};

        int flip(double prob) /* Flip a biased coin - true if heads */
	  {
	    return (randomPerc() <= prob);
	  };

	unsigned wheel(const std::vector<double>& proba)
	  {
	    unsigned i;
	    double perc=randomPerc(),cumProba=0; 
	    for(i=0;i<proba.size();i++)
	      {
		cumProba+=proba[i];
		 if(perc<cumProba)
		  return i;
	      }
	    if(perc==cumProba && perc==1.0) return i-1;
	    return i;
	  };

        double noise(double mu ,double sigma)
	/* normal noise with specified mean & std dev: mu & sigma */
	{
	  return((randomNormalDeviate()*sigma) + mu);
	};

	int randomBinomial(double p,int n)
	  {
	    int s = 0;
	    for (int i = 0; i < n; i++) {
	      if (randomPerc() < p) {
		s++;
	      }
	    }
	    return(s);
	  };

	int randomPoisson(double m)
	  {
	    double bound = exp(-1.0 * m);
	    int count = 0;

	    for (double product = 1.0;
		 product >= bound;
		 product *= randomPerc() ) {
	      count++;
	    }
	    return(count - 1);
	  };

	int randomGeometric(double p)
	  {
	    int count;
	    double q=1-p;
	    for (count = 1; randomPerc() < q; count++);
	    return count;
	  };

	long rnd(long arg)
	  {
	    return( rndLoHi(0,arg-1) );
	  };

	long rndLoHi(long low, long high)
	  /* Pick a random integer between low and high */
	  {
	    long i;

	    if(low >= high)
	      i = low;
	    else
	      {
		i = (long)((randomPerc() * (high - low + 1)) + low);
		if(i > high) i = high;
	      }
	    return(i);
	  };

	double rndRealLoHi(double lo ,double hi)
	  /* real random number between specified limits */
	  {
	    return((randomPerc() * (hi - lo)) + lo);
	  };

	double randomPerc(void)
	  {
	    count++;
/* 	    if(methode==mACG) */
/* 	      { */
/* 		  return randGenerator.asDouble(); */
/* 	      } */
	    // methode SMG
	    // Fetch a single random number between 0.0 and 1.0 
	    // Subtractive Method See Knuth, D. (1969), v. 2 for details 
	    if(count==ULONG_MAX) 
	      std::cerr<<"number of random numbers drawn exceed long size!\n";
	    jrand++;
	    if(jrand >= 55)
	      {
		jrand = 1;
		advanceRandom();
	      }
	    return(oldrand[jrand]);
	  };

	unsigned getSeed(){return seed;};
        void advanceRandom(void);
};

extern Alea alea;

class RandomIndex : public std::vector<unsigned>
{

public:

  RandomIndex(){};
  RandomIndex(std::vector<unsigned> vec)
    {
      for(std::vector<unsigned>::iterator i=vec.begin();i!=vec.end();i++)
	add(*i);
    }

  void add(unsigned i) { push_back(i);};
  unsigned draw() 
    {      
      if(empty()) 
	throw SDGException(this,"RandomIndex::draw(): empty!!");
      int i=alea.rnd(size());
      unsigned idx=operator[](i);
      operator[](i)=back();
      pop_back();
      return idx;
    };
  unsigned get() 
    {      
      if(empty()) 
	throw SDGException(this,"RandomIndex::get(): empty!!");
      int i=alea.rnd(size());
      unsigned idx=operator[](i);
      return idx;
    };
  void del(unsigned i)
    {
      std::vector<unsigned>::iterator idx;
      for(idx=begin();idx!=end() && (*idx)!=i;idx++);
      if(idx!=end())
	{
	  (*idx)=back();
	  pop_back();
	}
    }
};
#endif
