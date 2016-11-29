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
// RANDOM.C - contains random number generator and related utilities.
//----------------------------------------------------------------------------
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <new>

#include "Alea.h"

Alea alea;
//*****************************************************************************
void Alea::advanceRandom(void)
/* Create next batch of 55 random numbers */
{
    int j1;
    double new_random;

    for(j1 = 0; j1 < 24; j1++)
    {
        new_random = oldrand[j1] - oldrand[j1+31];
        if(new_random < 0.0) new_random = new_random + 1.0;
        oldrand[j1] = new_random;
    }
    for(j1 = 24; j1 < 55; j1++)
    {
        new_random = oldrand [j1] - oldrand [j1-24];
        if(new_random < 0.0) new_random = new_random + 1.0;
        oldrand[j1] = new_random;
    }
}
//----------------------------------------------------------------------------
void Alea::initRandomNormalDeviate(void)
/* initialization routine for randomnormaldeviate */
{
    rndcalcflag = 1;
}
//---------------------------------------------------------------------------
void Alea::init(void)
{
    for(int j1=0; j1<=54; j1++) oldrand[j1] = 0.0;
    jrand=0;
    count=0;
}
//---------------------------------------------------------------------------
void Alea::initParam(SDGString param )
{
  if(param.beforematch("=")=="SEED")
    {
      SDGString sseed;
      SDGString val=param.aftermatch("=");
      long int lenmatch;
      if(val.match("/",lenmatch))
	{
	  sseed=val.beforematch("/");
// 	  if(val=="ACG") 
// 	    {
// 	      methode=mACG;
// 	    }
// 	  if(val=="SMG") 
// 	    {
// 	      methode=mSMG;
// 	    }
	}
      else sseed=val;
      if(sseed == "TIME") initByTime();
      else init(atoi(sseed.start()));
    }
}
//---------------------------------------------------------------------------
void Alea::showParam(std::ostream& out)
{
  out<<"\nRandom generator parameters:"<<std::endl;
  out<<"----------------------------"<<std::endl;
  out<<"  SEED="<<seed<<std::endl;
//   out<<"  RNG=";
//   switch(methode)
//     {
//         case mACG : out<<"ACG"<<std::endl; break;
//         case mSMG : out<<"SMG"<<std::endl; break;
//     }
}
//---------------------------------------------------------------------------
void Alea::initByTime(void)
{
  time_t  now;
  time(&now);
  init(now);
}
//----------------------------------------------------------------------------
void Alea::init(unsigned s)
{
// randGenerator.setSeed(s);
//   if(methode==mSMG)
//     {
//       for(int j1=0; j1<=54; j1++) oldrand[j1] = 0.0;
//       jrand=0;
//       warmupRandom(randGenerator.asDouble());
//     }
  count=0;
  seed=s;
  for(int j1=0; j1<=54; j1++) oldrand[j1] = 0.0;
  jrand=0;
  warmupRandom((double)seed);
}
//----------------------------------------------------------------------------
double Alea::randomNormalDeviate(void)
/* random normal deviate after ACM algorithm 267 / Box-Muller Method */
{
    double t, rndx1;

    if(rndcalcflag)
    {
        rndx1 = sqrt(- 2.0*log(randomPerc()));
        t = 6.2831853072 * randomPerc();
        rndx2 = rndx1 * sin(t);
        rndcalcflag = 0;
        return(rndx1 * cos(t));
    }
    else
    {
        rndcalcflag = 1;
        return(rndx2);
    }
}
//----------------------------------------------------------------------------
void Alea::warmupRandom(double random_seed)
/* Get random off and running */
{
    int j1, ii;
    double new_random, prev_random;

    oldrand[54] = random_seed;
    new_random = 0.000000001;
    prev_random = random_seed;
    for(j1 = 1 ; j1 <= 54; j1++)
    {
        ii = (21*j1)%54;
        oldrand[ii] = new_random;
        new_random = prev_random-new_random;
        if(new_random<0.0) new_random = new_random + 1.0;
        prev_random = oldrand[ii];
    }

    advanceRandom();
    advanceRandom();
    advanceRandom();

    jrand = 0;
}
//*****************************************************************************

