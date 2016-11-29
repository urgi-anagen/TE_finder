/*
 some useful hash functions
*/
#ifndef HASH_H
#define HASH_H

#include <math.h>
#include <string>
#include <functional>

using std::string;

template <class string > struct HashDNA 
: public std::unary_function< string , unsigned >
{

  int nuc_val[256];

  HashDNA(void)
  {  
    nuc_val[(unsigned)'A']=0;
    nuc_val[(unsigned)'C']=1;
    nuc_val[(unsigned)'G']=2;
    nuc_val[(unsigned)'T']=3;
  };

  unsigned operator()(const std::string& word) const 
  {
    const char* x=word.c_str();
    int i=0;
    unsigned int h=0;
    while (*x != 0) h+=(unsigned)pow((double)4,i++)*nuc_val[(unsigned)*x++];
    return h;
  };
};

inline unsigned int hashpjw(const char* x) // From Dragon book, p436
{
  unsigned int h = 0;
  unsigned int g;

  while (*x != 0)
  {
    h = (h << 4) + *x++;
    if ((g = h & 0xf0000000) != 0)
      h = (h ^ (g >> 24)) ^ g;
  }
  return h;
};

inline unsigned int multiplicativehash(int x)
{
  // uses a const close to golden ratio * pow(2,32)
  return ((unsigned)x) * 2654435767U;
};


inline unsigned int foldhash(double x)
{
  union { unsigned int i[2]; double d; } u;
  u.d = x;
  unsigned int u0 = u.i[0];
  unsigned int u1 = u.i[1]; 
  return u0 ^ u1;
};
#endif
