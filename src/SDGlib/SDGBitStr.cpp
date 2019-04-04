/*
 *
 * Definition d'une classe SDGBitStr
 *
 *
 */

#include "SDGBitStr.h"
//***************************************************************
std::istream& operator>>(std::istream &in, SDGBitStr& bs)
{
  char c=0;
  int i=0;

  while ((!(in.eof())) && (c != '\n'))
	{
	  in.get(c);	 
	  bs.set(i++,0,atoi(&c));
	}
  return in;
}
//-----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream &out, const SDGBitStr& bs)
{
  for(int i=0;i<(__SDGBitStr::NB_BIT * bs.getPointer()->nb_bytes);i++)
    out<<bs.get(i);
  return out;
}
//---------------------------------------------------------------
unsigned long __SDGBitStr::get(int i, int l) const
{
  if(!byte) throw SDGBitStr::EmptyException(this,"get():  empty object !!!") ;
  if(--l<0) throw SDGBitStr::EmptyException(this,"get():  length <= 0 !!!") ;
  

  // bytes number
  int iisin = i/NB_BIT;
  if(iisin > nb_bytes)
    throw SDGBitStr::OutOfRangeException(this,"get(): start out of range !!");
  
  int jisin = (i+l)/NB_BIT;
  if(jisin > nb_bytes)
    throw SDGBitStr::OutOfRangeException(this,"get(): end out of range !!");
  
  unsigned long val=0;

  // Relative bit positions in the bytes
  int i1 = i - (iisin*NB_BIT);
  int j1 = (i+l) - (jisin*NB_BIT);
  
  if(iisin == jisin) //same byte?
    {
      BYTE mask = 1;
      mask<<=(j1-i1+1);
      mask--;
      mask<<=i1;
      val=(byte[iisin]&mask)>>i1;
    }
  else //different byte
    {
      BYTE mask_end = 1;
      mask_end<<=(j1+1);
      mask_end--;
      //   BYTE val_end = byte[jisin]&mask_end;
      BYTE mask_start = 1;
      mask_start<<=(NB_BIT-i1);
      mask_start--;
      mask_start<<=i1;
      BYTE val_start=(byte[iisin]&mask_start)>>i1;
      val=val_start;
      int bit_pos=NB_BIT-i1;
      for(int w=iisin+1;w<=jisin;w++)
	{
	  val|=byte[w]<<bit_pos;
	  bit_pos+=NB_BIT;
	}
    }
  return val;
}
//--------------------------------------------------------------
void __SDGBitStr::set(int i, int l, unsigned long val)
{
  if(--l<0) throw SDGBitStr::EmptyException(this,"set():  length <= 0 !!!") ;

  // value control
   if( l>= LONG_BIT )
    throw SDGBitStr::OutOfRangeException(this,
       	      "set(): length out of the range of an unsigned long int!!");

   if(((unsigned long)(pow(2,l+1)-1)&val)<val)
    throw SDGBitStr::OutOfRangeException(this,"set(): value out of range !!");

  // bytes number
  int iisin = i/NB_BIT;
  if(iisin>=nb_bytes) resize(iisin+1);
  int jisin = (i+l)/NB_BIT;
  if(jisin>=nb_bytes) resize(jisin+1);
  
  // Relative bit positions in the bytes
  int i1 = i - (iisin*NB_BIT);
  int j1 = (i+l) - (jisin*NB_BIT);
  
  // same byte?
  if(iisin == jisin)
    {
      BYTE mask = 1;
      mask<<=(j1-i1+1);
      mask--;
      BYTE mskval=val&mask;

      mask<<=i1;
      byte[iisin]=byte[iisin]&~mask;
      
      mskval<<=i1;
      byte[iisin]=byte[iisin]|mskval;
    }
  else
    {
      // start byte
      BYTE mask = 1;
      mask<<=(NB_BIT-i1);
      mask--;
      BYTE mskval=val&mask;
      mask<<=i1;
      byte[iisin]=byte[iisin]&~mask;
      
      mskval<<=i1;
      byte[iisin]=byte[iisin]|mskval;

      val>>=NB_BIT-i1;
      for(int w=iisin+1;w<jisin;w++)
	{
	  byte[w]=0;
	  byte[w]=byte[w]|(val&255);
	  val>>=NB_BIT;
	}

      // end byte
      mask = 1;
      mask<<=(j1+1);
      mask--;
      byte[jisin]=byte[jisin]&~mask;
      
      byte[jisin]=byte[jisin] | (val&mask);
      
    }
}
//---------------------------------------------------------------
void recombine(SDGBitStr& parent1, SDGBitStr& parent2, SDGBitStr& child1, SDGBitStr& child2, int jcross)
{
    if (!parent1.getPointer()->byte)
        throw SDGBitStr::EmptyException(&parent1, "recombine():  empty parent1 object !!!");
    if (!parent2.getPointer()->byte)
        throw SDGBitStr::EmptyException(&parent2, "recombine():  empty parent2 object !!!");

    if (parent1.getPointer()->nb_bytes < parent2.getPointer()->nb_bytes)
        parent1.resize(parent2.getPointer()->nb_bytes);

    if (parent1.getPointer()->nb_bytes > parent2.getPointer()->nb_bytes)
        parent2.resize(parent1.getPointer()->nb_bytes);

    int nbbytes = parent1.getPointer()->nb_bytes;

    int k = (jcross / __SDGBitStr::NB_BIT);
    if (k > nbbytes)
        throw SDGBitStr::OutOfRangeException(&parent1, "recombine(): recombine point out of range !!");


    memmove(child1.getPointer()->byte,
            parent1.getPointer()->byte, k * sizeof(__SDGBitStr::BYTE));

    memmove(child2.getPointer()->byte,
            parent2.getPointer()->byte, k * sizeof(__SDGBitStr::BYTE));

    if (jcross - (k * __SDGBitStr::NB_BIT) == 0) k--;
    else {
        __SDGBitStr::BYTE mask = 1;
        mask <<= (jcross - (k * __SDGBitStr::NB_BIT));
        mask--;
        child1.getPointer()->byte[k] =
                (parent1.getPointer()->byte[k] & mask) |
                (parent2.getPointer()->byte[k] & (~mask));

        child2.getPointer()->byte[k] =
                (parent1.getPointer()->byte[k] & (~mask)) |
                (parent2.getPointer()->byte[k] & mask);
    }

    memmove(&child1.getPointer()->byte[k + 1],
            &parent2.getPointer()->byte[k + 1],
            (nbbytes - k - 1) * sizeof(__SDGBitStr::BYTE));

    memmove(&child2.getPointer()->byte[k + 1],
            &parent1.getPointer()->byte[k + 1],
            (nbbytes - k - 1) * sizeof(__SDGBitStr::BYTE));
}








