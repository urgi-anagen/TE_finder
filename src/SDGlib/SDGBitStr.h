/**
 * Class <code>SDGBitStr</code>
 *
 * The class SDGBitStr derive from SDGData.
 * It allows to encode and decode interger value on a dynamically allocated 
 * bit string.
 *
 * @author  Hadi Quesneville
 * @version 0.1, 23/06/99
 * @since   SDG1.0
 *
 */

/*
 * Eric Coissac : 29/07/99
 *    modification de la structure pour permettre la gestion par
 *    une SDGReference
 */

#ifndef SDGBITSTR_H
#define SDGBITSTR_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <climits>
#include <math.h>

#include <SDGReference.h>
#include <SDGError.h>

#ifndef LONG_BIT
#define LONG_BIT 32
#endif

class SDGBitStr;


//****************************************************************
class __SDGBitStr
{
  friend void recombine(SDGBitStr&, SDGBitStr&,
			SDGBitStr&, SDGBitStr&, int pos);
  
  friend std::ostream& operator<<(std::ostream &out,const SDGBitStr& bs);
  friend std::istream& operator>>(std::istream &fin,SDGBitStr& bs);

  typedef char BYTE;
  static const int NB_BIT=8*sizeof(BYTE);
  
  int nb_bytes;
  BYTE* byte;
  
 public:
  
  __SDGBitStr(void){ byte=NULL;init();};
  __SDGBitStr(const __SDGBitStr&);
  virtual ~__SDGBitStr() { if(byte) free(byte);};

  virtual void* clone() const
    { 
      return (void*)new __SDGBitStr(*this);
    };
  
  /**
   * Initialize a bit string with a given number of byte 
   *  @param nbbytes the number of allocated bytes
   */

  void init(int nbbytes=1);

  /**
   * Resize a bit string with a given number of byte 
   *  @param nbbytes the number of allocated bytes
   */

  void resize(int nbbytes);
  __SDGBitStr& operator=(const __SDGBitStr&); 
  
  /**
   * Decode an integer value
   *  @param pos the starting position
   *  @param length the number of bit to read
   */

  unsigned long get(int pos, int length=1) const;

  /**
   * Encode an integer value
   *  @param pos the starting position
   *  @param length the number of bit to write
   *  @param val the integer value to encode
   */

  void set(int pos, int length=1, unsigned long val=1);
  

};


class SDGBitStr : public SDGReference<__SDGBitStr>
{

  friend void recombine(SDGBitStr&, SDGBitStr&,
			SDGBitStr&, SDGBitStr&, int pos);
  
  friend std::ostream& operator<<(std::ostream &out,const SDGBitStr& bs);
  friend std::istream& operator>>(std::istream &fin,SDGBitStr& bs);

 public:

  SDGBitStr(void):SDGReference<__SDGBitStr>(){};
  SDGBitStr(const __SDGBitStr& ori):SDGReference<__SDGBitStr>(ori){};

  SDGBitStr(const SDGBitStr& ori):SDGReference<__SDGBitStr>(ori) {};
  SDGBitStr(const SDGBitStr* ori):SDGReference<__SDGBitStr>(*ori) {};

  /**
   * Initialize a bit string with a given number of byte 
   *  @param nbbytes the number of allocated bytes
   */

  void init(int nbbytes=1)
    {
      getMutablePointer()->init(nbbytes);
    };

  /**
   * Resize a bit string with a given number of byte 
   *  @param nbbytes the number of allocated bytes
   */

  void resize(int nbbytes)
    {
      getMutablePointer()->resize(nbbytes);
    };

  
  /**
   * Decode an integer value
   *  @param pos the starting position
   *  @param length the number of bit to read
   */

  unsigned long get(int pos, int length=1) const
    {
      return getPointer()->get(pos,length);
    };

  /**
   * Encode an integer value
   *  @param pos the starting position
   *  @param length the number of bit to write
   *  @param val the integer value to encode
   */

  void set(int pos, int length=1, unsigned long val=1)
    {
      getMutablePointer()->set(pos,length,val);
    };

   /****************************************
   *
   * Classes d'exceptions
   *
   ****************************************/
  
    class Exception : public SDGException 
      { 
	public :
	Exception(const void *o = NULL, const char *m= NULL, int e=0) :
	  SDGException(o,m,e) {};
      };

    class MemoryException : public Exception 
      { 
	public :
	MemoryException(const void *o = NULL, const char *m= NULL, int e=0) :
	  Exception(o,m,e) {};
      };

    class EmptyException : public Exception 
      { 
	public :
	EmptyException(const void *o = NULL, const char *m= NULL, int e=0) :
	  Exception(o,m,e) {};
      };

    class OutOfRangeException : public Exception 
      { 
	public :
	OutOfRangeException(const void *o = NULL, const char *m= NULL, int e=0) :
	  Exception(o,m,e) {};
      };
};



//==============================================================
inline __SDGBitStr::__SDGBitStr(const __SDGBitStr& bitstr)
{ 
        nb_bytes=bitstr.nb_bytes;
        if( !(byte=(char *)malloc(nb_bytes*sizeof(BYTE))) )
	    throw SDGBitStr::MemoryException(this,
		       "__SDGBitStr(const __SDGBitStr&): Malloc error !!!") ;
        memmove(byte,bitstr.byte,nb_bytes*sizeof(BYTE));
}
//--------------------------------------------------------------
inline void __SDGBitStr::init(int nbbytes)
{ 
  if(byte) free(byte);
  nb_bytes=nbbytes;
  if(nb_bytes==0) byte=NULL;
  else
    if( !(byte=(char *)malloc(nb_bytes*sizeof(BYTE))) ) 
      throw SDGBitStr::MemoryException(this,"init(): Malloc error !!!") ;
  for(int i=0;i<nb_bytes;i++) byte[i]=0;
}
//--------------------------------------------------------------
inline void __SDGBitStr::resize(int nbbytes)
{ 
  if( !(byte=(char *)realloc(byte,nbbytes*sizeof(BYTE))) ) 
    throw SDGBitStr::MemoryException(this,"resize(): realloc error !!!") ;
  for(int i=nb_bytes;i<nbbytes;i++) byte[i]=0;
  nb_bytes=nbbytes;
  if(nb_bytes==0) byte=NULL;
}
//--------------------------------------------------------------
inline __SDGBitStr& __SDGBitStr::operator=(const __SDGBitStr& bitstr)
{
  if(byte) free(byte);
  nb_bytes=bitstr.nb_bytes;
  if(nb_bytes==0) byte=NULL;
  else
    if( !(byte=(char *)malloc(nb_bytes*sizeof(BYTE))) ) 
  throw SDGBitStr::MemoryException(this,"operator=(const __SDGBitStr&): Malloc error !!!") ;

  memmove(byte,bitstr.byte,nb_bytes*sizeof(BYTE));
  return *this;
}
#endif
