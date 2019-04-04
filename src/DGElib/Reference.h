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
// Reference.h 
//----------------------------------------------------------------------------

#ifndef REFERENCE_H
#define REFERENCE_H

#include <stddef.h>
#include <iostream>
#include <SDGError.h>

template <class T>
class Reference{

  struct Ref
  {
    T* ptr;
    unsigned count;
    
    Ref(T* p):ptr(p),count(1){};
  };

 protected:
 
  Ref* ref;

 public:


  Reference()
    {
      ref=NULL;
    };

  Reference(T* p)
    {
      ref=new Ref(p);
    };

  Reference(Ref* r)
    {
      if(r)
	r->count++;
      ref=r;
    };
  
  Reference(const Reference<T>& robj)
    {
      if(robj.ref)
	robj.ref->count++;
      ref=robj.ref;
    };

  virtual ~Reference()
    {
      clear();
    };

  virtual void clear()
    {
      if(ref && --(ref->count)==0)
	{
	  if(ref->ptr) delete ref->ptr;
	  delete ref;
	}
      ref=NULL;
    }

  virtual Reference<T>  operator=(T* p)
    {      
      clear();
      ref=new Ref(p);
      return *this;
    };

  virtual  Reference<T>  operator=(const Reference<T>& robj)
    {
      if(robj.ref)
	robj.ref->count++;
      clear();
      ref=robj.ref;
      return *this;
    };

  bool isAllocated()
    {
      if(!ref) return false;
      if(ref->ptr) return true;
      else return false;
    };

  bool isAllocated() const
    {
      if(!ref) return false;
      if(ref->ptr) return true;
      else return false;
    };

    T *getPtr() {
        if (!ref)
            throw SDGException(this,
                               "Reference<T>::getPtr():  empty reference !!!");
        return ref->ptr;
    };

    const T *getPtr() const {
        if (!ref)
            throw SDGException(this,
                               "Reference<T>::getPtr():  empty reference !!!");
        return ref->ptr;
    };

    T &getObj() {
        if (!ref)
            throw SDGException(this,
                               "Reference<T>::getObj():  empty reference !!!");
        return *(ref->ptr);
    };

    const T &getObj() const {
        if (!ref)
            throw SDGException(this,
                               "const Reference<T>::getObj():  empty reference !!!");
        return *(ref->ptr);
    };

    T *operator->() {
        if (!ref)
            throw SDGException(this,
                               "Reference<T>::operator->():  empty reference !!!");
        return ref->ptr;
    };

    const T *operator->() const {
        if (!ref)
            throw SDGException(this,
                               "Reference<T>::operator->():  empty reference !!!");
        return ref->ptr;
    };

    Reference<T> cloneRef() {
        if (ref) {
            return Reference((T *) ref->ptr->clone());
        } else
            return Reference();
    };

    void makeExclusif() {
        Reference<T> old = *this;
        clear();
        *this = old.cloneRef();
    };

 friend int operator==(const Reference<T>& r1, const Reference<T>& r2)
    {
      if(*(r1.ref->ptr)==*(r2.ref->ptr)) return 1;
      return 0;	 
    };

 friend int operator!=(const Reference<T>& r1, const Reference<T>& r2)
    {
      if(*(r1.ref->ptr)!=*(r2.ref->ptr)) return 1;
      return 0;	 
    };

 friend int operator>(const Reference<T>& r1, const Reference<T>& r2)
    {
      if(*(r1.ref->ptr)>*(r2.ref->ptr)) return 1;
      return 0;	 
    };

 friend int operator>=(const Reference<T>& r1, const Reference<T>& r2)
    {
      if(*(r1.ref->ptr)>=*(r2.ref->ptr)) return 1;
      return 0;	 
    };

 friend int operator<(const Reference<T>& r1, const Reference<T>& r2)
    {
      if(*(r1.ref->ptr)<*(r2.ref->ptr)) return 1;
      return 0;	 
    };

 friend int operator<=(const Reference<T>& r1, const Reference<T>& r2)
    {
      if(*(r1.ref->ptr)<=*(r2.ref->ptr)) return 1;
      return 0;	 
    };


 friend std::ostream& operator<<(std::ostream& os,const Reference<T> &robj)
   {
     if(robj.ref)
       {
	 os<<"[count="<<robj.ref->count<<"]";
	 //	 return operator<<(os,*(robj.ref->ptr));
	 return os;
       }
     else
       os<<"(ref NULL!!)";
   };
 
 friend std::istream& operator>>(std::istream& is,const Reference<T> &robj)
   {
      if(!robj.ref)  
	throw SDGException(&robj,
	      "Reference<T>::operator>>():  empty reference !!!") ;

	 return operator>>(is,*(robj.ref->ptr));
   };
};
#endif







