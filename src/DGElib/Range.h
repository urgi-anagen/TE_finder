/**
 * \file Range.h
 * \brief Header file for the class Range
 */

/*! \mainpage Documentation of the C++ code in the REPET package
 *
 * \section intro_sec Introduction
 *
 * The REPET package integrates bioinformatics programs in order to tackle biological issues at the genomic scale.
 * It is distributed under the CeCILL license  and deposited to the Agence de Protection des Programmes (APP).
 * Code is written in Python and C++, but this documentation only deals with C++ code.
 *
 * The C++ part of the package contains three main programs, BLASTER, GROUPER and MATCHER.
 * These programs take advantage of three libraries, BLRlib, DGElib and SDGlib.
 * Several tools are also available.
 *
 * \section install_sec Installation
 *
 * To install the whole REPET package, enter into the main directory and type "make".
 * Then, follow the instructions, but you may type "make compile" and "make install".
 *
 * If you are only interested in the C++ binaries, enter into the "TE_finder/" directory.
 * Then, type "make" and "make install". All the binaries will be in the "bin/" directory.
 *
 */

#ifndef RANGE_H
#define RANGE_H

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <string>
#include <ostream>


typedef unsigned long ulong;
/**
 * \class Range
 * \brief Implement a pair of coordinate
 *
 * This class has two attributes, start and end.
 */
class Range
{
  struct Less
  {
  	Less(void){};
  	int operator () (const Range& a, const Range& b) const
  	{
  		if(a<b) return true;
  		return false;
  	};
  };

  struct Greater
  {
    Greater(void){};
    int operator () (const Range& a, const Range& b) const
    {
      if(a>b) return true;
      return false;
    };
  };

  struct GreaterLength
  {
    GreaterLength(void){};
    int operator () (const Range& a, const Range& b) const
    {
      if(a.getLength()>b.getLength())
	return true;
      return false;
    };
  };

	public:
		static const Less less;
		static const Greater greater;
		static const GreaterLength greaterLength;

	protected:

		unsigned long start,end;

	public:
		Range( const Range *o );
		Range( const Range &o );
		Range( ulong s=0, ulong e=0 );

    virtual ~Range(void);
    virtual void *clone(void) const;

    void set( ulong s=0, ulong e=0 )
    {
    	start = s;
    	end = e;
    };

    ulong getStart( void ) const { return start; };

    void setStart( ulong s ){ start=s; };

    ulong getEnd( void ) const { return end; };

    void setEnd( ulong e ){ end=e; };

    int isPlusStrand( void ) const { return start<end; };

    ulong getMin( void ) const { return std::min( start, end ); };

    ulong getMax( void ) const { return std::max( start, end ); };

    unsigned long getLength( void ) const
    {
    	return getMax()-getMin()+1;
    };

    const void view( void ) const
    {
    	std::cout<<"range: start="<<getStart()<<", end="<<getEnd()<<std::endl;
    }

    void reverse( void )
    {
    	ulong tmp = start;
    	start = end;
    	end = tmp;
    };

    bool empty( void )
    {
    	return ( start==0 && end==0 );
    };

    bool empty( void ) const
    {
    	return ( start==0 && end==0 );
    };

    bool overlap( const Range& r ) const
    {
    	ulong s1 = std::min(start,end);
    	ulong e1 = std::max(start,end);
    	ulong s2 = std::min(r.start,r.end);
    	ulong e2 = std::max(r.start,r.end);

    	if( s1>=s2 && s1<=e2 ) return true;
    	if( s2>=s1 && s2<=e1 ) return true;

    	return false;
    };

    bool overlap( const Range& r )
    {
    	ulong s1 = std::min(start,end);
    	ulong e1 = std::max(start,end);
    	ulong s2 = std::min(r.start,r.end);
    	ulong e2 = std::max(r.start,r.end);

    	if( s1>=s2 && s1<=e2 ) return true;
    	if( s2>=s1 && s2<=e1 ) return true;

    	return false;
    };

    bool isContained(const Range& r) const
    // object is contained in r
    {
    	ulong s1=std::min(start,end);
    	ulong e1=std::max(start,end);
    	ulong s2=std::min(r.start,r.end);
    	ulong e2=std::max(r.start,r.end);
    	if(s2<=s1 && e2>=e1) return true;
    	return false;
    };

    bool isContained(const Range& r)
    // object is contained in r
    {
    	ulong s1=std::min(start,end);
    	ulong e1=std::max(start,end);
    	ulong s2=std::min(r.start,r.end);
    	ulong e2=std::max(r.start,r.end);
    	if(s2<=s1 && e2>=e1) return true;
    	return false;
    };

    bool isStrictlyContained(const Range& r) const
    // object is strictly contained in r
    {
    	ulong s1=std::min(start,end);
    	ulong e1=std::max(start,end);
    	ulong s2=std::min(r.start,r.end);
    	ulong e2=std::max(r.start,r.end);
    	if(s2<s1 && e2>e1) return true;
    	return false;
    };

    bool isStrictlyContained(const Range& r)
    // object is strictly contained in r
    {
    	ulong s1=std::min(start,end);
    	ulong e1=std::max(start,end);
    	ulong s2=std::min(r.start,r.end);
    	ulong e2=std::max(r.start,r.end);
    	if(s2<s1 && e2>e1) return true;
    	return false;
    };

    bool isIncluded(const Range& r) const
    {
    	// object contains r
    	ulong s1=std::min(start,end);
    	ulong e1=std::max(start,end);
    	ulong s2=std::min(r.start,r.end);
    	ulong e2=std::max(r.start,r.end);

    	if(s1<=s2 && e1>=e2) return true;

    	return false;
    };

    bool isIncluded(const Range& r)
    {
    	// object contains r
    	ulong s1=std::min(start,end);
    	ulong e1=std::max(start,end);
    	ulong s2=std::min(r.start,r.end);
    	ulong e2=std::max(r.start,r.end);

    	if(s1<=s2 && e1>=e2) return true;

    	return false;
    };

    long distance(const Range& r)
    {
    	return std::min(r.start,r.end)-std::min(start,end);
    };

    long distance(const Range& r) const
    {
    	return std::min(r.start,r.end)-std::min(start,end);
    };

    Range diff(const Range& r)
    {
    	// this is adjusted according to argument
    	Range new_range;
    	if( overlap(r) )
    	{
    		ulong is=std::min(start,end);
    		ulong ie=std::max(start,end);
    		ulong js=std::min(r.start,r.end);
    		ulong je=std::max(r.start,r.end);

    		if(is<js)
    		{
    			if(ie<je)
    			{
    				if(isPlusStrand())
    				{
    					start=is;
    					end=js-1;
    				}
    				else
    				{
    					start=js-1;
    					end=is;
    				}
    			}
    			else // r is included
    			{
    				if(isPlusStrand())
    				{
    					start=is;
    					end=js-1;
    					new_range.set(je+1,ie);
    				}
    				else
    				{
    					start=js-1;
    					end=is;
    					new_range.set(ie,je+1);
    				}
    			}
    		}
    		else //is>=js
    		{
    			if(ie<=je)
    			{
    				set();
    			}
    			else
    			{
    				if(isPlusStrand())
    				{
    					start=je+1;
    					end=ie;
    				}
    				else
    				{
    					start=ie;
    					end=je+1;
    				}
    			}
    		}
    	} // if
    	return new_range;
    };

    /*
     * \fn void merge( const Range& r, bool keep_strand=true )
     * \brief merge two Range instances
     */
    void merge( const Range& r, bool keep_strand=true )
    {
    	ulong s1 = std::min(start,end);
    	ulong e1 = std::max(start,end);
    	ulong s2 = std::min(r.start,r.end);
    	ulong e2 = std::max(r.start,r.end);

    	if( ! keep_strand )
    	{
    		if( r.isPlusStrand() )
    		{
    			start = std::min(s1,s2);
    			end = std::max(e1,e2);
    		}
    		else
    		{
    			start = std::max(e1,e2);
    			end = std::min(s1,s2);
    		}
    	}
    	else
    	{
    		if( isPlusStrand() )
    		{
    			start = std::min(s1,s2);
    			end = std::max(e1,e2);
    		}
    		else
    		{
    			start = std::max(e1,e2);
    			end = std::min(s1,s2);
    		}
    	}
    };

    friend bool operator<( const Range& r1, const Range& r2 )
    {
    	if(std::min(r1.start,r1.end)<std::min(r2.start,r2.end))
    		return true;
    	else
    		if(std::min(r1.start,r1.end)==std::min(r2.start,r2.end)
    		&& std::max(r1.start,r1.end)<std::max(r2.start,r2.end))
    			return true;
    	return false;
    };

    friend bool operator>( const Range& r1, const Range& r2 )
    {
    	if(std::max(r1.start,r1.end)>std::max(r2.start,r2.end))
    		return true;
    	else
    		if(std::max(r1.start,r1.end)==std::max(r2.start,r2.end)
    		&& std::min(r1.start,r1.end)>std::min(r2.start,r2.end))
    			return true;
    	return false;
    };

    friend bool operator==( const Range& r1, const Range& r2 )
    {
    	//	return r1.overlap(r2);
    	if(r1.start==r2.start && r1.end==r2.end)
    		return true;
    	return false;
    };

};
const Range emptyRange;

#endif /* RANGESEQ_H */
