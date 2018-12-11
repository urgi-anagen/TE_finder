/**
 * \file RangeAlignSet.h
 * \brief Header file for the class RangeAlignSet
 */

#ifndef RANGEALIGNSET_H
#define RANGEALIGNSET_H

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <string>
#include <list>
#include <SDGError.h>
#include "RangeAlign.h"

/**
 * \class RangeAlignSet
 * \brief Implement a chain of fragments
 *
 * This class handles a list of Range objects.
 */
class RangeAlignSet : public RangeAlign
{
  std::list<Range> range_set;
  unsigned length_set;
  bool included;

	protected:

  void calcLengthSet( void )
  {
  	length_set=0;
  	for(std::list<Range>::iterator i=range_set.begin();
  	i!=range_set.end();i++)
  		length_set+=i->getLength();
  };

	public:
		RangeAlignSet(void):range_set(),length_set(0),included(false){};
		RangeAlignSet(const RangeAlign& r,std::list<Range> rs)
		: RangeAlign(r), range_set(rs), included(false)
		{
			if(range_set.empty())
				range_set.push_back(Range(*this));
			calcLengthSet();
		}
		RangeAlignSet(const RangeAlign& r)
		: RangeAlign(r), range_set(), included(false)
		{
			range_set.push_back(Range(*this));
			calcLengthSet();
		}

  ~RangeAlignSet(void) {};

  void reset( void )
  {
  	range_set.clear();
  	length_set = 0;
  	included = false;
  	RangeAlign::set();
  };

  bool isEmpty( void )
  {
  	return( range_set.size() == 0 );
  }

  /*
   * \fn friend bool operator==( const RangeAlignSet& ras1, const RangeAlignSet& ras2 )
   * \brief overload equal operator
   * \note need to sort both range_set first!
   */
  friend bool operator==( const RangeAlignSet& ras1, const RangeAlignSet& ras2 )
  {
	  if (ras1.range_set.size() != ras2.range_set.size())
		  return false;
	  std::list<Range>::const_iterator it1 = ras1.range_set.begin();
	  std::list<Range>::const_iterator it2 = ras2.range_set.begin();
	  while (*it1 == *it2 && it1!=ras1.range_set.end() && it2!=ras2.range_set.end()) {
		  it1++;
		  it2++;
	  }
	  if (it1 != ras1.range_set.end() || it2 != ras2.range_set.end())
		  return false;
	  if (ras1.included != ras2.included)
		  return false;
	  return true;
  };

  /*
   * \fn unsigned long getStart( void ) const
   * \brief return the start coordinate
   */
  unsigned long getStart( void ) const
  {
  	unsigned long start;
  	std::list<Range>::const_iterator it = range_set.begin();
  	start = it->getStart();
  	bool plusStrand = isPlusStrand();
  	while( it != range_set.end() )
  	{
  		if( plusStrand )
  			start = std::min( start, it->getStart() );
  		else
  			start = std::max( start, it->getStart() );
  		++ it;
  	}
  	return( start );
  }

  /*
   * \fn unsigned long getEnd( void ) const
   * \brief return the end coordinate
   */
  unsigned long getEnd( void ) const
  {
  	unsigned long end;
  	std::list<Range>::const_iterator it = range_set.begin();
  	end = it->getEnd();
  	bool plusStrand = isPlusStrand();
  	while( it != range_set.end() )
  	{
  		if( plusStrand )
  			end = std::max( end, it->getEnd() );
  		else
  			end = std::min( end, it->getEnd() );
  		++ it;
  	}
  	return( end );
  }

  /*
   * \fn bool isPlusStrand( void )
   * \brief return true if all Range instances are on plus strand, false if all are on minus strand, and throw an exception otherwise
   */
  bool isPlusStrand( void ) const
  {
  	unsigned nbPlusStrand = 0;
  	unsigned nbMinusStrand = 0;
  	for( std::list<Range>::const_iterator it=range_set.begin(); it!=range_set.end(); it++ )
  	{
  		if( it->isPlusStrand() )
  			nbPlusStrand ++;
  		else
  			nbMinusStrand ++;
  	}
  	if( nbPlusStrand == range_set.size() )
  		return true;
  	else if( nbMinusStrand == range_set.size() )
  		return false;
  	else
  		throw SDGException( this, "RangeAlignSet::isPlusStrand: ranges on different strands!" );
  }

  /*
   * \fn void sortUsingStrand( void )
   * \brief sort range_set with Range::less if plus strand and Range::greater if minus strand
   */
  void sortUsingStrand( void )
  {
  	if( isPlusStrand() )
  		range_set.sort( Range::less );
  	else if( ! isPlusStrand() )
  		range_set.sort( Range::greater );
  }

  void setIncluded( bool b ) { included = b; };
  bool getIncluded( void ) { return included; };

  /*
   * \fn unsigned getLengthSet( void )
   * \brief return attribute 'length_set'
   */
  unsigned getLengthSet( void )
  {
  	if(length_set==0)
  		calcLengthSet();
  	return length_set;
  };

  const std::list<Range>& getRangeSet(void) const {return range_set;};
  unsigned getRangeSetSize(void) const {return range_set.size();};

  /*
   * \fn void reverse( void )
   * \brief reverse the instance
   */
  void reverse( void )
  {
  	this->Range::reverse();
  	for( std::list<Range>::iterator i=range_set.begin(); i!=range_set.end(); i++ )
  		i->reverse();
  	sortUsingStrand();
  };

  /*
   * \fn unsigned overlap_length(const RangeAlignSet& r) const
   * \brief return the length of the overlap between the instance and r
   */
  unsigned overlap_length( const RangeAlignSet& r ) const
  {
  	unsigned overlap = 0;
  	if( getNumChr() != r.getNumChr() )
  		return false;
  	for( std::list<Range>::const_iterator i=range_set.begin(); i!=range_set.end(); i++ )
  	{
  		bool found = false;
  		for( std::list<Range>::const_iterator j=r.range_set.begin(); j!=r.range_set.end(); j++ )
  		{
  			if( i->overlap( *j ) )
  			{
  				found = true;
  				if( i->isIncluded(*j) )
  					overlap += j->getLength();
  				else if( i->isContained( *j ) )
  					overlap += i->getLength();
  				else if( i->getMin() > j->getMin() )
  					overlap += j->getMax() - i->getMin();
  				else
  					overlap += i->getMax() - j->getMin();
  			}
  			else if( found ) break;
  		}
  	}
  	return overlap;
  };

  /*
   * \fn void view( void )
   * \brief send the data of the instance to stdout
   */
  void view( void )
  {
  	std::cout<<getNumChr()
  			<<"\t"<<getStart()
  			<<"\t"<<getEnd()
  			<<"\trange set:";
  	for( std::list<Range>::const_iterator i=range_set.begin(); i!=range_set.end(); i++ )
  		std::cout<<i->getStart()
  		<<".."
  		<<i->getEnd()
  		<<",";
  	std::cout<<std::endl;
  };

  /*
   * \fn bool isStrictlyContained( RangeAlignSet r )
   * \brief return true if each Range instance of the object is strictly contained in at least one Range instance of r
   * \param r RangeAlignSet object
   */
  bool isStrictlyContained( const RangeAlignSet& r )
  {
	  unsigned nbContainedRanges = 0;
	  for( std::list<Range>::const_iterator i=range_set.begin(); i!=range_set.end(); i++ )
	  {
	  	for( std::list<Range>::const_iterator j=r.range_set.begin(); j!=r.range_set.end(); j++ )
	  	{
	  		if( i->isStrictlyContained( *j ) )
	  		{
	  			nbContainedRanges ++;
	  			continue;
	  		}
	  	}
	  }
	  if( nbContainedRanges == range_set.size() )
	  	return true;
	  else
		  return false;
  }

  /*
   * \fn bool isContained( RangeAlignSet r )
   * \brief return true if each Range instance of the object is contained in at least one Range instance of r
   * \param r RangeAlignSet object
   */
  bool isContained( const RangeAlignSet& r )
  {
	  unsigned nbContainedRanges = 0;
	  for( std::list<Range>::const_iterator i=range_set.begin(); i!=range_set.end(); i++ )
	  {
	  	for( std::list<Range>::const_iterator j=r.range_set.begin(); j!=r.range_set.end(); j++ )
	  	{
	  		if( i->isContained( *j ) )
	  		{
	  			nbContainedRanges ++;
	  			continue;
	  		}
	  	}
	  }
	  if( nbContainedRanges == range_set.size() )
	  	return true;
	  else
		  return false;
  }

  /*
   * \fn bool doesItContain( const RangeAlignSet& r )
   * \brief return true if each Range instance of the object contains at least one Range instance of r
   */
  bool doesItContain( const RangeAlignSet& r )
  {
  	unsigned nbRangesItContains = 0;
  	for( std::list<Range>::const_iterator i=range_set.begin(); i!=range_set.end(); i++ )
  	{
  		for( std::list<Range>::const_iterator j=r.range_set.begin(); j!=r.range_set.end(); j++ )
  		{
  			if( i->isIncluded( *j ) )
  			{
  				nbRangesItContains ++;
  				continue;
  			}
  		}
  	}
  	if( nbRangesItContains == range_set.size() )
  		return true;
  	else
  		return false;
  }

  /*
   * \fn void merge_set( void )
   * \brief merge the Range instances recorded in 'range_set'
   */
  void merge_set( void )
  {
  	sortUsingStrand();
  	std::list<Range>::iterator i=range_set.begin();
  	std::list<Range>::iterator j=range_set.begin();j++;
  	while( i!=range_set.end() && j!=range_set.end() )
  	{
  		while( i!=range_set.end() && j!=range_set.end() && i->overlap(*j) )
  		{
  			if( i != j )
  			{
  				i->merge( *j );
  				j = range_set.erase( j );
  			}
  			else
  				j++;
  		}
  		if( i!=range_set.end() && j!=range_set.end() )
  		{
  			i++;
  			j++;
  		}
  	}
  	calcLengthSet();
  };

  /*
   * \fn void merge( RangeAlignSet r )
   * \brief merge self with a given RangeAlignSet object
   * \param r RangeAlignSet object
   * \note if self and r have different 'included' and self doesn't contain r, 'included' of self is put to 'false'
   * \note if self has 'included=false', r has 'included=true' and self is contained in r, 'included' of self is put to 'true'
   */
  void merge( RangeAlignSet r )
  {
  	if( ( getNumChr() != r.getNumChr() ) || ( getNameSeq() != r.getNameSeq() ) )
  		return;

  	if( isPlusStrand() != r.isPlusStrand() )
  		r.reverse();

  	if( included == r.included )
  	{
  		range_set.insert( range_set.end(), r.range_set.begin(), r.range_set.end() );
  		merge_set();
  	}

  	else  // one has 'included=true' and not the other
  	{
  		sortUsingStrand();
  		r.sortUsingStrand();
  		if( range_set == r.range_set )
  			included = true;
  		else
  		{
  			if( ! doesItContain( r ) )
  				included = false;
  			if( ! included && isContained( r ) )
  				included = true;
  			range_set.insert( range_set.end(), r.range_set.begin(), r.range_set.end() );
  			merge_set();
  		}
  	}
  };

  /*
   * \fn hasSameRanges( RangeAlignSet& ras2 )
   * \brief return true if ras2 has the same Ranges as self, false otherwise
   */
  bool hasSameRanges( RangeAlignSet& ras2 )
  {
	  sortUsingStrand();
	  ras2.sortUsingStrand();
	  if( range_set.size() != ras2.range_set.size() )
		  return false;
	  std::list<Range>::const_iterator it1 = range_set.begin();
	  std::list<Range>::const_iterator it2 = ras2.range_set.begin();
	  while( *it1 == *it2 )
	  {
		  it1 ++;
		  it2 ++;
	  }
	  if( it1 != range_set.end() )
		  return false;
	  return true;
  };

};

#endif /* RANGEALIGNSET_H */
