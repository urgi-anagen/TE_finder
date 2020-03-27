/**
 * \file RangePairSet.h
 * \brief Header file for the class RangePairSet
 */

#ifndef RANGEPAIRSET_H
#define RANGEPAIRSET_H

#include <iostream>
#include <SDGString.h>
#include <utility>
#include <list>
#include <vector>
#include <map>
#include "RangeAlignSet.h"
#include "RangePair.h"

/**
 * \class RangePairSet
 * \brief Implement a chain of matches
 *
 * This class handles a list of RangePair objects.
 */
class RangePairSet: public RangePair
{
	protected:
		std::list<RangePair> path;

	public:

		RangePairSet(void){};

		RangePairSet(const std::list<RangePair>& rp_list)
		{
			setRpsFromRpList(rp_list);
		};

		RangePairSet(const RangePair& rp) : RangePair(rp)
		{
			path.push_back(rp);
		};

		RangeAlignSet getRangeAlignSetQ(void) const
		{
			std::list<Range> rs;
			for(std::list<RangePair>::const_iterator r=path.begin();
			r!=path.end();r++)
				rs.push_back(Range(r->getRangeQ()));
			RangeAlignSet ras(getRangeQ(),rs);
			ras.merge_set();
			return ras;
		};

		RangeAlignSet getRangeAlignSetS(void) const
		{
			std::list<Range> rs;
			for(std::list<RangePair>::const_iterator r=path.begin();
			r!=path.end();r++)
				rs.push_back(Range(r->getRangeS()));
			RangeAlignSet ras(getRangeS(),rs);
			ras.merge_set();
			return ras;
		};

		void clear(void)
		{
		  RangePair::clear();
		  path.clear();
		}

   /*
    * \fn long getNumQuery( void )
    * \brief return the identifier of the query of the first RangePair
    */
   long getNumQuery( void ) const
   {
  	 std::list<RangePair>::const_iterator i=path.begin();
  	 return i->getRangeQ().getNumChr();
   }

   /*
    * \fn long getNumSubject( void )
    * \brief return the identifier of the subject of the first RangePair
    */
   long getNumSubject( void ) const
   {
  	 std::list<RangePair>::const_iterator i=path.begin();
  	 return i->getRangeS().getNumChr();
   }

   /*
    * \fn unsigned getNbRangePairs( void )
    * \brief return the number of RangePair instance in the list
    */
   unsigned getNbRangePairs( void )
   {
   	return path.size();
   }

  std::list<RangePair>::iterator begin(){ return path.begin(); };
  std::list<RangePair>::iterator end(){ return path.end(); };

  std::list<RangePair>::const_iterator begin() const { return path.begin(); };
  std::list<RangePair>::const_iterator end() const { return path.end(); };

  void view( void ) const;
  void viewWithLabel( void );

  void computeScoreWithDynaProg( double mism=0.0, double gapo_p=0.0, double gape_p=0.0 );
  void computeScoreWithLengthAndId();

  void setRpsFromRpList(const std::list<RangePair> rp_list);
  void updateQueryFromPathList(void);
  void addPath(const RangePair& rp);

  void setQSName(std::string query_name, std::string subject_name);

  void write( std::ostream& out, unsigned id,
  		const std::string& nameQ, const std::map<long,std::string>& nameS ) const;
  void write( std::ostream& out, unsigned id,
                const std::string& nameQ, const std::string& nameS ) const;
  void writeGFF3( std::ostream& out, unsigned id,
                    const std::string& nameQ, const std::string& nameS, const std::string& source ) const;
  void writeGFF3( std::ostream& out, unsigned id,
  		const std::string& nameQ, const std::map<long,std::string>& nameS, const std::string& source ) const;
  void writeBED( std::ostream& out, const std::string& nameQ,  const std::map<long,std::string>& nameS, const std::string& color) const;
  void writeRpsAttr(std::ostream& out, unsigned id, const std::string& nameQ, const std::string& nameS) const;
  void writeRpsAttr(std::ostream& out, unsigned id, const std::string& nameQ, const std::map<long,std::string>& nameS) const;

  bool diffQ( const RangePairSet& r );
  void mergeQ(RangePairSet& rpsOther);
  bool split( const RangePairSet& r, std::list<RangePairSet>& lrp_out );
  bool inserted( const RangePairSet& r );
  unsigned overlapQ_length( const RangePairSet& r ) const;


  void setPathDirectly(const std::list<RangePair>& rp_list);
  const std::list<RangePair> getPath(){return path;};



  void readtxt(std::istream& in);

  friend bool operator==( const RangePairSet &rp1, const RangePairSet &rp2 );
  friend bool operator!=( const RangePairSet &rp1, const RangePairSet &rp2 );
	
  static bool compareCoordsOnQuery(const RangePairSet &rps1, const RangePairSet &rps2);

  void orientSubjects(void);
  private: 
	void cleanConflictsOnOverlappingQuery(RangePairSet &rpsOther);
	bool doOrientOnPlusStrand(void);
};
#endif
