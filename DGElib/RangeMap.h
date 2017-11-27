/***
 *
 * RangeMap.h
 *
 ***/

#ifndef RANGEMAP_H
#define RANGEMAP_H

#include <stdlib.h>
#include <fstream>
#include <map>
#include <list>

#include <SDGString.h>
#include <SDGFastaIstream.h>
#include <SDGFastaOstream.h>

//#include "SkipList.h"
#include "RangeSeq.h"
#include "SDGBioSeqDB.h"
//#include "Galign.h"

class RangeMap
: public std::map<std::string,std::list<RangeSeq> >
{
  unsigned countRange;
 public:


  RangeMap():
    std::map<std::string,std::list<RangeSeq> >(),countRange(0)
    {    };

  void add(RangeSeq r)
    {
      (*this)[r.getChr()].push_back(r);
      countRange++;
    };

  unsigned getCountRange(void){return countRange;};
  void load(const SDGString& filename);
  void save(const SDGString& filename);
  void saveSet(const SDGString& filename);
  void view(void);
  unsigned size(void);
  void cut(unsigned int length, unsigned int over);
  void sort(void);
  void diff( RangeMap& m, int verbose=0 );
  void selectInclude(RangeMap& m,RangeMap& mapout);
  void selectOverlap(RangeMap& m,RangeMap& mapout);
  void merge(void);
  void extend(unsigned size_flank);
  void selectSrcSeq(const SDGString& outfname,const SDGBioSeqDB& db);
  void writeSeq(const SDGString& outfname,const SDGBioSeqDB& db) const;
  void writeCutSeq( const SDGString& outfname, const SDGBioSeqDB& db, int verbose=0 );
  void writeFlank53Seq(const SDGString& outfname,const SDGBioSeqDB& db,
		       unsigned len=100);
  void writeFlank5Seq(const SDGString& outfname,const SDGBioSeqDB& db,
		      unsigned len=100);
  void writeFlank3Seq(const SDGString& outfname,const SDGBioSeqDB& db,
		      unsigned len=100);
  void entropfilt(const SDGString& db, double thres) ;
};
#endif
