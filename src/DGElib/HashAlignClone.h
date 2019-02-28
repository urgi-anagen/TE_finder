#ifndef HASHALIGNCLONE_H
#define HASHALIGNCLONE_H

#include <SDGMemBioSeq.h>
#include <string>
#include <list>
#include <map>
#include "DiagClone.h"
#include "Range.h"
#include "FragAlign.h"
#include "HashAlign.h"

class HashAlignClone : public HashAlign
{
 private:
  SDGBioSeq sequence1;
  SDGBioSeq sequence2;
  std::multimap<unsigned long,RangePair> sort_frag;
  std::list<RangePair> rp_list;
  std::list<RangePairSet> path;

  SDGString datafile1,datafile2;
  

 public:
  HashAlignClone(void);
  ~HashAlignClone(void);

  void setSeq(SDGBioSeq seq1,SDGBioSeq seq2)
    { sequence1=seq1;sequence2=seq2; };

  void search(unsigned word_size=6);

  void addToDiagMap(const Range& r, int key, std::map<int, std::list<Range> >& diag_map);
};
#endif

