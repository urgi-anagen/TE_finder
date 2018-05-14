#ifndef HASHER_H
#define HASHER_H

#include <SDGMemBioSeq.h>
#include <SDGSubBioSeq.h>
#include <SDGFastaIstream.h>
#include <SDGBioSeqDB.h>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <map>
#include <tuple>
#include <unordered_map>

#include "RangePair.h"
#include "FragAlign.h"
#include "FastExtAlign.h"

#include "HashDNASeq.h"
#include "../duster/Duster.h"

class Test_Hasher;


class Hasher : public Duster
{
	friend class Test_Hasher;
	void matchKmers(const SDGBioSeq& sequence,
			    unsigned start, unsigned end, unsigned numseq, bool repeat,
				std::vector< Diag >& diag_map, std::vector< Diag >& diag_map_comp);
	void diagSearch(const SDGBioSeq& sequence, std::vector< Diag >& diag_map, std::vector< Diag >& diag_map_comp,
	  		unsigned connect_dist, unsigned kmer_size);

	int extend_len;

	SDGBioSeqDB subject_db;

	typedef std::map<Key, std::list<RangePair> > MapAlign;
	typedef std::map<Key, std::list<RangePairSet> > MapPath;
	MapAlign map_align;
	MapPath map_path;


 public:

  Hasher(unsigned w=10, unsigned msk=100, unsigned bw=2, unsigned wd=1, unsigned fd=1, unsigned minsize=20,unsigned step=1,
		  int ext_len=-1):
			  Duster(w,msk,bw,wd,fd,minsize,step),
			  extend_len(ext_len)
    {
      if(extend_len==-1) extend_len=2*kmer_size;
    };
  void load(const SDGString& filenameS, unsigned kmer_size, unsigned kmask, unsigned bkmer_size, unsigned mkmer_size, double count_cutoff, double diversity_cutoff,
		  unsigned min_count, bool & valid_idx_file)
	  {
		  Duster::load(filenameS,kmer_size, kmask, bkmer_size,mkmer_size , count_cutoff, diversity_cutoff, min_count,valid_idx_file, true);
		  subject_db.load(filenameS);
	  };
  void search(const SDGBioSeq& sequence, unsigned start, unsigned end, unsigned numseq, bool repeat);
  void fragAlign(double match,double mism, double gopen,
		 double gext, unsigned over, bool join);
  void write(const SDGBioSeq& sequence, unsigned min_size=6,std::ostream& out=std::cout);
  void write_align(const SDGBioSeq& sequence, unsigned min_size,std::ostream& out);
  void print_frag(const SDGBioSeq& sequence, std::ostream& out=std::cout);
  void print(const SDGBioSeq& sequence, unsigned min_size=6,std::ostream& out=std::cout);
  void extend(const SDGBioSeq& sequence, const SDGBioSeq& comp_sequence, unsigned min_size=6, unsigned verbose=0);

  struct comp {
    int operator () (const std::pair<unsigned,RangePair> & a, const std::pair<unsigned,RangePair> & b)
    {return a.first < b.first;};
  };
};


#endif






