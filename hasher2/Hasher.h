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
	  		unsigned connect_dist, unsigned kmer_size, unsigned min_frag_size,  unsigned verbose);


	SDGBioSeqDB subject_db;
	std::list< RangePair > frag;



 public:

  Hasher(unsigned w=10, unsigned msk=100, unsigned bw=2, unsigned wd=1, unsigned fd=1, unsigned minsize=20,unsigned step=1):
			  Duster(w,msk,bw,wd,fd,minsize,step)
    {};
  void load(const SDGString& filenameS, unsigned kmer_size, unsigned kmask, unsigned bkmer_size, unsigned mkmer_size, double count_cutoff, double diversity_cutoff,
		  unsigned min_count, bool & valid_idx_file)
	  {
		  Duster::load(filenameS,kmer_size, kmask, bkmer_size,mkmer_size , count_cutoff, diversity_cutoff, min_count,valid_idx_file, true);
		  subject_db.load(filenameS);
	  };
  void search(const SDGBioSeq& sequence, unsigned start, unsigned end, unsigned numseq, unsigned connect_dist, unsigned min_frag_size,
		  bool repeat, unsigned verbose);
  void write_align(const SDGBioSeq& sequence, std::ostream& out);

};


#endif






