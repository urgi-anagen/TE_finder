#ifndef HASHER_H
#define HASHER_H

#include <SDGMemBioSeq.h>
#include <SDGSubBioSeq.h>
#include <SDGFastaIstream.h>
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
				std::vector< std::list<Diag> >& diag_map, std::vector< std::list<Diag> >& diag_map_comp);
	void diagSearch(const SDGBioSeq& sequence, std::vector< std::list<Diag> >& diag_map, std::vector< std::list<Diag> >& diag_map_comp,
	  		unsigned connect_dist, unsigned kmer_size, unsigned min_frag_size, std::ostream& out, unsigned verbose);


	std::vector<std::string> subject_names;



 public:

  Hasher(unsigned w=10, unsigned msk=100, unsigned bw=2, unsigned wd=1, unsigned fd=1, unsigned minsize=20,unsigned step=1):
			  Duster(w,msk,bw,wd,fd,minsize,step)
    {};
  void load(const SDGString& filenameS, unsigned kmer_size, unsigned kmask, unsigned bkmer_size, unsigned mkmer_size, double count_cutoff, double diversity_cutoff,
		  unsigned min_count, bool & valid_idx_file)
	  {
		  Duster::load(filenameS,kmer_size, kmask, bkmer_size,mkmer_size , count_cutoff, diversity_cutoff, min_count,valid_idx_file, true);
		  std::string line;
		   std::ifstream myfile (filenameS);
		   if (myfile.is_open())
		    {
		      while ( getline (myfile,line) )
		      {
		    	  if(line[0]=='>')
					{
					  subject_names.push_back(line.substr(1));
					}
		      }
		      myfile.close();
		    }
	  };
  void search(const SDGBioSeq& sequence, unsigned start, unsigned end, unsigned numseq, unsigned connect_dist, unsigned min_frag_size,
		  bool repeat, std::ostream& out, unsigned verbose);

};


#endif






