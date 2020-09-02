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
#include "BLRBioSeqDB.h"

#include "HashDNASeq.h"

class Test_Hasher;


class Hasher : public HashDNASeq
{
	friend class Test_Hasher;

	void matchKmers(const SDGBioSeq& sequence,
			    unsigned start, unsigned end, unsigned numseq, bool repeat,
				std::vector< std::list<Diag> >& diag_map, std::vector< std::list<Diag> >& diag_map_comp);

	void diagSearch(const SDGBioSeq& sequence, unsigned numseqQ, std::vector< std::list<Diag> >& diag_map,
	        std::vector< std::list<Diag> >& diag_map_comp,
	        unsigned connect_dist, unsigned kmer_size, unsigned min_frag_size,
	        std::list< RangePair >& frag, unsigned verbose);

    void fragMerge(std::list< RangePair >& frag);

    static RangePair rangePairFactory(const unsigned numseqQ, unsigned int qstart, unsigned int qend,
                               const unsigned numseqS, unsigned int sstart, unsigned int send,
                               unsigned int score, unsigned step_q, unsigned id) {
        double kmer_density = (double(score) / ((send - sstart + 1) / step_q));
        unsigned kmer_score = (double) (send - sstart + 1) * kmer_density;
        if(numseqS>100000){std::cout<<"error subject num:"<<numseqS<<std::endl;}
        RangePair rp(numseqQ,qstart,qend,numseqS,sstart,send,kmer_score,0,kmer_density, id);
        return rp;
    }


	std::vector<std::string> subject_names;



 public:

  Hasher(unsigned w=10, unsigned msk=100, unsigned bw=2, unsigned wd=1, unsigned fd=1, unsigned minsize=20,unsigned step=1):
          HashDNASeq(w, msk, bw, wd, fd, minsize, step)
    {};
  void load(const SDGString& filenameS, unsigned kmer_size, unsigned kmask, unsigned bkmer_size, unsigned mkmer_size, double count_cutoff, double diversity_cutoff,
		  unsigned min_count, bool & valid_idx_file, bool first_iter)
	  {
		  HashDNASeq::load(filenameS, kmer_size, kmask, bkmer_size, mkmer_size , count_cutoff, diversity_cutoff, min_count, valid_idx_file, first_iter);
		  std::string line;
		   std::ifstream myfile (filenameS);
          subject_names.clear();
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

  void search(const SDGBioSeq& sequence, unsigned start, unsigned end, unsigned numseq, unsigned connect_dist,
              unsigned min_frag_size, bool repeat, std::list< RangePair >& fmerged, unsigned verbose);
  unsigned fragStat(const std::list< RangePair >& frag, double quantile, unsigned& coverage);
  void fragLenFilter(std::list< RangePair >& frag, unsigned min_len);
  void fragAlignWrite(std::list< RangePair >& frag, const SDGString& qfilename, const SDGString& sfilename, std::ostream& out);
  void fragSeqWrite(const std::list< RangePair >& frag, const SDGString& fasta_filename, SDGFastaOstream& out);
};


#endif






