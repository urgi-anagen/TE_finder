#ifndef HASHDNASEQ_H
#define HASHDNASEQ_H

#include <SDGError.h>
#include <SDGMemBioSeq.h>
#include <SDGFastaIstream.h>
#include <fstream>
#include <string>
#include <cmath>
#include <list>
#include <tuple>
#include <unordered_map>
#include "HashFuncDNASeq.h"


class Test_HashDNASeq;


struct Info_kmer
{
	  unsigned count;
	  double expected_count;
	  double entropy;
	  double diversity;
	  double good_kmer;
	  unsigned hash_key;

	  Info_kmer(unsigned h=0, unsigned c=0, double ec=0, double e=0, double d=0, double g=0):
		  count(c),
		  expected_count(ec),
		  entropy(e),
		  diversity(d),
		  good_kmer(g),
		  hash_key(h)
		  {};

    bool operator< (const Info_kmer& k) const
    {
		if(count<k.count) return true;
		return false;
    };
};

struct Compare_Info_kmer {
  bool operator() (const Info_kmer& k1, const Info_kmer& k2)
  {
      return k1 < k2;
  };
};


class HashDNASeq
{
	friend class Test_HashDNASeq;

protected:

  struct KmerSpos //Store a kmer position in a sequence
    {
      unsigned pos;
      unsigned numSeq;

      KmerSpos(unsigned p=0, unsigned n=0): pos(p), numSeq(n) {
          if (n>1000000)
              throw Unsigned_Out_Of_Range("KmerSpos - sequence number setup out of range: n=",n);
          if(p>1000000000)
              throw Unsigned_Out_Of_Range("KmerSpos - sequence position setup out of range: p=",p);
          pos=p;
          numSeq=n;
      };

      friend int operator< (const KmerSpos& w1, const KmerSpos& w2)
      {
		if(w1.numSeq<w2.numSeq) return 1;
		if(w1.numSeq>w2.numSeq) return 0;
		if(w1.numSeq==w2.numSeq && w1.pos<w2.pos) return 1;
		return 0;
      };

      friend int operator== (const KmerSpos& w1, const KmerSpos& w2)
      {
		if(w1.numSeq!=w2.numSeq) return 0;
		if(w1.pos!=w2.pos) return 0;
		return 1;
      };      
      
      friend std::ostream& operator<<(std::ostream& os, const KmerSpos& w)
      {
		os<<"("<<w.numSeq<<","<<w.pos<<")";
		return os;
      };

  };

  struct Diag // Store the diagonal of a kmer match
    {
      long diag;
      KmerSpos wpos;

      Diag(long d = 0, unsigned p = 0, unsigned n = 0) {
          if (d > 100000000)
              throw Long_Out_Of_Range("Diag setup error: d=",d);
          diag = d;
          wpos = KmerSpos(p, n);
      };

      friend std::ostream& operator<<(std::ostream& os, const Diag& d)
      {
		os<<"["<<d.wpos<<","<<d.diag<<"]";
		return os;
      };

      friend int operator< (const Diag& d1, const Diag& d2)
      {
  		if(d1.wpos.numSeq<d2.wpos.numSeq) return 1;
  		if(d1.wpos.numSeq>d2.wpos.numSeq) return 0;
  		if(d1.wpos.numSeq==d2.wpos.numSeq)
  		{
  			if(d1.diag<d2.diag) return 1;
  			if(d1.diag==d2.diag)
  				if(d1.wpos.pos<d2.wpos.pos) return 1;
  		}
  		return 0;
      };
  };


  HashFuncDNASeq hseq,bhseq,mhseq,nhseq;
  //std::vector< std::vector<KmerSpos>::iterator > hash2wpos, hash_ptr;
  std::unordered_map < unsigned, std::vector<KmerSpos>::iterator > hash2wpos, hash_ptr;
  std::vector<KmerSpos> kmer_pos;



  std::vector< SDGString > subjectName;
  unsigned kmer_size, bkmer_size, mkmer_size, wdist, fdist, min_size, step_q, max_key, nfrag;
  unsigned nbseqQ,nbseqS;

  struct Key : std::pair<long,long>
  {
      Key(long i,long j)
      { 
    	  first=i; second=j;
      };
  }; 
  void kmer_counts(const SDGString& filenameS, unsigned kmer_size, unsigned bkmer_size, unsigned mkmer_size,
  		std::vector<unsigned>& wcount, unsigned& nb_kmer,
  		std::vector<unsigned>& bcount, unsigned& nb_bkmer,
  		std::vector<unsigned>& mcount, unsigned& nb_mkmer,
  		std::vector<unsigned>& ncount, unsigned& nb_nuc);
  void kmer_prob(unsigned wsize, unsigned bwsize, unsigned mwsize, unsigned mask, unsigned mask_hole_length,
		  const std::vector<unsigned>& wcount, unsigned nb_kmer,
		  const std::vector<unsigned>& bcount, unsigned nb_bkmer,
		  const std::vector<unsigned>& mcount, unsigned nb_mkmer,
		  const std::vector<unsigned>& ncount, unsigned nb_nuc,
		  std::list< Info_kmer >& list_infokmer);
  void kmer_count_percentiles(const std::list< Info_kmer >& list_infokmer, double cutoff_count,
		  Info_kmer& kmer_threshold);
  void kmer_entropy_percentiles(const std::list< Info_kmer >& list_infokmer, double cutoff_entropy,
		  Info_kmer& kmer_threshold);
  void kmer_diversity_percentiles(const std::list< Info_kmer >& list_infokmer, double cutoff_diversity,
		  Info_kmer& kmer_threshold);
  void kmer_goodkmer_percentiles(const std::list< Info_kmer >& list_infokmer);
  void kmer_filter(const std::list< Info_kmer >& list_infokmer, const Info_kmer& kmer_threshold, unsigned min_count,
  			std::vector<unsigned>& wcount, bool first_iter);
  void kmer_ssr_filter(unsigned wsize, std::vector<unsigned>& wcount);

  bool read_idx(const SDGString& filename, double count_cutoff, double diversity_cutoff, unsigned min_count, unsigned kmask, unsigned mask_hole_length);
  void save_idx(const SDGString& filename, double count_cutoff, double diversity_cutoff, unsigned min_count, unsigned kmask, unsigned mask_hole_length, const std::vector<unsigned>& wcount);
  unsigned hashSeqCount(const SDGBioSeq& seq, unsigned wsize, std::vector<unsigned>& wcount);
  unsigned hashSeqBackgroundCount(const SDGBioSeq& seq, unsigned wsize, std::vector<unsigned>& wcount);
  unsigned hashSeqModelCount(const SDGBioSeq& seq, unsigned wsize, std::vector<unsigned>& wcount);
  unsigned hashSeqNucCount(const SDGBioSeq& seq, std::vector<unsigned>& wcount);
  void hashSeqPos(const SDGBioSeq& seq, const std::vector<unsigned>& wcount);
  void matchKmers(const SDGBioSeq& sequence,
		    unsigned start, unsigned end, unsigned numseq, bool repeat,
			std::vector< Diag >& diag_map);
  void diagSearch(std::vector< Diag >& diag_map,
  		unsigned connect_dist, unsigned kmer_size,
  		std::vector< std::pair<unsigned,unsigned> >& frag);


 public:

  HashDNASeq(unsigned w=10, unsigned msk=100, unsigned mask_hole_length=1, unsigned bw=2, unsigned wd=1, unsigned fd=1, unsigned minsize=20, unsigned step=1):
    hseq(w,msk,mask_hole_length),
    bhseq(bw,0,0),
    mhseq(w/2,0,0),
    nhseq(1,0,0),
    kmer_size(w),
    bkmer_size(bw),
    mkmer_size(w/2),
    wdist(wd),
    fdist(fd),
    min_size(minsize),
    step_q(step),
    nbseqQ(0), 
    nbseqS(0)
    {
      if(hseq.getEffectiveKmerSize()>16)
		throw SDGException(NULL,"HashDNASeq: Kmer size must be <= 16 !!");
      max_key=(unsigned)pow(4,hseq.getEffectiveKmerSize());
    };

    virtual unsigned getEffectiveKmerSize() {return hseq.getEffectiveKmerSize();};

    virtual void load(const SDGString& filenameS, unsigned kmer_size, unsigned mask, unsigned mask_hole_length, unsigned bkmer_size, unsigned mkmer_size, double count_cutoff, double diversity_cutoff,
		  unsigned min_count,
		  bool & valid_idx_file, bool first_iter);
  void kmer_analysis(const SDGString& filenameS, unsigned kmer_size, unsigned mask, unsigned mask_hole_length, unsigned bkmer_size, unsigned mkmer_size,
		  double count_cutoff, double diversity_cutoff,
		  std::vector<unsigned>& wcount,
		  unsigned& nb_kmer,
		  std::list< Info_kmer >& list_infokmer, Info_kmer& kmer_threshold);

    virtual void search(const SDGBioSeq& seq, unsigned start, unsigned end,unsigned numseq, bool repeat,
		  std::vector< std::pair<unsigned,unsigned> >& fmerged);
};


#endif






