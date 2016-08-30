#ifndef HASHSEARCH_H
#define HASHSEARCH_H

#include <SDGMemBioSeq.h>
#include <SDGFastaIstream.h>
#include <SDGBioSeqDB.h>
#include <fstream>
#include <string>
#include <cmath>
#include <list>

#include "RangePair.h"
#include "FragAlign.h"
#include "FastExtAlign.h"

class Hasher
{
  struct HashDNASeq // Compute hashing value of a word
  {
    unsigned nuc_val[256];
    unsigned word_size;
    unsigned maskw;
      
    HashDNASeq(unsigned w): word_size(w)
    {  
      if(w>16) 
	throw SDGException(NULL,"HashDNASeq: Word size must be <= 16 !!");

      maskw = 1;
      maskw<<=((word_size*2-1)+1);
      maskw--;
	
      for(unsigned i=0;i<256;i++)
    	  nuc_val[i]=0;
      
      nuc_val[(unsigned)'A']=0;
      nuc_val[(unsigned)'C']=1;
      nuc_val[(unsigned)'G']=2;
      nuc_val[(unsigned)'T']=3; 
      
      nuc_val[(unsigned)'a']=0;
      nuc_val[(unsigned)'c']=1;
      nuc_val[(unsigned)'g']=2;
      nuc_val[(unsigned)'t']=3; 
    };

   unsigned operator()(std::string::const_iterator p)
    {
      unsigned h=0;
      for(unsigned i=0;i<word_size;i++)
		{
		  h<<=2;
		  h|=nuc_val[(unsigned)*p++];
		}
      return h&maskw;
    };
  };

  struct WordSpos //Store a word position in a sequence
    {
      unsigned pos;
      unsigned numSeq;

      WordSpos(unsigned p=0, unsigned n=0): pos(p), numSeq(n) {};

      friend int operator< (const WordSpos& w1, const WordSpos& w2)
      {
		if(w1.numSeq<w2.numSeq) return 1;
		if(w1.numSeq>w2.numSeq) return 0;
		if(w1.numSeq==w2.numSeq && w1.pos<w2.pos) return 1;
		return 0;
      };

      friend int operator== (const WordSpos& w1, const WordSpos& w2)
      {
		if(w1.numSeq!=w2.numSeq) return 0;
		if(w1.pos!=w2.pos) return 0;
		return 1;
      };      
      
      friend std::ostream& operator<<(std::ostream& os, const WordSpos& w)
      {
		os<<"("<<w.numSeq<<","<<w.pos<<")";
		return os;
      };

  };

  struct Diag // Store the diagonal of a word match
    {
      int diag;
      WordSpos wpos;

      Diag(int d=0, unsigned p=0, unsigned n=0): diag(d), wpos(p,n) {};

      friend std::ostream& operator<<(std::ostream& os, const Diag& d)
      {
		os<<"["<<d.wpos<<","<<d.diag<<"]";
		return os;
      };

      friend int operator< (const Diag& d1, const Diag& d2)
      {
		if(d1.diag<d2.diag) return 1;
		if(d1.diag==d2.diag && d1.wpos<d2.wpos) return 1;
		return 0;
      };
  };


  HashDNASeq hseq;
  std::vector< std::vector<WordSpos>::iterator > hash2wpos, hash_ptr;
  std::vector<unsigned> word_count;
  std::vector<WordSpos> word_pos;

  std::vector< Diag > diag_map, diag_map_comp;

  std::vector< SDGString > subjectName;
  unsigned word_size, step_q, step_s, wdist, max_key, nb_word, nfrag;
  double cutoff;
  unsigned nbseqQ,nbseqS;
  bool repeat_mode;

  SDGBioSeq sequence, comp_sequence;
  SDGBioSeqDB subject_db;

  std::vector< std::pair<unsigned,RangePair> > frag;

  struct Key : std::pair<long,long>
  {
      Key(long i,long j)
      { 
    	  first=i; second=j;
      };
  }; 
  typedef std::map<Key, std::list<RangePair> > MapAlign;
  typedef std::map<Key, std::list<RangePairSet> > MapPath;
  MapAlign map_align;
  MapPath map_path;

  bool read_idx(const SDGString& filename);
  void save_idx(const SDGString& filename);
  void hashSeqCount(const SDGBioSeq& seq);
  void hashSeqPos(const SDGBioSeq& seq);
  void matchWords(void);
  void matchWords_repeat(void);
  void diagSearch(void);

  int extend_len;

 public:

  Hasher(unsigned w=6,unsigned wd=1, unsigned n=100, double cut=0.0, bool rep_mode=false, int ext_len=-1):
    hseq(w),
    word_size(w), 
    step_q(w),
    step_s(1),
    wdist(wd),
    nb_word(0), 
    nfrag(n), 
    cutoff(cut),
    nbseqQ(0), 
    nbseqS(0),
    repeat_mode(rep_mode),
    extend_len(ext_len)
    {
      if(w>16) 
		throw SDGException(NULL,"Hasher: Word size must be <= 16 !!");
      max_key=(unsigned)pow(4,word_size);
      hash2wpos.resize(max_key+1);
      word_count.resize(max_key);
      if(extend_len==-1) extend_len=2*word_size;
    };

  void load(const SDGString& filenameS);
  void fragAlign(double match,double mism, double gopen,
		 double gext, unsigned over, bool join);

  void search(const SDGBioSeq& seq);
  void write(unsigned min_size=6,std::ostream& out=std::cout);
  void write_align(unsigned min_size,std::ostream& out);
  void print_frag(std::ostream& out=std::cout);
  void print(unsigned min_size=6,std::ostream& out=std::cout);
  void extend(unsigned min_size=6);

  struct comp {
    int operator () (const std::pair<unsigned,RangePair> & a, const std::pair<unsigned,RangePair> & b)
    {return a.first < b.first;};
  };

};


#endif






