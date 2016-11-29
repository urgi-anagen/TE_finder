#ifndef HASHSEARCH_H
#define HASHSEARCH_H

#include <SDGMemBioSeq.h>
#include <SDGFastaIstream.h>
#include <string>
#include <cmath>
#include <list>

#include "RangePair.h"
#include "FragAlign.h"

class HashSearch
{
  struct HashDNASeq
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

   unsigned operator()(std::string::const_iterator p, unsigned l, unsigned h=0)
    {
      for(unsigned i=0;i<l;i++)
	{
	  h<<=2;
	  h|=nuc_val[(unsigned)*p++];
	}
      return h&maskw;
    };
  };

  struct WordSpos
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

  struct Diag
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

  std::vector< SDGString > queryName, subjectName;
  unsigned word_size, max_key, nb_word, nfrag;
  double cutoff;
  unsigned nbseqQ,nbseqS;
  bool repeat_mode;

  SDGBioSeq sequence;

  std::vector< std::pair<unsigned,RangePair> > frag;
  std::list<RangePair> rp_list;
  std::list<RangePairSet> path;

  void hashSeqCount(SDGBioSeq& seq);
  void hashSeqPos(SDGBioSeq& seq);
  void matchWords(void);
  void matchWords_repeat(void);
  void diagSearch(void);

 public:

  HashSearch(unsigned w=6,unsigned n=100, double cut=0.0, bool rep_mode=false):
    hseq(w),
    word_size(w), 
    nb_word(0),
    nfrag(n), 
    cutoff(cut),
    nbseqQ(0), 
    nbseqS(0),
    repeat_mode(rep_mode)
    {
      if(w>16) 
	throw SDGException(NULL,"HashSearch: Word size must be <= 16 !!");
      max_key=(unsigned)pow(4,word_size);
      hash_ptr.resize(max_key+1);
      hash2wpos.resize(max_key+1);
      word_count.resize(max_key);
    };

  void load(SDGString& filenameS);
  void fragAlign(double match,double mism, double gopen,
		 double gext, unsigned over);

  void search(SDGBioSeq& seq);
  void write(unsigned min_size=6,std::ostream& out=std::cout, bool name=true);
  void print(unsigned min_size=6,std::ostream& out=std::cout, bool name=true);
  std::list<RangePair> getRangePairList(void) { return rp_list;};

  struct comp {
    int operator () (const std::pair<unsigned,RangePair> & a, const std::pair<unsigned,RangePair> & b)
    {return a.first < b.first;};
  };

};


#endif






