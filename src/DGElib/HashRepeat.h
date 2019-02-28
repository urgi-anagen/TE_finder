#ifndef HASHREPEAT_H
#define HASHREPEAT_H

#include <limits>
#include <SDGString.h>
#include <SDGMemBioSeq.h>
#include <vector>
#include <list>
#include <unistd.h>
//#include "gnuplot_i.h"
#include "Range.h"
#include "FragAlign.h"
#include "FastExtAlign.h"
//#include "RangeGroups.h"

class HashRepeat
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

   unsigned operator()(std::string::const_reverse_iterator p)
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


  struct Diag
    {
      int diag;
      unsigned wpos;

      Diag(int d=0, unsigned p=0): diag(d), wpos(p) {};

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
  unsigned max_key;
  std::vector< std::vector<unsigned>::iterator > hash2wpos, hash_ptr;
  std::vector<unsigned> word_count;
  std::vector<unsigned> word_pos;

  std::vector< Diag > diag_map, diag_map_comp;

  std::vector< SDGString > queryName;
  unsigned word_size, nb_word, extend_len;
  unsigned nbseqQ;

  SDGBioSeq sequence, sequence_comp;

  std::vector< std::pair<unsigned,RangePair> > frag;
  std::list<RangePair> rp_list;
  std::list<RangePairSet> path;

  void prepare(double cutoff=0.01);
  void matchWords(void);
  void diagSearch(void);
  void fragAlign(double match,double mism, double gopen,
		 double gext, unsigned over);
  void extend(void);

  //void reputer(unsigned errors=1);

  // gnuplot_ctrl* h ;
  SDGString datafile1,datafile2,datafile3;

 public:

  HashRepeat(unsigned w=6, unsigned ext_len=20):
    hseq(w),
    max_key(0),
    word_size(w), 
    nb_word(0),
    extend_len(ext_len),
    nbseqQ(0) 
    // h(NULL)
    {
      if(w>=16) 
	throw SDGException(NULL,"HashSearch: Word size must be < 16 !!");
      max_key=(unsigned)pow(4,word_size);
      hash_ptr.resize(max_key+1);
      hash2wpos.resize(max_key+1);
      word_count.resize(max_key);
    };


 public:
  ~HashRepeat(void)
    {
/*       if(h) */
/* 	{ */
/* 	  gnuplot_close(h); */

/* 	  SDGString rm_cmd="rm "+datafile1+" "+datafile2+" "+datafile3; */
/* 	  system(rm_cmd); */
/* 	} */
    };
  void search(SDGBioSeq& seq);
  //void search_reputer(SDGBioSeq& seq,unsigned errors);
  void write(std::ostream& out=std::cout);
  void print(std::ostream& out=std::cout);
  void writeMap(std::ostream& out=std::cout);
  unsigned writeSet(std::ostream& out=std::cout,unsigned count=0); 
  std::list<RangePairSet> getRangePairList(void) { return path;};
/*   void initPlot(SDGString filename) */
/*     { */
/*       h = gnuplot_init() ; */
/*       gnuplot_cmd(h,"set terminal postscript landscape color"); */

/*       std::ostringstream cmd; */
/*       cmd<<"set output \" "<<filename<<"\" "<<std::endl; */
/*       gnuplot_cmd(h,(char*)cmd.str().c_str()); */
      
/*       std::ostringstream tmpfilename1; */
/*       tmpfilename1<<"gnuplot_"<<getpid()<<"_1.tmp"; */
/*       datafile1=tmpfilename1.str(); */

/*       std::ostringstream tmpfilename2; */
/*       tmpfilename2<<"gnuplot_"<<getpid()<<"_2.tmp"; */
/*       datafile2=tmpfilename2.str(); */

/*       std::ostringstream tmpfilename3; */
/*       tmpfilename3<<"gnuplot_"<<getpid()<<"_3.tmp"; */
/*       datafile3=tmpfilename3.str(); */
/*     }; */
/*   void plot(const SDGString& title); */
//  double compress(unsigned min_size=10, double cover=0.95);

  struct comp {
    int operator () (const std::pair<unsigned,RangePair> & a, const std::pair<unsigned,RangePair> & b)
    {return a.first < b.first;};
  };

};


#endif






