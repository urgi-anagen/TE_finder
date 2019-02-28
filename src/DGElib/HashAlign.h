#ifndef HASHREPEAT_H
#define HASHREPEAT_H

#include <SDGMemBioSeq.h>
#include <string>
#include <list>
#include <map>
//#include "gnuplot_i.h"
#include "Diag.h"
#include "Range.h"
#include "FragAlign.h"

class HashAlign
{
 private:
  SDGBioSeq sequence1;
  SDGBioSeq sequence2;
  std::multimap<unsigned long,RangePair> sort_frag;
  std::list<RangePair> rp_list;
  std::list<RangePairSet> path;

  //gnuplot_ctrl* h ;
  SDGString datafile1,datafile2;

 public:
  //  HashAlign(void):h(NULL){};
  HashAlign(void){};
  //  HashAlign(SDGBioSeq& seq1,SDGBioSeq& seq2)
  //  :sequence1(seq1),sequence2(seq2),h(NULL)
  //  {};
  HashAlign(SDGBioSeq& seq1,SDGBioSeq& seq2)
    :sequence1(seq1),sequence2(seq2)
    {};
  ~HashAlign(void)
    {
/*       if(h) */
/* 	{ */
/* 	  gnuplot_close(h); */

/* 	  SDGString rm_cmd="rm "+datafile1+" "+datafile2; */
/* 	  system(rm_cmd); */
/* 	} */
    };

  void setSeq(SDGBioSeq seq1,SDGBioSeq seq2)
    { sequence1=seq1;sequence2=seq2; };

  void fragAlign(double match,double mism, double gopen,
		 double gext, unsigned over,unsigned nfrag=0);
  void search(unsigned word_size=6);
  void write(std::ostream& out=std::cout);
  void print(unsigned min_size=100, std::ostream& out=std::cout);
  int score(void)
    {
      int sc=0;
      for(std::list<RangePairSet>::iterator i=path.begin();
	  i!=path.end();i++)
	{
	  sc+=i->getScore();
	}
      return sc;
    }
  std::list<RangePairSet> getMatches(void){return path;};

  void load_matches(SDGString filename);
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
/*     }; */
/*   void plot(void); */
};


#endif

