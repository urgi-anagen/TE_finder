#include <unistd.h> //getopt
#include <fstream>
#include <list>
#include <RangePairSet.h>
#include <SDGFastaIstream.h>
#include <SDGMemBioSeq.h>

#include "Lalign.h"
#include "HashRepeat.h"

int match=10,mismh=16,gap_extend=32,gap_open=32;
unsigned word_len=6, min_LTRsize=100, min_LTRTEsize=1000, min_TIRsize=10, min_TIRTEsize=100,max_TIRTEsize=5000,max_LTRTEsize=20000, TSDsize=10,extend=20;
bool reputer=false;
    
void help(void)
{
  std::cerr<<"usage: TSDsearch"
	   <<" [<options>] <fasta sequence database>"<<std::endl
	   <<" options:"<<std::endl
	   <<"   -h, --help:\n\t this help"<<std::endl
	   <<"   -r, --reputer:\n\t reputer search"<<std::endl
	   <<"   -w, --word:\n\t word length, default: "<<word_len<<std::endl
	   <<"   -m, --match:\n\t match bonus (>0), default: "
	   <<match<<std::endl
	   <<"   -d, --mismatch:\n\t mismatch penalty (>0), default: "
	   <<mismh<<std::endl
	   <<"   -g, --gapopen:\n\t gap open penalty (>0), default: "
	   <<gap_open<<std::endl
	   <<"   -e, --gapextend:\n\t gap extend penalty (>0), default: "
	   <<"   -x, --extlen:\n\t length from extremities in %, default: "
	   <<extend<<std::endl;
};

int main(int argc, char* argv[])
{
  try{
    
    int c;
    while (1)
      {
// 	static struct option long_options[] =
// 	{
// 	  {"help",no_argument, 0, 'h'},
// 	  {"reputer",no_argument, 0, 'r'},
// 	  {"word",required_argument, 0, 'w'},
// 	  {"match",required_argument, 0, 'm'},
// 	  {"mismatch",required_argument, 0, 'd'},
// 	  {"gapopen",required_argument, 0, 'g'},
// 	  {"gapextend",required_argument, 0, 'e'},
// 	  {"extend",required_argument, 0, 'x'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;
	
// 	c = getopt_long (argc, argv, "hx:rw:m:d:g:e:",
// 			 long_options, &option_index);
	c = getopt (argc, argv, "hx:rw:m:d:g:e:");
	
	/* Detect the end of the options. */
	if (c == -1)
	  break;
	
	switch (c)
	  {
	  case 'h':
	    {
	      help();
	      return 0;
	    }
	  case 'x':
	    {
	      extend=atoi(optarg);
	      break;
	    }
	  case 'w':
	    {
	      word_len=atoi(optarg);
	      break;
	    }
	  case 'm':
	    {
	      match=atoi(optarg);
	      break;
	    }
	  case 'd':
	    {
	      mismh=atoi(optarg);
	      break;
	    }
	  case 'g':
	    {
	      gap_open=atoi(optarg);
	      break;
	    }
	  case 'e':
	    {
	      gap_extend=atoi(optarg);
	      break;
	    }
	  case 'r':
	    {
	      reputer=true;
	      break;
	    }
	  case '?':
	    help();
	    return 1;
	  default:
	    abort ();
	  }
      }

    if(argc==optind)
      {
	help();
	return 1;
      }

    SDGString filename=argv[optind];
    /* Print any remaining command line arguments (not options). */
    if (++optind < argc)
      {
	help();
	std::cout<<"non-option ARGV-elements: "<<std::endl;
	while (optind < argc)
	  std::cout<<argv[optind++]<<std::endl;
	return 1;
      }
     
    Lalign al(1);
    al.setMismatch(match,mismh);
    al.setGap(gap_open, gap_extend);

    std::list<RangePairSet> path;
    HashRepeat hrep(word_len,extend);

    SDGFastaIstream in(filename);
    SDGBioSeq s;

    SDGString filenameout1=filename.afterlast("/")+".itr";
    SDGFastaOstream out1(filenameout1);
    SDGString filenameout2=filename.afterlast("/")+".ltr";
    SDGFastaOstream out2(filenameout2);
    SDGString filenameout3=filename.afterlast("/")+".out";
    SDGFastaOstream out3(filenameout3);
    while(in)
      {
	in>>s;
	std::cout<<s.getDE()<<" len="<<s.length()<<std::endl;
//	if(reputer)
//	  hrep.search_reputer(s,1);
//	else
	  hrep.search(s);
	path=hrep.getRangePairList();
	for(std::list<RangePairSet>::iterator i=path.begin();
	    i!=path.end();i++)
	  {
	    if(i->getRangeQ().isPlusStrand() 
	       && i->getRangeQ().getLength()>=min_LTRsize
	       && i->getRangeQ().getEnd()-i->getRangeS().getStart()
	       >=min_LTRTEsize
	       && i->getRangeQ().getEnd()-i->getRangeS().getStart()
	       <=max_LTRTEsize)
	      {
		if(i->getRangeS().getStart()<TSDsize)
		  continue;
		SDGBioSeq s1=s.subseq(i->getRangeS().getStart()-TSDsize-1,TSDsize);
		SDGBioSeq s2=s.subseq(i->getRangeQ().getEnd()-1,TSDsize);
		if(s1.length()!=TSDsize || s2.length() !=TSDsize)
		  continue;
		al.setSeq(s1,s2);
		al.align();
		if(al.getScore()>0
		   && (unsigned)al.getEndSeq1()==TSDsize
		   && (unsigned)al.getStartSeq2()==1)
		  {
		    SDGString desc="LTR-TE";
		    desc+=" TE-size="+SDGString(std::fabs((int)(i->getRangeS().getStart())-i->getRangeQ().getEnd())+1);
		    desc+=" LTR-size="+SDGString(i->getRangeQ().getLength());
		    desc+=" TSD-size="+SDGString(s1.subseq(al.getStartSeq1()-1).length());
		    out3<<s.getDE()<<" len="<<s.length();
		    out3<<" ==>"<<desc<<std::endl;
		    out3<<" TSD5:"<<s1.subseq(al.getStartSeq1()-1).toString()
			     <<" "<<i->getRangeS().getStart()-TSDsize
		      +al.getStartSeq1()-1
			     <<".."<<i->getRangeS().getStart()-TSDsize
		      +al.getEndSeq1()-1<<std::endl
			     <<" LTR5:"<<i->getRangeS().getStart()
			     <<".."<<i->getRangeS().getEnd()
			     <<" LTR3:"<<i->getRangeQ().getStart()
			     <<".."<<i->getRangeQ().getEnd()<<std::endl
			     <<" TSD3:"<<s2.subseq(0,al.getEndSeq2()).toString()
			     <<" "<<i->getRangeQ().getEnd()+al.getStartSeq2()
			     <<".."<<i->getRangeQ().getEnd()+al.getEndSeq2()
			     <<std::endl;

		    std::cout<<" ==>"<<desc<<std::endl;
		    std::cout<<" TSD5:"<<s1.subseq(al.getStartSeq1()-1).toString()
			     <<" "<<i->getRangeS().getStart()-TSDsize
		      +al.getStartSeq1()-1
			     <<".."<<i->getRangeS().getStart()-TSDsize
		      +al.getEndSeq1()-1<<std::endl
			     <<" LTR5:"<<i->getRangeS().getStart()
			     <<".."<<i->getRangeS().getEnd()
			     <<" LTR3:"<<i->getRangeQ().getStart()
			     <<".."<<i->getRangeQ().getEnd()<<std::endl
			     <<" TSD3:"<<s2.subseq(0,al.getEndSeq2()).toString()
			     <<" "<<i->getRangeQ().getEnd()+al.getStartSeq2()
			     <<".."<<i->getRangeQ().getEnd()+al.getEndSeq2()
			     <<std::endl;
		    SDGBioSeq sub=s.subseq(i->getRangeS().getStart()-1,
				   i->getRangeQ().getEnd()
				   -i->getRangeS().getStart()+1);
		    sub.setDE(sub.getDE()+" "+desc+" ==> "+i->getRangeS().getStart()
			      +".."+i->getRangeQ().getEnd());
		    out2<<sub;
		  } 

	      }
	    if(!i->getRangeQ().isPlusStrand()	       
	       && i->getRangeQ().getLength()>=min_TIRsize
	       && i->getRangeQ().getStart()-i->getRangeS().getStart()
	       >=min_TIRTEsize
	       && i->getRangeQ().getStart()-i->getRangeS().getStart()
	       <=max_TIRTEsize)
	      {
		if(i->getRangeS().getStart()<TSDsize)
		  continue;
		SDGBioSeq s1=s.subseq(i->getRangeS().getStart()-TSDsize-1,TSDsize);
		SDGBioSeq s2=s.subseq(i->getRangeQ().getStart()-1,TSDsize);

		if(s1.length()!=TSDsize || s2.length() !=TSDsize)
		  continue;
		al.setSeq(s1,s2);
		al.align();
		if(al.getScore()>0 
		   && (unsigned)al.getEndSeq1()==TSDsize
		   && (unsigned)al.getStartSeq2()==1)
		  {
		    SDGString desc="TIR-TE";
		    desc+=" TE-size="+SDGString(std::fabs((int)(i->getRangeS().getStart())-i->getRangeQ().getStart())+1);
		    desc+=" TIR-size="+SDGString(i->getRangeQ().getLength());
		    desc+=" TSD-size="+SDGString(s1.subseq(al.getStartSeq1()-1).length());
		    out3<<s.getDE()<<" len="<<s.length();
		    out3<<" ==>"<<desc<<std::endl;
		    out3<<" TSD5:"<<s1.subseq(al.getStartSeq1()-1).toString()
			     <<" "<<i->getRangeS().getStart()-TSDsize
		      +al.getStartSeq1()-1
			     <<".."<<i->getRangeS().getStart()-TSDsize
		      +al.getEndSeq1()-1<<std::endl
			     <<" TIR5:"<<i->getRangeS().getStart()
			     <<".."<<i->getRangeS().getEnd()
			     <<" TIR3:"<<i->getRangeQ().getStart()
			     <<".."<<i->getRangeQ().getEnd()<<std::endl
			     <<" TSD3:"<<s2.subseq(0,al.getEndSeq2()).toString()
			     <<" "<<i->getRangeQ().getStart()+al.getStartSeq2()
			     <<".."<<i->getRangeQ().getStart()+al.getEndSeq2()
			     <<std::endl;

		    std::cout<<" ==>"<<desc<<std::endl;
		    std::cout<<" TSD5:"<<s1.subseq(al.getStartSeq1()-1).toString()
			     <<" "<<i->getRangeS().getStart()-TSDsize
		      +al.getStartSeq1()-1
			     <<".."<<i->getRangeS().getStart()-TSDsize
		      +al.getEndSeq1()-1<<std::endl
			     <<" TIR5:"<<i->getRangeS().getStart()
			     <<".."<<i->getRangeS().getEnd()
			     <<" TIR3:"<<i->getRangeQ().getStart()
			     <<".."<<i->getRangeQ().getEnd()<<std::endl
			     <<" TSD3:"<<s2.subseq(0,al.getEndSeq2()).toString()
			     <<" "<<i->getRangeQ().getStart()+al.getStartSeq2()
			     <<".."<<i->getRangeQ().getStart()+al.getEndSeq2()
			     <<std::endl;
		    SDGBioSeq sub=s.subseq(i->getRangeS().getStart()-1,
				   i->getRangeQ().getStart()
				   -i->getRangeS().getStart()+1);
		    sub.setDE(sub.getDE()+" "+desc+" ==> "+i->getRangeS().getStart()
			      +".."+i->getRangeQ().getStart());
		    out1<<sub;

		  } 

	      }
	  }
      }

    in.close();
  }
  catch(SDGException e)
    {
      std::cerr << e.message << std::endl;
    }
  return 0;
}






