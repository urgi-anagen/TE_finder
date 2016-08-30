#include <unistd.h> //getopt
#include <sstream>
#include <fstream>

#include <SDGFastaIstream.h>
#include <SDGFastaOstream.h>
#include <SDGMemBioSeq.h>

#include "Lalign.h"

int match=10,mismh=16,gap_extend=32,gap_open=32;
double id_thres=0.8;
unsigned max_align_len=5000, min_len=100;
    
void help(void)
{
  std::cerr<<"usage: ltrsearch"
      <<" [<options>] <fasta sequence database>"<<std::endl
      <<" options:"<<std::endl
      <<"   -h, --help:\n\t this help"<<std::endl
      <<"   -m, --match:\n\t match bonus (>0), default: "
      <<match<<std::endl
      <<"   -d, --mismatch:\n\t mismatch penalty (>0), default: "
      <<mismh<<std::endl
      <<"   -g, --gapopen:\n\t gap open penalty (>0), default: "
      <<gap_open<<std::endl
      <<"   -e, --gapextend:\n\t gap extend penalty (>0), default: "
      <<gap_extend<<std::endl
      <<"   -i, --identity:\n\t min identity threshold, default: "
      <<id_thres<<std::endl
      <<"   -l, --length:\n\t min itr length, default: "<<min_len<<std::endl;
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
// 	  {"match",required_argument, 0, 'm'},
// 	  {"mismatch",required_argument, 0, 'd'},
// 	  {"gapopen",required_argument, 0, 'g'},
// 	  {"gapextend",required_argument, 0, 'e'},
// 	  {"identity",required_argument, 0, 'i'},
// 	  {"length",required_argument, 0, 'l'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;
	
// 	c = getopt_long (argc, argv, "hm:d:g:e:l:i:",
// 			 long_options, &option_index);
	
	c = getopt(argc, argv, "hm:d:g:e:l:i:");
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
	  case 'l':
	    {
	      min_len=atoi(optarg);
	      break;
	    }
	  case 'i':
	    {
	      id_thres=atof(optarg);
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
     

    Lalign al;
    al.setMismatch(match,mismh);
    al.setGap(gap_open, gap_extend);
    
    SDGFastaIstream in(filename);
    SDGString filenameout=filename.afterlast("/");
    filenameout=filenameout+".ltr";
    SDGFastaOstream out(filenameout);
    
    int count_ltr=0;
    int count_all=0;
    while(in)
      {
	SDGBioSeq seq;
	std::cout<<"load sequence #"<<++count_all;
	in>>seq;
	unsigned len=seq.length();
	std::cout<<" "<<seq.getDE()<<"\tlength="<<len<<std::endl;

	SDGBioSeq s1,s2;

	unsigned align_len=len/2<max_align_len?len/2:max_align_len;
	s1=seq.subseq(0,align_len);
	s2=seq.subseq(len-align_len,len);
	al.setSeq(s1,s2);	  
	al.align();
	unsigned length=al.getEndSeq1()-al.getStartSeq1()+1;

	unsigned r1s=al.getStartSeq1();
	unsigned r1e=al.getEndSeq1();
	unsigned r2s=len-align_len+al.getStartSeq2();
	unsigned r2e=len-align_len+al.getEndSeq2();

	std::cout<<std::endl;
	al.view();

	if(al.getScore()==0) continue;

	std::cout<<"\n\tQuality="<<(double)al.getScore()/length
	    <<"\n\t("<<r1s
	    <<","<<r1e<<")"
	    <<"..("<<r2s
	    <<","<<r2e<<")"
	    <<" Length ltr="<<length
	    <<"\n\tLength TE="<<r2e-r1s+1
	    <<std::endl;

	if( length>=min_len && al.getIdentity()>id_thres 
	    && r1s<5 && len-r2e<5)
	  {
	    std::ostringstream header;
	      header<<seq.getDE()
	      <<" LTR("<<r1s
	      <<","<<r1e<<")"
	      <<"..("<<r2s
	      <<","<<r2e<<")"
	      <<" Length ltr="<<length;

	    seq.setDE(header.str());
	    out<<seq;
	    std::cout<<"\t\t==> saved !";
	    count_ltr++;
	  }
	std::cout<<std::endl;
      }
    std::cout<<"\nSequences treated:"<<count_all<<std::endl;
    std::cout<<"Sequences found:"<<count_ltr
	<<" ("<<(double)count_ltr/count_all<<")"<<std::endl;
  }
  catch(SDGException e)
    {
      std::cerr << e.message << std::endl;
    }
  
  return 0;
}












