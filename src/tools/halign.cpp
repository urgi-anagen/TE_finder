#include <unistd.h> //getopt
#include <fstream>
#include <SDGFastaIstream.h>
#include <SDGMemBioSeq.h>

#include "HashAlign.h"
double match=1,mismh=0.8,gap_extend=0.4,gap_open=1.6;
unsigned word_len=13, nb_frag=100, min_size=100, overlap=0;
SDGString psfilename="halign_out.ps",infilename;
    
void help(void)
{
  std::cerr<<"usage: halign"
      <<" [<options>] <fasta sequence> <fasta db>"<<std::endl
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
      <<"   -c, --overlap:\n\t authorized overlap, default: "
      <<overlap<<std::endl
      <<"   -w, --word:\n\t word length, default: "<<word_len<<std::endl
      <<"   -n, --nfrag:\n\t max number of fragment to align, default: "
	<<nb_frag<<std::endl
      <<"   -s, --minsize:\n\t min alignment size to report, default: "
	<<min_size<<std::endl
//       <<"   -o, --psout:\n\t filename for postscript output,"
//       <<" default:"<<psfilename<<std::endl
      <<"   -f, --fmatches:\n\t filename for input matches,"
      <<" default: no input"<<std::endl;
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
// 	  {"overlap",required_argument, 0, 'c'},
// 	  {"word",required_argument, 0, 'w'},
// 	  {"nfrag",required_argument, 0, 'n'},
// 	  {"minsize",required_argument, 0, 's'},
// 	  {"psout",required_argument, 0, 'o'},
// 	  {"fmatch",required_argument, 0, 'f'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;
	
// 	c = getopt_long (argc, argv, "hm:d:g:e:c:w:n:s:o:f:",
// 			 long_options, &option_index);
	c = getopt (argc, argv, "hm:d:g:e:c:w:n:s:f:");
	
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
	      match=atof(optarg);
	      break;
	    }
	  case 'd':
	    {
	      mismh=atof(optarg);
	      break;
	    }
	  case 'g':
	    {
	      gap_open=atof(optarg);
	      break;
	    }
	  case 'e':
	    {
	      gap_extend=atof(optarg);
	      break;
	    }
	  case 'c':
	    {
	      overlap=atoi(optarg);
	      break;
	    }
	  case 'w':
	    {
	      word_len=atoi(optarg);
	      break;
	    }
	  case 'n':
	    {
	      nb_frag=atoi(optarg);
	      break;
	    }
	  case 's':
	    {
	      min_size=atoi(optarg);
	      break;
	    }
	  case 'o':
	    {
	      psfilename=optarg;
	      break;
	    }
	  case 'f':
	    {
	      infilename=optarg;
	      break;
	    }
	  case '?':
	    help();
	    return 1;
	  default:
	    abort ();
	  }
      }

    std::cout<<infilename<<std::endl;
    SDGString filename1,filename2;
    if(infilename=="")
      {
	if(argc==optind)
	  {
	    help();
	    return 1;
	  }

	filename1=argv[optind];
	filename2=argv[++optind];
      }
    /* Print any remaining command line arguments (not options). */
    if (++optind < argc)
      {
	help();
	std::cout<<"non-option ARGV-elements: "<<std::endl;
	while (optind < argc)
	  std::cout<<argv[optind++]<<std::endl;
	return 1;
      }
     
    HashAlign hal1,hal2;
//     if(psfilename!="")
//       {
// 	hal1.initPlot(psfilename);
// 	hal2.initPlot(psfilename);
//       }

    if(infilename!="")
      {
	hal1.load_matches(infilename);
	hal1.fragAlign(match,mismh,gap_open,
		      gap_extend,overlap,nb_frag);
	hal1.print(min_size,std::cout);
	//    hal1.write(std::cout);
// 	if(psfilename!="")
// 	  hal1.plot();
	std::cout<<"ok!\n"<<std::endl;
      }
    else
      {
	SDGFastaIstream in1(filename1);
	SDGBioSeq s1;
	if(in1)
	  in1>>s1;
	std::cout<<s1.getDE()<<" len:"<<s1.length()<<" read!"<<std::endl;
	
	SDGFastaIstream in2(filename2);
	while(in2)
	  {
	    SDGBioSeq s2;
	    if(in2)
	      in2>>s2;
	    std::cout<<s2.getDE()<<" len:"<<s2.length()<<" read!"<<std::endl;

	    hal1.setSeq(s1,s2);
	    hal1.search(word_len);
	    hal1.fragAlign(match,mismh,gap_open,
			  gap_extend,overlap,nb_frag);
	    unsigned sc1=hal1.score();

	    s2=s2.complement();
	    hal2.setSeq(s1,s2);
	    hal2.search(word_len);
	    hal2.fragAlign(match,mismh,gap_open,
			  gap_extend,overlap,nb_frag);
	    unsigned sc2=hal2.score();

	    if(sc1>=sc2)
	      {
		hal1.print(min_size,std::cout);
		//    hal1.write(std::cout);
// 		if(psfilename!="")
// 		  hal1.plot();
	      }
	    else
	      {
		hal2.print(min_size,std::cout);
		//    hal2.write(std::cout);
// 		if(psfilename!="")
// 		  hal2.plot();
	      }
	    std::cout<<"ok!\n"<<std::endl;
	  }
      }
  }
  catch(SDGException e)
    {
      std::cerr << e.message << std::endl;
    }
  return 0;
}












