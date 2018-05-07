#include <getopt.h>
#include <fstream>
#include <SDGString.h>
#include <SDGFastaIstream.h>
#include <SDGMemBioSeq.h>

#include "Hasher.h"

double match=1,mismh=0.8,gap_extend=0.4,gap_open=1.6,
  filter_cutoff=0.0;
unsigned word_len=10, nb_frag=100000, min_size=20,overlap=0;
int ext_len=-1;
SDGString outfilename="";
bool repeat_mode=false;
    
void help(void)
{
  std::cerr<<"usage: hasher"
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
      <<"   -c, --overlap:\n\t authorized overlap percentage (<1.0), default: "
      <<overlap<<std::endl
      <<"   -w, --word:\n\t word length, default: "<<word_len<<std::endl
      <<"   -n, --nfrag:\n\t max number of fragment to align, default: "
	<<nb_frag<<std::endl
      <<"   -s, --minsize:\n\t min alignment size to report, default: "
	<<min_size<<std::endl
      <<"   -r, --repeat:\n\t repeat mode, default: FALSE"
      <<std::endl
      <<"   -o, --file_out:\n\t filename for output,"
      <<" default: None"<<std::endl
      <<"   -f, --fmatches:\n\t filename for input matches,"
      <<" default: no input"<<std::endl;
};

int main(int argc, char* argv[])
{
  try{
    
    if(argc==1)
      {
	help();
	exit(EXIT_SUCCESS);
      } 
    int c;
    while (1)
      {
	static struct option long_options[] =
	{
	  {"help",no_argument, 0, 'h'},
	  {"match",required_argument, 0, 'm'},
	  {"mismatch",required_argument, 0, 'd'},
	  {"gapopen",required_argument, 0, 'g'},
	  {"gapextend",required_argument, 0, 'e'},
	  {"overlap",required_argument, 0, 'c'},
	  {"word",required_argument, 0, 'w'},
	  {"nfrag",required_argument, 0, 'n'},
	  {"minsize",required_argument, 0, 's'},
	  {"repeat",no_argument, 0, 'r'},
	  {"file_out",required_argument, 0, 'o'},
	  {"filter",required_argument, 0, 'f'},
	  {0, 0, 0, 0}
	};
	/* `getopt_long' stores the option index here. */
	int option_index = 0;
	
	c = getopt_long (argc, argv, "hm:d:g:e:c:w:n:s:ro:f:",
			 long_options, &option_index);
	
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
	  case 'r':
	    {
	      repeat_mode=true;
	      break;
	    }
	  case 'o':
	    {
	      outfilename=optarg;
	      break;
	    }
	  case 'f':
	    {
	      filter_cutoff=atof(optarg);
	      break;
	    }
	  case '?':
	    help();
	    return 1;
	  default:
	    abort ();
	  }
      }

    SDGString filename1,filename2;
    std::ofstream out;

    filename1=argv[optind];
    if(!repeat_mode) filename2=argv[++optind];

    /* Print any remaining command line arguments (not options). */
    if (++optind < argc)
      {
	help();
	std::cout<<"non-option ARGV-elements: "<<std::endl;
	while (optind < argc)
	  std::cout<<argv[optind++]<<std::endl;
	return 1;
      }
     
    Hasher hsrch(word_len,2,nb_frag, filter_cutoff, repeat_mode, ext_len);
    if(repeat_mode) hsrch.load(filename1);
    else hsrch.load(filename2);

    if(outfilename!="")
    	out.open(outfilename);

    std::list<RangePair> rp_list;
	
    SDGFastaIstream in(filename1);
    while(in)
      {
		SDGBioSeq s;
		if(in)
		  in>>s;
		std::cout<<s.getDE()<<" len:"<<s.length()<<" read!"<<std::endl;

		hsrch.search(s);
//	 	std::cout<<"Found fragments:"<<std::endl;
//	 	hsrch.print_frag(std::cout);
		hsrch.fragAlign(match,mismh,gap_open,gap_extend,overlap,false);
	// 	std::cout<<"Aligned fragments:"<<std::endl;
	// 	hsrch.print(min_size,std::cout);
	//	hsrch.extend(min_size);
	//	std::cout<<"Extended fragments:"<<std::endl;
		hsrch.print(min_size,std::cout);
		if(outfilename!="")
			hsrch.write_align(min_size,out);
		std::cout<<"ok!\n"<<std::endl;
      }
	out.close();
  }
  catch(SDGException e)
    {
      std::cerr << e.message << std::endl;
    }
//   catch(...)
//     {
//       std::cerr <<"catched an unknown exception !"<< std::endl;
//     }
  return 0;
}












