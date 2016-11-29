#include <unistd.h> //getopt
#include <fstream>
#include <SDGFastaIstream.h>
#include <SDGMemBioSeq.h>

#include "HashRepeat.h"
unsigned word_len=6, extend=10;
bool reputer=false;
SDGString psfilename="hrepeat_out.ps";
    
void help(void)
{
  std::cerr<<"usage: hrepeat"
      <<" [<options>] <fasta sequence database>"<<std::endl
      <<" options:"<<std::endl
      <<"   -h, --help:\n\t this help"<<std::endl
      <<"   -r, --reputer:\n\t reputer search"<<std::endl
      <<"   -x, --extend:\n\t max size for boundaries extension, default: "
      <<extend<<std::endl
	   <<"   -w, --word:\n\t word length, default: "<<word_len<<std::endl;
//       <<"   -o, --psout:\n\t filename for postscript output,"
//       <<" default:"<<psfilename<<std::endl;
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
// 	  {"extend",required_argument, 0, 'x'},
// 	  {"word",required_argument, 0, 'w'},
// 	  {"psout",required_argument, 0, 'o'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;
	
// 	c = getopt_long (argc, argv, "hrx:w:o:",
// 			 long_options, &option_index);
	c = getopt(argc, argv, "hrx:w:");
	
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
	  case 'r':
	    {
	      reputer=true;
	      break;
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
	  case 'o':
	    {
	      psfilename=optarg;
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
     

    SDGFastaIstream in(filename);
    SDGBioSeq s;
    HashRepeat hrep(word_len,extend);
 //    if(psfilename!="")
//       hrep.initPlot(psfilename);

    SDGString outfilename=filename+".hrepeat.set";
    std::ofstream out(outfilename);
    unsigned count=0;
    while(in) 
      {
	in>>s;
	std::cout<<s.getDE()<<" len="<<s.length()<<std::endl;
//	if(reputer)
//	  hrep.search_reputer(s,1);
//	else
	  hrep.search(s);

//	float compress=hrep.compress();
//	std::cout<<"**"<<s.getDE()<<" compress ratio="<<compress<<std::endl;
	hrep.write();
	count=hrep.writeSet(out,count);
// 	if(psfilename!="")
// 	  hrep.plot(s.getDE()+" compress="+SDGString(compress));
	std::cout<<"ok!"<<std::endl;
      }
    out.close();
    in.close();
  }
  catch(SDGException e)
    {
      std::cerr << e.message << std::endl;
    }
  return 0;
}

