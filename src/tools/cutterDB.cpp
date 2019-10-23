#include <unistd.h> //getopt
#include <iostream>
#include <cstdlib>
#include <Cutter.h>

int length=50000,over=100,word=11,verbose=0;

void help(void)
{
  std::cerr<<"usage: cutterDB"
	   <<" [<options>] <fasta database sequence>"<<std::endl
	   <<" options:"<<std::endl
	   <<"   -h, --help:\n\t this help"<<std::endl
	   <<"   -l, --length:\n\t length, default: "<<length<<std::endl
	   <<"   -o, --overlap:\n\t overlap, default: "<<over<<std::endl
	   <<"   -w, --word:\n\t N stretch word length (0 no detection), default: "<<word<<std::endl
	   <<"   -v, --verbose:\n\t verbosity level, default: "<<verbose<<std::endl;
};

//-------------------------------------------------------
int main(int argc, char* argv[])
{ 
  try
    {
      while (1)
	{

 	  int c = getopt (argc, argv, "hl:o:w:v:");
	  
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
	    case 'l':
	      {
		length=atoi(optarg);
		break;
	      }
	    case 'o':
	      {
		over=atoi(optarg);
		break;
	      }
	    case 'w':
	      {
		word=atoi(optarg);
		break;
	      }
	    case 'v':
	      {
		verbose=atoi(optarg);
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

      Cutter cutter;            //!< parameter of bank cut-out

      cutter.setLength(length);
      cutter.setOver(over);
      cutter.setWord(word);
      if(!cutter.check(filename,verbose))
	{
	  SDGString bqName=cutter.cutDB(filename,verbose);
	  if(verbose>0)
	    std::cout<<"\n\nbank '"<<bqName<<"' created!!\n";
	}
    }
  catch(SDGException e)
    {
      std::cout<<e.message<<std::endl;
    }
}
