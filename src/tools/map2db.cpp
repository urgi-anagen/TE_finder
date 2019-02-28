#include <iostream>
#include <unistd.h> //getopt

#include <SDGError.h>
#include <SDGString.h>

#include "RangeMap.h"

unsigned flank_size=0;
bool merge=false;
unsigned verbose=0;

void help( void )
{
  std::cerr<<std::endl<<"usage: map2db"
  		<<" [<options>] <mapfile> <fasta sequences>"<<std::endl
  		<<" options:"<<std::endl
      <<"   -h, --help:\n\t this help"<<std::endl
      <<"   -m, --merge:\n\t merge overlapping range, default: false"<<std::endl
      <<"   -s, --size:\n\t flank size to extend (>0), default: "<<flank_size<<std::endl
      <<"   -v, --verbose:\n\t verbosity level (0/1), default: "<<verbose<<std::endl
      <<std::endl;
};

int main( int argc, char* argv[] )
{
  try
  {
  	int c;
  	while( true )
  	{
// 	static struct option long_options[] =
// 	{
// 	  {"help",no_argument, 0, 'h'},
// 	  {"merge",no_argument, 0, 'm'},
// 	  {"size",required_argument, 0, 's'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;

// 	c = getopt_long (argc, argv, "hms:",
// 			 long_options, &option_index);

  		c = getopt(argc, argv, "hms:v:");

  		/* Detect the end of the options. */
  		if(c == -1)
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
  			merge=true;
  			break;
  		}
  		case 's':
  		{
  			flank_size=atoi(optarg);
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

  	SDGString map_filename = argv[optind];
    SDGString fasta_filename = argv[++optind];

    /* Print any remaining command line arguments (not options). */
    if (++optind < argc)
    {
    	help();
    	std::cout<<"non-option ARGV-elements: "<<std::endl;
    	while (optind < argc)
    		std::cout<<argv[optind++]<<std::endl;
    	return 1;
    }

    RangeMap range_container;

    if( verbose > 0 )
    	std::cout<<"load map... "<<std::flush;
    range_container.load(map_filename);
    if( verbose > 0 )
    	std::cout<<"done!"<<std::endl<<std::flush;

    if( verbose > 0 )
    {
    	range_container.view();
    	unsigned countRange = range_container.getCountRange();
    	std::cout<<countRange<<" found"<<std::endl;
    }

    if( flank_size > 0 )
      range_container.extend( flank_size );

    if( merge )
      range_container.merge();

    SDGString dbname( map_filename );
    dbname = dbname.afterlast("/")+".flank_size"+SDGString(flank_size)+".fa";
    if( verbose > 0 )
    	std::cout<<"write sequences"<<std::endl<<std::flush;
    range_container.writeSeq( dbname, fasta_filename );
  }

  catch(SDGException e)
  {
  	std::cerr<<e.message<<std::endl;
  	exit(EXIT_FAILURE);
  }
  exit(EXIT_SUCCESS);
};
