/***
 *
 * BLRGrouperParameter.cpp
 *
 ***/

#include <BLRGrouperParameter.h>
#include <SDGFastaIstream.h>
#include <unistd.h> //getopt
#include <sstream>

//**************************************************************************
void BLRGrouperParameter::view(std::ostream& out) const
{
  BLRJoinParameter::view(out);
  out<<"============================================="<<std::endl;
  out<<"GROUPER parameter:"<<std::endl;
  out<<"-----------------"<<std::endl;
  out<<" graphfilter="<<graphfilter<<std::endl;
  out<<" coverage="<<coverage<<std::endl;
  out<<" sizefilter="<<sizefilter<<std::endl;
  out<<" includefilter="<<includefilter<<std::endl;
  out<<" nb sets="<<nb_sets<<std::endl;
  out<<"============================================="<<std::endl;
}
//**************************************************************************
void BLRGrouperParameter::help(void)
{
  std::cout<<"usage: grouper [ options ]"<<std::endl;
  BLRJoinParameter::help();
  std::cout<<"-G:\n\tmin coverage length for connecting groups in a cluster, default="<<graphfilter<<std::endl;
  std::cout<<"-C:\n\tcoverage for group construction, default="<<coverage<<std::endl;
  std::cout<<"\tinterpreted in % if <= 1; in bp otherwise"<<std::endl;
  std::cout<<"-Z:\n\tgroup size filter, default="<<sizefilter<<std::endl;
  std::cout<<"-X:\n\tkeep groups where at least X members are not included in others, default="<<includefilter<<std::endl;
  std::cout<<"-S:\n\tsplit input data in 'n' sets, default="<<nb_sets<<std::endl;
  std::cout<<"-v:\n\tverbose, default="<<verbose<<std::endl;
}
//**************************************************************************
void BLRGrouperParameter::parseOptArg(int numarg, char* tabarg[])
{
  while(1)
    {
      if(numarg==1)
	{
	  help();
	  exit(EXIT_FAILURE);
	}

// 	static struct option long_options[] =
// 	{
// 	  {"help",no_argument, 0, 'h'},
// 	  {"match",required_argument, 0, 'm'},
// 	  {"subject_bank",required_argument, 0, 's'},
// 	  {"query_bank",required_argument, 0, 'q'},
// 	  {"e-value_filter",required_argument, 0, 'E'},
// 	  {"identity_filter",required_argument, 0, 'I'},
// 	  {"length_filter",required_argument, 0, 'L'},
// 	  {"sizefilter",required_argument, 0, 'Z'},
// 	  {"join",no_argument, 0, 'j'},
// 	  {"id_tolerance",required_argument, 0, 'i'},
// 	  {"gap_penality",required_argument, 0, 'g'},
// 	  {"distance_penality",required_argument, 0, 'd'},
// 	  {"authorized_overlap",required_argument, 0, 'c'},
// 	  {"prefix_file",required_argument, 0, 'b'},
// 	  {"prefix_file_no_stamp",required_argument, 0, 'B'},
// 	  {"include",no_argument, 0, 'X'},
// 	  {"graphfilter",required_argument, 0, 'G'},
// 	  {"coverage",required_argument, 0, 'C'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;

// 	int c = getopt_long (numarg, tabarg,
// 			     "hm:s:q:I:E:L:Z:ji:g:d:b:B:c:XG:C:",
// 			 long_options, &option_index);
	int c = getopt (numarg, tabarg,
			     "hm:s:q:I:E:L:M:Z:ji:g:d:b:B:c:X:G:C:p:t:v:S:");

	/* Detect the end of the options. */
	if (c == -1)
	  break;

	switch (c)
	  {
	  case 'h':
	    {
	      help();
	      exit(EXIT_FAILURE);
	    }
	  case 'm':
	    {
	      match_filename=optarg;
	      break;
	    }
	  case 'b':
	    {
	      prefix_filename=optarg;
	      time_stamp=true;
	      break;
	    }
	  case 'B':
	    {
	      prefix_filename=optarg;
	      time_stamp=false;
	      break;
	    }
	  case 'q':
	    {
	      bankQuery_name=optarg;
	      break;
	    }
	  case 's':
	    {
	      bank_name=optarg;
	      break;
	    }
	  case 'j':
	    {
	      join_frag=true;
	      break;
	    }
	  case 'i':
	    {
	      idtol=atof(optarg);
	      break;
	    }
	  case 'g':
	    {
	      gap_pen=atof(optarg);
	      break;
	    }
	  case 'd':
	    {
	      dist_pen=atof(optarg);
	      break;
	    }
	  case 'c':
	    {
	      overlap=atoi(optarg);
	      break;
	    }
	  case 'E':
	    {
	      eval_filter=atof(optarg);
	      break;
	    }
	  case 'I':
	    {
	      id_filter=atof(optarg);
	      break;
	    }
	  case 'L':
	    {
	      len_filter=atoi(optarg);
	      break;
	    }
	  case 'Z':
	    {
	      sizefilter=atoi(optarg);
	      break;
	    }
	  case 'C':
	    {
	      coverage=atof(optarg);
	      break;
	    }
	  case 'G':
	    {
	      graphfilter=atoi(optarg);
	      break;
	    }
	  case 'X':
	    {
	      includefilter=atoi(optarg);
	      break;
	    }
 	  case 'p':
	    {
	      path_filename=optarg;
	      loadPath=true;
	      break;
	    }
	  case 't':
	  	{
	  		nbthread=atoi(optarg);
	  	    break;
	  	}
	  case 'v':
	    {
		  verbose=atoi(optarg);
		  break;
	    }
	  case 'S':
	    {
		  nb_sets=atoi(optarg);
		  break;
	    }
	  case '?':
	    help();
	    exit(EXIT_FAILURE);
	  default:
	    abort ();
	  }
      }

    /* Print any remaining command line arguments (not options). */
    if (optind < numarg)
      {
	help();
	std::cout<<"non-option ARGV-elements: "<<std::endl;
	while (optind < numarg)
	  std::cout<<tabarg[optind++]<<std::endl;
	exit(EXIT_FAILURE);
      }

//    if( includefilter >= sizefilter )
//    {
//    	std::cout<<"ERROR: include filter (-X "<<includefilter<<") should be strictly lower than size filter"<<"(-Z "<<sizefilter<<")"<<std::endl;
//    	exit( EXIT_FAILURE );
//    }

    if(bank_name=="<not set>")
      bank_name=bankQuery_name;
    if(bankQuery_name=="<not set>")
      bankQuery_name=bank_name;
   
    if(prefix_filename=="<not set>")
      {
	if (path_filename =="<not set>")
		prefix_filename=match_filename.afterlast("/");
	else
		prefix_filename=path_filename.afterlast("/");
      }


       std::ifstream filebkQ(bankQuery_name);
    if(!filebkQ)
      {
	char *p_dirs=getenv("BLASTER_BANK_DIRS");
	if(p_dirs!=NULL)
	  {
	    std::cout<<"search bank:"<<bankQuery_name<<"..."<<std::flush;
	    std::string dirs=p_dirs;
	    std::istringstream stream(dirs);
	    SDGString path, pathfile;
	    while(!filebkQ && stream>>path)
	      {
		std::cout<<path<<std::endl;
		filebkQ.clear();
		pathfile=path+"/"+bankQuery_name;
		filebkQ.open(pathfile);
	      }
	    if(filebkQ)
	      {
		bankQuery_name=path+"/"+bankQuery_name;
		std::cout<<"find:"<<bankQuery_name<<std::endl;
	      }
	  }
      }

    std::ifstream filebkS(bank_name);
    if(!filebkS)
      {
	char *p_dirs=getenv("BLASTER_BANK_DIRS");
	if(p_dirs!=NULL)
	  {
	    std::cout<<"search bank:"<<bank_name<<"..."<<std::flush;
	    std::string dirs=p_dirs;
	    std::istringstream stream(dirs);
	    SDGString path, pathfile;
	    while(!filebkS && stream>>path)
	      {
		filebkS.clear();
		pathfile=path+"/"+bank_name;
		filebkS.open(pathfile);
	      }
	    if(filebkS)
	      {
		bank_name=path+"/"+bank_name;
		std::cout<<"find:"<<bank_name<<std::endl;
	      }
	  }
      }

    if(!filebkS)
      std::cout<<"Error bank:"<<bank_name<<" not found!!"<<std::endl;

    if(!filebkQ)
      std::cout<<"Error bank:"<<bankQuery_name<<" not found!!"<<std::endl;

    if(!filebkQ || !filebkS) exit(EXIT_FAILURE);

    filebkQ.close();
    filebkS.close();

}
