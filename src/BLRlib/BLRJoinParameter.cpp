/***
 *
 * BLRJoinParameter.cpp
 *
 ***/

#include <BLRJoinParameter.h>
#include <SDGFastaIstream.h>

#include <unistd.h> //getopt
#include <sstream>

//**************************************************************************
void BLRJoinParameter::write(const SDGString& filename) const
{
  std::ofstream fout(filename);
  view(fout);
  fout.close();
}
//**************************************************************************
void BLRJoinParameter::view(std::ostream& out) const
{
  out<<"============================================="<<std::endl; 
  out<<"JOIN parameter:"<<std::endl;
  out<<"-----------------"<<std::endl;
  out<<"Match input file: "<<match_filename<<std::endl;
  out<<"Output prefix: "<<prefix_filename<<std::endl;
  out<<"Subject bank: "<<bank_name<<std::endl;
  out<<"Query bank: "<<bankQuery_name<<std::endl;
  out<<"e-value filter: "<<eval_filter<<std::endl; 
  out<<"identity filter: "<<id_filter<<std::endl; 
  out<<"length filter: "<<len_filter<<std::endl; 
  out<<"nb thread: "<<nbthread<<std::endl;

  if(join_frag)
    {
      out<<"Join fragments"<<std::endl;
      out<<"Gap penalty to join sequence:"<<gap_pen<<std::endl;
      out<<"Distance penalty to join sequence:"<<dist_pen<<std::endl;
      out<<"Authorized overlap to join sequence:"<<overlap<<std::endl;
      out<<"Identity tolerance:"<<idtol<<std::endl;
    }
  else
    {
      out<<"No join fragments"<<std::endl;
    }
}
//**************************************************************************
void BLRJoinParameter::help(void)
{  
  std::cout<<"Options are:"<<std::endl;
  std::cout<<"-m:\n\tinput file with matches (format='align')"<<std::endl;
  std::cout<<"-q:\n\tinput file with queries (format='fasta'), default="<<bankQuery_name<<std::endl;
  std::cout<<"-s:\n\tinput file with subjects (format='fasta'), default="<<bank_name<<std::endl;
  std::cout<<"-j:\n\tjoin matches, default=false"<<std::endl;
  std::cout<<"-i:\n\tidentity tolerance to join matches, default="<<idtol<<std::endl;
  std::cout<<"-g:\n\tgap penalty to join matches, default="<<gap_pen<<std::endl;
  std::cout<<"-d:\n\tdistance penalty to join matches, default="<<dist_pen<<std::endl;
  std::cout<<"-c:\n\tauthorized overlap to join matches, default="<<overlap<<std::endl;
  std::cout<<"-E:\n\tE-value filter, default="<<eval_filter<<std::endl;
  std::cout<<"-I:\n\tidentity filter, default="<<id_filter<<std::endl;
  std::cout<<"-L:\n\tminimum length filter, default="<<len_filter<<std::endl;
  std::cout<<"-b:\n\toutput filename prefix, default="<<prefix_filename<<std::endl;
  std::cout<<"-B:\n\toutput filename prefix (no time stamp), default="<<prefix_filename<<std::endl;
  std::cout<<"-p:\n\tinput file with matches (format='path') "<<path_filename<<std::endl;
  std::cout<<"-t:\n\tinput number of thread to use, default: 1"<<std::endl;
}
//**************************************************************************
void BLRJoinParameter::parseOptArg(int numarg, char* tabarg[])
{ 
	std::cout<<" "<<std::endl;
	std::cout<<" Parse opt arg"<<std::endl;
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
// 	  {"join",no_argument, 0, 'j'},
// 	  {"id_tolerance",required_argument, 0, 'i'},
// 	  {"gap_penality",required_argument, 0, 'g'},
// 	  {"distance_penality",required_argument, 0, 'd'},
// 	  {"authorized_overlap",required_argument, 0, 'c'},
// 	  {"prefix_file",required_argument, 0, 'b'},
// 	  {"prefix_file_no_stamp",required_argument, 0, 'B'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;
	
// 	int c = getopt_long (numarg, tabarg, 
// 			     "hm:s:q:I:E:L:ji:g:d:b:B:c:",
// 			 long_options, &option_index);

      int c = getopt (numarg, tabarg, 
		      "hm:s:q:I:E:L:ji:g:d:b:B:c:p:t:");
	
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

    if(bank_name=="<not set>")
      bank_name=bankQuery_name;
    if(bankQuery_name=="<not set>")
      bankQuery_name=bank_name;
    
    if(prefix_filename=="<not set>")
      {
	  prefix_filename=match_filename;
      }

    std::ifstream filebkQ(bankQuery_name);
    if(!filebkQ && bankQuery_name!="<not set>")
      {
	char *p_dirs=getenv("BLASTER_BANK_DIRS");
	if(p_dirs!=NULL)
	  {
	    std::cout<<"run_test_search_wSW bank:"<<bankQuery_name<<"..."<<std::flush;
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
    if(!filebkS && bank_name!="<not set>")
      {
	char *p_dirs=getenv("BLASTER_BANK_DIRS");
	if(p_dirs!=NULL)
	  {
	    std::cout<<"run_test_search_wSW bank:"<<bank_name<<"..."<<std::flush;
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
    
    if(!filebkS && bank_name!="<not set>")
      {
	std::cout<<"Error bank:"<<bank_name<<" not found!!"<<std::endl;
	exit(EXIT_FAILURE);
      }
    if(!filebkQ && bankQuery_name!="<not set>")
      {
	std::cout<<"Error bank:"<<bankQuery_name<<" not found!!"<<std::endl;
	exit(EXIT_FAILURE);
      }

    filebkQ.close();
    filebkS.close();
     
}















