/***
 *
 * BLRMatcherThreadsParameter.cpp
 *
 ***/

#include "BLRMatcherThreadsParameter.h"
#include <SDGFastaIstream.h>
#include <unistd.h> //getopt
#include <sstream>
//**************************************************************************
void BLRMatcherThreadsParameter::write(const SDGString& filename) const
{
  std::ofstream fout(filename);
  view(fout);
  fout.close();
}
//**************************************************************************
void BLRMatcherThreadsParameter::view(std::ostream& out) const
{
  BLRJoinParameter::view(out);
  out<<"============================================="<<std::endl; 
  out<<"MATCHER parameter:"<<std::endl;
  out<<"-----------------"<<std::endl;
  if(clean_before)
        out<<"Filter conflicting subjects before joining"<<std::endl;
  if(clean_after)
    out<<"Filter conflicting subjects after joining"<<std::endl;
  else
    out<<"No Filtering of conflicting subjects after joining"<<std::endl; 
  if (merge)
    out<<"Merge option enabled"<<std::endl;
  else
    out<<"No merge option enabled"<<std::endl;
  out<<" nb sets="<<nb_sets<<std::endl;
  out<<"============================================="<<std::endl; 
}
//**************************************************************************
void BLRMatcherThreadsParameter::help(void)
{  
  std::cout<<"usage: matcher.threads [option]"<<std::endl;
  BLRJoinParameter::help();
  std::cout<<"-a:\n\tall conflicting subject, default=FALSE"<<std::endl;
  std::cout<<"-x:\n\tclean conflicts after join, default=FALSE"<<std::endl;
  std::cout<<"-M:\n\tmerge (use it with clean after join option), default=FALSE"<<std::endl;
  std::cout<<"-S:\n\tsplit input data in 'n' sets, default="<<nb_sets<<std::endl;
  std::cout<<"-v:\n\tverbose, default=0"<<std::endl;
}
//**************************************************************************
void BLRMatcherThreadsParameter::parseOptArg(int numarg, char* tabarg[])
{ 
  while(1)
    {
      if(numarg==1)
	{
	  help();
	  exit(EXIT_FAILURE);
	}


	int c = getopt (numarg, tabarg, 
 			     "hm:Maxs:q:I:E:L:ji:g:d:b:B:c:v:t:S:");

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
      case 'a':
      {
          clean_before=false;
          clean_after=false;
          break;
      }
	  case 'x':
	    {
          clean_before=false;
	      clean_after=true;
	      break;
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
	  case 'v':
	    {
	      verbose=atoi(optarg);
	      break;
	    }
    case 'M':
	    {
	      merge=true;
	      break;
	    }
	case 't':
	  	{
	  		nbthread=atoi(optarg);
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
    if(!filebkS && bank_name!="<not set>")
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
