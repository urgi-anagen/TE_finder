/***
 *
 * BLRBlasterParameter.cpp
 *
 ***/

#include <sys/stat.h>
#include "BLRBlasterParameter.h"
#include <SDGFastaIstream.h>
#include <unistd.h> //getopt
#include <sstream>

//**************************************************************************
/*!- Function start or restart blaster
 *- Call by blaster.cpp with number of argument(numarg) and array of this
 *- find in command line and return false if error occur and true if not. */
bool BLRBlasterParameter::start (int numarg, char* tabarg[])
{
  BLRBlasterParameter last;   //! Set argument of last run

  // Parse the command line
  parseOptArg(numarg, tabarg);
  if(verbose>0)
	  view(std::cout);

  if(bank_name=="<not set>")
    {
      std::cerr<<"ERROR: Bank not set"<<std::endl;
      throw SDGException(NULL,"parameter parse error",-1);
    }

  std::ifstream filebkQ(bankQuery_name);
  if(!filebkQ)
    {
      char *p_dirs=getenv("BLASTER_BANK_DIRS");
      if(p_dirs!=NULL)
	{
      if(verbose>0)
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
	      if(verbose>0)
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
      if(verbose>0)
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
	      if(verbose>0)
	    	  std::cout<<"find:"<<bank_name<<std::endl;
	    }
	}
    }

  if(!filebkS)
      std::cout<<"Error bank:"<<bank_name<<" not found!!"<<std::endl;

  if(!filebkQ)
      std::cout<<"Error bank:"<<bankQuery_name<<" not found!!"<<std::endl;

  if(!filebkQ || !filebkS)
      throw SDGException(NULL,"parameter parse error",-1);

  filebkQ.close();
  filebkS.close();

  if(bankQuery_cut=="<not set>")
    bankQuery_cut=bankQuery_name+getCutExtention();
  if(bank_cut=="<not set>")
    bank_cut=bank_name+getCutExtention();

  if(blaster_filename=="<not set>")
    {
      if(bankQuery_name==bank_name)
	blaster_filename=bankQuery_name.afterlast("/");
      else
	blaster_filename=bankQuery_name.afterlast("/")
	  +"_vs_"+bank_name.afterlast("/");

      blaster_filename=blaster_filename+"-"+type_blast;
    }

  loadPrevious(last);

  if( ! cutter.check( bank_name, verbose ) )
    cutter.cutDB( bank_name, verbose );
  if( bankQuery_name != bank_name )
    if( ! cutter.check( bankQuery_name, verbose ) )
      cutter.cutDB( bankQuery_name, verbose );

  numdo=blastcount();

  if(restart)
  {
	  std::ostringstream cmd;
	  struct stat stFileInfo;
	  int intStat;
	  intStat = stat( blaster_filename+".raw", &stFileInfo );
	  if( intStat == 0 )
	  {
		  cmd<<"rm "<<blaster_filename<<".raw";
		  system( cmd.str().c_str() );
		  cmd.str("");
	  }
	  intStat = stat( blaster_filename+".seq_treated", &stFileInfo );
	  if( intStat == 0 )
	  {
		  cmd<<"rm "<<blaster_filename<<".seq_treated";
		  system( cmd.str().c_str() );
	  }
  }

  save();

  return true;
}
//**************************************************************************
void BLRBlasterParameter::loadPrevious(BLRBlasterParameter& last)
{
  try
    {
      last.load(parameter_filename);
      last.parameter_filename=parameter_filename;
      last.blaster_filename=blaster_filename;
      parameter_file=true;
      first_run=false;
    }
  catch(SDGException& e)
    {
      last.reset();
      first_run=true;
    }
}
//**************************************************************************
void BLRBlasterParameter::write(const SDGString& filename) const
{
  std::ofstream fout(filename);
  view(fout);
  fout.close();
}
//**************************************************************************
void BLRBlasterParameter::view(std::ostream& out) const
{
  out<<"============================================="<<std::endl;
  out<<"BLASTER parameter:"<<std::endl;
  out<<"-----------------"<<std::endl;
  out<<"Blaster output prefix: "<<blaster_filename<<std::endl;
  out<<"Subject bank: "<<bank_name<<std::endl;
  out<<"Query bank: "<<bankQuery_name<<std::endl;
  if(all_by_all)
    out<<"all_by_all: TRUE"<<std::endl;
  else
    out<<"all_by_all: FALSE"<<std::endl;
  if(is_wuBlast)
    out<<"WU-BLAST"<<std::endl;
  if(is_ncbiBlast)
    out<<"NCBI-BLAST"<<std::endl;
  if(is_ncbiBlastPlus)
    out<<"NCBI-BLAST-PLUS"<<std::endl;
  out<<"Blast name: "<<type_blast<<std::endl;
  out<<"Blast parameter: "<<option_blast<<std::endl;
  out<<"Nb concurrent processes: "<<nb_process<<std::endl;
  out<<"Cutter bank:\tlength="<<cutter.getLength()<<std::endl;
  out<<"\t \t over="<<cutter.getOver()<<std::endl;
  out<<"\t \t word="<<cutter.getWord()<<std::endl;
  out<<"\t \t extention="<<cutter.getExtention()<<std::endl;
  out<<"Bank cut: "<<bank_cut<<std::endl;
  out<<"Query cut : "<<bankQuery_cut<<std::endl;
  out<<"sensitivity level: "<<sensitivity_lvl<<std::endl;
  out<<"e-value filter: "<<eval_filter<<std::endl;
  out<<"identity filter: "<<id_filter<<std::endl;
  out<<"length filter: "<<len_filter<<std::endl;
  if(prepare)
    out<<"**** prepare only ! *****"<<std::endl;
  if(cleanTmpFiles)
	  out<<"clean temporary files"<<std::endl;
  out<<"============================================="<<std::endl;
}
//**************************************************************************
//---------------------------------------------------
//! file checks and count the number of blast to do -
//---------------------------------------------------
int BLRBlasterParameter::blastcount()
{
  std::ifstream bankin(bank_cut);       //! bank file in fasta format
  // inspect status of file
  if(bankin.bad())
    {
      std::cerr<<"error wrong file name bank_cut:"<<bank_cut<<"!!"<<std::endl;
      return 0;
    }
  bankin.close();

  std::ifstream Queryin(bankQuery_cut); //! query file in fasta format
  if(Queryin.bad())
    {
      std::cerr<<"error wrong file name bankQuery_cut:"<<bankQuery_cut<<"!!"<<std::endl;
      return 0;
    }
  unsigned count_seq=0;
  char buff[256];
  while (Queryin)
    {
	  Queryin.getline(buff,256);
	  if (*(buff) == '>')
		  count_seq=atol(strtok(&buff[1]," "));
    }
  Queryin.close();
  return count_seq;
}
//**************************************************************************
void BLRBlasterParameter::help(void)
{
  std::cout<<"usage: blaster [option]"<<std::endl;
  std::cout<<"Options are:"<<std::endl;
  std::cout<<"-q:\n\tquery bank name, default="<<bankQuery_name<<std::endl;
  std::cout<<"-s:\n\tsubject bank name, default="<<bank_name<<std::endl;
  std::cout<<"-a:\n\tall_by_all, default=FALSE"<<std::endl;
  std::cout<<"-X:\n\tuse NCBI-BLASTPLUS, default=First blast executable program found (blastplus, ncbi blast, wublast)"<<std::endl;
  std::cout<<"-N:\n\tuse NCBI-BLAST, default=First blast executable program found (blastplus, ncbi blast, wublast)"<<std::endl;
  std::cout<<"-W:\n\tuse WU-BLAST, default=First blast executable program found (blastplus, ncbi blast, wublast)"<<std::endl;
  std::cout<<"-n:\n\tblast name, default="<<type_blast<<std::endl;
  std::cout<<"-p:\n\tblast options, default="<<option_blast<<std::endl;
  std::cout<<"-l:\n\tfragment cutter length, default="<<cutter.getLength()<<std::endl;
  std::cout<<"-o:\n\tfragment cutter overlap, default="<<cutter.getOver()<<std::endl;
  std::cout<<"-w:\n\tminimum cutter word, default="<<cutter.getWord()<<std::endl;
  std::cout<<"-e:\n\textention of cutted file, default="<<cutter.getExtention()<<std::endl;
  std::cout<<"-S:\n\tsensitivity (low:0->high:4), default="<<sensitivity_lvl<<std::endl;
  std::cout<<"-E:\n\tE-value filter, default="<<eval_filter<<std::endl;
  std::cout<<"-I:\n\tidentity filter, default="<<id_filter<<std::endl;
  std::cout<<"-L:\n\tlength filter, default="<<len_filter<<std::endl;
  std::cout<<"-b:\n\toutput filename prefix, default="<<blaster_filename<<std::endl;
  std::cout<<"-B:\n\toutput filename prefix (no time stamp), default="<<blaster_filename<<std::endl;
  std::cout<<"-r:\n\tforce re-run all, default=Off"<<std::endl;
  std::cout<<"-P:\n\tonly prepare the banks, default=Off"<<std::endl;
  std::cout<<"-c:\n\tclean temporary files, default=Off"<<std::endl;
  std::cout<<"-v:\n\tverbose, default=0"<<std::endl;
}
//**************************************************************************
void BLRBlasterParameter::parseOptArg(int numarg, char* tabarg[])
{
  while(1)
    {
      if(numarg==1)
	{
	  help();
	  throw SDGException(NULL,"no parameter!",-1);
	}

// 	static struct option long_options[] =
// 	{
// 	  {"help",no_argument, 0, 'h'},
// 	  {"subject_bank",required_argument, 0, 's'},
// 	  {"query_bank",required_argument, 0, 'q'},
// 	  {"all_by_all",no_argument, 0, 'a'},
// 	  {"wu_blast",no_argument, 0, 'W'},
// 	  {"blast_plus",no_argument, 0, 'X'},
// 	  {"blast",no_argument, 0, 'N'},	
// 	  {"blast_prg_name",required_argument, 0, 'n'},
// 	  {"blast_parameter",required_argument, 0, 'p'},
// 	  {"cut_length",required_argument, 0, 'l'},
// 	  {"cut_overlap",required_argument, 0, 'o'},
// 	  {"cut_word",required_argument, 0, 'w'},
// 	  {"cut_extention",required_argument, 0, 'e'},
// 	  {"sensitivity_lvl",required_argument, 0, 'S'},
// 	  {"e-value_filter",required_argument, 0, 'E'},
// 	  {"identity_filter",required_argument, 0, 'I'},
// 	  {"length_filter",required_argument, 0, 'L'},
// 	  {"prefix_file",required_argument, 0, 'b'},
// 	  {"prefix_file_no_stamp",required_argument, 0, 'B'},
// 	  {"rerun",no_argument, 0, 'r'},
// 	  {"prepare",no_argument, 0, 'P'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;

// 	int c = getopt_long (numarg, tabarg,
// 			     "hs:q:aXNWn:p:l:o:w:e:I:S:E:L:b:B:rP",
// 			 long_options, &option_index);
	int c = getopt (numarg, tabarg,
			"hs:q:aXNWn:p:l:o:w:e:I:S:E:L:b:B:rPcv:");

	/* Detect the end of the options. */
	if (c == -1)
	  break;

	switch (c)
	  {
	  case 'h':
	    {
	      help();
	      exit(EXIT_SUCCESS);
	    }
	  case 'b':
	    {
	      blaster_filename=optarg;
	      time_stamp=true;
	      break;
	    }
	  case 'B':
	    {
	      blaster_filename=optarg;
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
	  case 'X':
	    {
	      is_ncbiBlastPlus=true;
	      break;
	    }
	  case 'N':
	    {
	      is_ncbiBlast=true;
	      break;
	    }
	  case 'W':
	    {
	      is_wuBlast=true;
	      break;
	    }
	  case 'a':
	    {
	      all_by_all=true;
	      break;
	    }
	  case 'n':
	    {
	      if (strcmp(optarg,"blastn")==0
		  || strcmp(optarg,"megablast")==0
		  || strcmp(optarg,"blastx")==0
		  || strcmp(optarg,"tblastn")==0
		  || strcmp(optarg,"blastp")==0
		  || strcmp(optarg,"tblastx")==0 )
		type_blast=optarg;
	      else
			{
			  std::cout<<"**"<<optarg
				  <<"** wrong choice for -n option, choose between:\n";
			  std::cout<<"   blastn"<<std::endl;
			  std::cout<<"   tblastx"<<std::endl;
			  std::cout<<"   megablast"<<std::endl;
			  throw SDGException(NULL,"parameter parse error",-1);
			}
	      break;
	    }
	  case 'p':
	    {
	      option_blast+=optarg;
	      option_blast+=" ";
	      break;
	    }
	  case 'l':
	    {
	      setLength(atoi(optarg));
	      break;
	    }
	  case 'o':
	    {
	      setOver(atoi(optarg));
	      break;
	    }
	  case 'w':
	    {
	      setWord(atoi(optarg));
	      break;
	    }
	  case 'e':
	    {
	      setExtention(optarg);
	      break;
	    }
	  case 'S':
	    {
	      sensitivity_lvl=atoi(optarg);
	      if(sensitivity_lvl>4)
		throw SDGException(NULL,"parse error: sensitivity value must be between 0 and 4",-1);
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
	  case 'r':
	    {
	      restart=true;
	      break;
	    }
	  case 'P':
	    {
	      prepare=true;
	      break;
	    }
	  case 'c':
	  {
		  cleanTmpFiles=true;
		  break;
	  }
	  case 'v':
	    {
	      verbose=atoi(optarg);
	      break;
	    }
	  case '?':
	  {
	    help();
	    throw SDGException(NULL,"parameter parse error",-1);
	    break;
	  }
	  default:
	    abort();
	  }
    }

    /* Print any remaining command line arguments (not options). */
    if (optind < numarg)
      {
	help();
	std::cout<<"non-option ARGV-elements: "<<std::endl;
	while (optind < numarg)
	  std::cout<<tabarg[optind++]<<std::endl;
	throw SDGException(NULL,"parameter parse error",-1);
      }

    if(bank_name=="<not set>")
      bank_name=bankQuery_name;
    if(bankQuery_name=="<not set>")
      bankQuery_name=bank_name;
}
//**************************************************************************
//------------------------
//! overload operator>>  -
//------------------------
std::istream& operator>>(std::istream& in, BLRBlasterParameter& para)
{
  char buff[256];
  para.reset();

  in.getline(buff,255);
  para.blaster_filename=buff;

  in.getline(buff,255);
  para.bank_name=buff;

  in.getline(buff,255);
  para.bankQuery_name=buff;

  in.getline(buff,255);
  para.type_blast=buff;

  in.getline(buff,255);
  para.option_blast=buff;

  in.getline(buff,255);
  para.bank_cut=buff;

  in.getline(buff,255);
  para.bankQuery_cut=buff;

  in.getline(buff,255);
  para.setLength(atol(buff));

  in.getline(buff,255);
  para.setOver(atoi(buff));

  in.getline(buff,255);
  para.setWord(atoi(buff));

  in.getline(buff,255);
  para.eval_filter=atof(buff);

  in.getline(buff,255);
  para.id_filter=atof(buff);

  in.getline(buff,255);
  para.len_filter=atoi(buff);

  in.getline(buff,255);
  para.numdo=atol(buff);

  in.getline(buff,255);
  para.nb_process=atoi(buff);

  in.getline(buff,255);
  if(buff[0]=='W')
    para.is_wuBlast=true;
    para.is_ncbiBlast=false;
    para.is_ncbiBlastPlus=false;  
  if(buff[0]=='N')
  	para.is_wuBlast=false;
    para.is_ncbiBlast=true;
    para.is_ncbiBlastPlus=false;  
  if(buff[0]=='X')
   	para.is_wuBlast=false;
    para.is_ncbiBlast=false;
    para.is_ncbiBlastPlus=true;
    
  return in;
}
//**************************************************************************
//------------------------
//! overload operator << -
//------------------------
std::ostream& operator<<(std::ostream& out, const BLRBlasterParameter& para)
{
  out<<para.blaster_filename<<std::endl;
  out<<para.bank_name<<std::endl;
  out<<para.bankQuery_name<<std::endl;
  out<<para.type_blast<<std::endl;
  out<<para.option_blast<<std::endl;
  out<<para.bank_cut<<std::endl;
  out<<para.bankQuery_cut<<std::endl;
  out<<para.getLength()<<std::endl;
  out<<para.getOver()<<std::endl;
  out<<para.getWord()<<std::endl;
  out<<para.eval_filter<<std::endl;
  out<<para.id_filter<<std::endl;
  out<<para.len_filter<<std::endl;
  out<<para.numdo<<std::endl;
  out<<para.nb_process<<std::endl;
  if(para.is_wuBlast)
    out<<"W"<<std::endl;
  if(para.is_ncbiBlast)
    out<<"N"<<std::endl;
  if(para.is_ncbiBlastPlus)
    out<<"X"<<std::endl;
  return out;
}
//**************************************************************************
//------------------------
//! overload operator == -
//------------------------
bool operator==(const BLRBlasterParameter& op1, const BLRBlasterParameter& op2)
{
  if ( op1.bank_cut==op2.bank_cut &&
       op1.bankQuery_cut==op2.bankQuery_cut &&
       op1.type_blast==op2.type_blast &&
       op1.option_blast==op2.option_blast)
    return true;
  else
    return false;
}
//**************************************************************************
//------------------------
//! overload operator != -
//------------------------
bool operator!=(const BLRBlasterParameter& op1,const BLRBlasterParameter& op2)
{
  if ( op1.bank_cut!=op2.bank_cut ||
       op1.bankQuery_cut!=op2.bankQuery_cut ||
       op1.type_blast!=op2.type_blast ||
       op1.option_blast!=op2.option_blast)
    return true;
  else
    return false;
}
