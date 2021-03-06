/**
 * \file BLRBlaster.h
 * \brief Header file for the class BLRBlaster
 */

#ifndef BLRBLASTER_H
#define BLRBLASTER_H

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <vector>
#include <list>
#include <map>

#include "SDGError.h"
#include "Reference.h"
#include "SDGString.h"
#include "SDGFastaIstream.h"
#include "SDGBioSeqDB.h"
#include "BLRBlasterParameter.h"
#include "BLRBlast.h"
#include "BLRNCBIBlast.h"
#include "BLRNCBIBlastPlus.h"
#include "BLRWUBlast.h"
#include "BLRMatchList.h"
#include "BLRBioSeqDB.h"
#include "RangeSeq.h"
#include "RangePair.h"

/**
 * \class BLRBlaster
 * \brief Class implementing BLASTER
 */
class BLRBlaster
{
 public:
  struct Key : std::pair<long,long>
  {
    Key(long i,long j)
    {
      first=i; second=j;
    };
  };

  typedef std::map<Key, std::list<RangePair> > MapAlign;

  BLRBlast *blasting;
  BLRWUBlast wubl;
  BLRNCBIBlast ncbibl;
  BLRNCBIBlastPlus ncbiblplus;
  BLRBlasterParameter para;

  unsigned min_len;

  unsigned count_treated;
  std::list<unsigned> waiting_seq;
  SDGString update_filename;

  void update(const std::list<unsigned>& nseq);

  SDGString listfilenamecut,bank_name,query_name;
  SDGBioSeqDB query_db,subject_db;
  BLRBioSeqDB queryCut_db,subjectCut_db;
  MapAlign map_align;
  MapAlign::iterator begin(){ return map_align.begin();};
  MapAlign::iterator end(){ return map_align.end();};
  std::map<std::string,long> name2numQ,name2numS;
  std::map<long,std::string> num2nameQ,num2nameS;
  bool same_db;

  void insert(RangePair& range);
  void add_raw(const BlastMatch& align);
  void glue( int verbose=0 );
  void clean_self( int verbose=0 );
  void save_align( int verbose=0 );
  void view_align(void);
  void load_raw( int verbose=0 );

public:

  /**
   * \fn BLRBlaster(BLRBlasterParameter& p)
   * \brief Constructor
   * \note initialize 'blasting', a pointer towards BLRNCBIBlast or BLRNCBIBlastPlus or BLRWUBlast
   * \param p an instance of BLRBlasterParameter
   */
  BLRBlaster(BLRBlasterParameter& p):wubl(p),ncbibl(p),ncbiblplus(p),para(p)
    {
      listfilenamecut=para.getBlasterFileName()+".raw";
      bank_name=para.getBankCut();
      query_name=para.getQuery();
      min_len=20;
      blasting = NULL;
      
      if(para.get_is_wuBlast())
      {
	  		blasting=&wubl;
	  }
      else if(para.get_is_ncbiBlastPlus())
      {
      		blasting=&ncbiblplus;
      }
      else if(para.get_is_ncbiBlast())
      {
        	blasting=&ncbibl;
      }
      else {
		  //in case no blasting program was selected , we try to auto find any available program (either blastplus, ncbi blast or wublast)
		  bool isBlastPlusInPath = system("makeblastdb >/dev/null 2>&1") == 256;
		  bool isWublastInPath = system("xdformat >/dev/null 2>&1") == 256;
		  bool isBlastInPath = system("formatdb >/dev/null 2>&1") == 512;
		  bool isWuBeforeBlastPlus = false;

		  if (isBlastPlusInPath and isWublastInPath)
		  {
			  FILE *p = popen("blastn -h", "r");
			  std::string output;
			  // char* buf; gcc warning -> change to following line
			  char buf[1024];
			  for (std::size_t count; (count = fread(buf, 1, sizeof(buf), p));)
				  output += std::string(buf, buf + count);
			  pclose(p); // 'output' variable now contains blastn -h output

			  std::transform(output.begin(), output.end(), output.begin(), ::tolower);

			  std::size_t found = output.find("washu");
			  if (found != std::string::npos)
			  {
				  isWuBeforeBlastPlus = true;
			  }
//			  std::system("blastn -h 2&> helpBlastn.log");
//			  std::ifstream file("helpBlastn.log");
//			  if (file)
//			  {
//				  std::stringstream buffer;
//				  buffer << file.rdbuf();
//				  file.close();
//				  remove("helpBlastn.log");
//
//				  std::string output;
//				  output = buffer.str();
//				  std::transform(output.begin(), output.end(), output.begin(), ::tolower);
//
//				  std::size_t found = output.find("washu");
//				  if (found != std::string::npos)
//				  {
//					  isWuBeforeBlastPlus = true;
//				  }
//			  }
//			  else {
//				  throw SDGException(NULL, "Fatal error, when trying to select blast executable", -1);
//			  }
		  }

		  if ((isBlastPlusInPath and !isWublastInPath) or (isBlastPlusInPath and isWublastInPath and !isWuBeforeBlastPlus)){
			  blasting=&ncbiblplus;
			  std::cout<<"Selecting NCBI Blastplus"<<std::endl;
		  }
		  else if ((!isBlastPlusInPath and isBlastInPath) or (isBlastPlusInPath and isWuBeforeBlastPlus and isBlastInPath)){
			  blasting=&ncbibl;
			  std::cout<<"Selecting NCBI Blast"<<std::endl;
		  }
		  else if ((!isBlastPlusInPath and !isBlastInPath and isWublastInPath) or (isBlastPlusInPath and isWuBeforeBlastPlus and !isBlastInPath)){
			  blasting=&wubl;
			  std::cout<<"Selecting WuBlast"<<std::endl;
		  }
		  else {
			  throw SDGException(NULL, "Fatal error, no blast executable found", -1);
		  }
      }
    };

  void run( int verbose=0 );

  void cleanTmpFiles( int verbose=0 );

};

#endif
