/*
 * \file blaster2.cpp
 */

#include <string>
#include "BLRBlaster.h"
#include "BLRBlasterParameter.h"

int main( int argc, char* argv[] )
{
  std::string cmd = "date";
  try
    {
	  std::cout<<"Beginning Blaster (version "<<VERSION<<")"<<std::endl<<std::flush;
	  system( cmd.c_str() );

	  BLRBlasterParameter para;
      if (!para.start(argc, argv))
    	  throw SDGException(NULL,"parameter parse error",-1);

      BLRBlaster blaster(para);
      if(para.getVerbose()>0)
    	  std::cout<<"Run Blaster:"<<std::endl<<std::flush;
      blaster.run( para.getVerbose() );

      if(para.getCleanTmpFiles())
    	  blaster.cleanTmpFiles();

  		system( cmd.c_str() );
      std::cout<<"End Blaster (version "<<VERSION<<")"<<std::endl;
      exit(EXIT_SUCCESS);
    }
  catch(SDGException e)
    {
      std::cerr<<"BLASTER Exception catch: "<<e.message<<std::endl;;
      std::cerr<<"Exit with value "<<EXIT_FAILURE<<std::endl;
  		system( cmd.c_str() );
      exit(EXIT_FAILURE);
    }
  catch(...)
    {
      std::cerr<<"BLASTER Exception catch: unknown exception!!"<<std::endl;;
      std::cerr<<"Exit with value "<<EXIT_FAILURE<<std::endl;
  		system( cmd.c_str() );
      exit(EXIT_FAILURE);
    }
}
