/****
 *
 * grouper.cpp
 *
 *
 ***/

#include <iostream>
#include <cstdlib>
#include "SDGError.h"
#include "Reference.h"
#include "SDGString.h"
#include "BlastMatch.h"
#include "../../grouper.threads/BLRGrouperParameter.h"
#include "BLRGroup.h"

int main(int argc, char* argv[])
{
   try{
   std::cout<<"Beginning Grouper (version "<<VERSION<<")"<<std::endl;

    BLRGrouperParameter para;
    para.parseOptArg(argc,argv);
    if(para.getVerbose()>0)
    	para.view(std::cout);

    if(para.getVerbose()>0)
    	std::cout<<"Initialization..."<<std::endl<<std::flush;
    BLRGroup groups(&para);
    if(para.getVerbose()>0)
    	std::cout<<"Initialization was done."<<std::endl<<std::flush;

    if(para.getVerbose()>0)
    	std::cout<<"Group building..."<<std::endl<<std::flush;
    groups.group( para.getVerbose() );
    if(para.getVerbose()>0)
    	std::cout<<"Groups were built."<<std::endl<<std::flush;
    if(para.getVerbose()>1)
    	groups.show_group();

    // Select group containing at least n members (new function 19.04.2006)
    if( para.getSizefilter() > 1 )
      groups.group_size_filter( para.getVerbose() );
    if( para.getIncludeFilter() )
      groups.include_filter( para.getVerbose() );

	if( para.getGraphfilter() >= 0 )
	{
		if(para.getVerbose()>0)
			std::cout<<"Cluster building (length filter="<<para.getGraphfilter()<<")..."
			<<std::endl<<std::flush;
		groups.cluster( para.getVerbose() );
		if(para.getVerbose()>0)
			std::cout<<"Cluster were built."<<std::endl<<std::flush;
	}

	if(para.getVerbose()>0)
		std::cout<<"Write the results..."<<std::endl<<std::flush;
    groups.save( para.getVerbose()-1 );
    if(para.getVerbose()>1)
        std::cout<<"Writing group description... ";
    SDGString filename=para.getPrefixFileName()+".group.c"
      +SDGString(para.getCoverage())+".txt";
    std::ofstream out(filename);
    groups.show_group(out);
    out.close();
    if(para.getVerbose()>1)
    	std::cout<<" done"<<std::endl;
    if(para.getVerbose()>0)
    	std::cout<<"Results were written."<<std::endl<<std::flush;

    std::cout<<"End Grouper (version "<<VERSION<<")"<<std::endl;
    exit(EXIT_SUCCESS);
    }// end try

   catch(SDGException e)
    {
      std::cout<<std::endl<<e.message;
      exit(EXIT_FAILURE);
    }
   catch(...)
   {
	   std::cout<<std::endl<<"unknown exception catch";
	   exit(EXIT_FAILURE);
   }
}
