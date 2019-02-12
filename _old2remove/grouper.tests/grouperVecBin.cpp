/**
 * \file grouperVecBin.cpp
 */

#include <iostream>
#include <cstdlib>
#include "SDGError.h"
#include "Reference.h"
#include "SDGString.h"
#include "BlastMatch.h"
#include "../../grouper.threads/BLRGrouperParameter.h"
#include "BLRGroupVecBin.h"

int main( int argc, char* argv[] )
{
	try
	{
		std::cout<<"Beginning Grouper (version "<<VERSION<<")"<<std::endl;
		clock_t begin, end;
                double time_spent;
                begin = clock();

		BLRGrouperParameter para;
		para.parseOptArg( argc, argv );
		if(para.getVerbose()>0)
			para.view(std::cout);

		if(para.getVerbose()>0)
			std::cout<<"Initialization..."<<std::endl<<std::flush;
		BLRGroup groups( &para );
		if(para.getVerbose()>0)
			std::cout<<"Initialization was done."<<std::endl<<std::flush;

		// Build groups via simple-link clustering using percentage of overlap
		// between matches and merge matches if necessary
		groups.group( para.getVerbose() );
		if(para.getVerbose()>1)
			groups.show_group();

		// Select group containing at least n members
		// (new function 19.04.2006)
		if( para.getSizefilter() > 1 )
			groups.group_size_filter( para.getVerbose() );

		// Keep groups containing at least n unique members
		// (new function 25.08.2009)
		if( para.getIncludeFilter() > 1 )
			groups.include_filter( para.getVerbose() );

		// Build clusters of groups using length of overlap
		if( para.getGraphfilter() > 0 )
			groups.cluster( para.getVerbose() );

		// Save the results in files
		groups.save( para.getVerbose() );

		std::cout<<"End Grouper (version "<<VERSION<<")"<<std::endl;
		end = clock();
                time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
                std::cout<<"Time_spent****: "<<time_spent<<std::endl;

		exit( EXIT_SUCCESS );
	}

	catch( SDGException& e )
	{
		std::cout<<std::endl<<e.message;
		exit( EXIT_FAILURE );
	}
	catch(...)
	{
		std::cout<<std::endl<<"unknown exception catch";
		exit( EXIT_FAILURE );
	}
}
