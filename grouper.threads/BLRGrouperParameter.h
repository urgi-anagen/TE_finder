/**
 * \file BLRGrouperParameter.h
 * \brief Header file for the class BLRGrouperParameter
 */

#ifndef BLRGROUPERPARAMETER_H
#define BLRGROUPERPARAMETER_H

#include <cstdlib>
#include <iostream>
#include <time.h>
#include <fstream>
#include <SDGString.h>
#include <list>
#include "BLRJoinParameter.h"

/**
 * \class BLRGrouperParameter
 * \brief Parameters for the GROUPER program
 */
class BLRGrouperParameter: public BLRJoinParameter
{	
	private:
	    double coverage;
		int graphfilter;
		unsigned includefilter;
		unsigned sizefilter;
		unsigned nb_sets;
		int verbose;

	public:
		//! -Constructor
		BLRGrouperParameter( void ){ reset(); };

		void reset( void )
		{
			BLRJoinParameter::reset();
			gap_pen=1;
      		dist_pen=10;
			coverage=0.95;
			graphfilter=0;
			sizefilter=1;
			includefilter=1;
			verbose=0;
			nb_sets=1;
			loadPath=false;
		};

		//! -Destructor
		virtual ~BLRGrouperParameter(void){};

	//! Methods
	public:
		void parseOptArg( int numarg, char *tabarg[] );
		void help( void );
		void view( std::ostream& out ) const;
		double getOverlap( void ) const { return overlap; };
		double getCoverage( void ) const { return coverage; };
		int getGraphfilter( void ) const { return graphfilter; };
		unsigned getIncludeFilter(void) const { return includefilter; };
		unsigned getSizefilter( void ) const { return sizefilter; };
		void setSizeFiler(unsigned sizeFilterParam){ sizefilter = sizeFilterParam; };
		unsigned getNbSets(void) {return nb_sets;}
		int getVerbose( void ) { return verbose; };
};

const static BLRGrouperParameter defaultBLRGrouperParameters;
#endif
