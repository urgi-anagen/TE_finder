/**
 * \file BLRNCBIBlast.h
 * \brief Header file for the class BLRNCBIBlast
 */

#ifndef BLRNCBIBLAST_H
#define BLRNCBIBLAST_H

#include <iostream>
#include <stdlib.h>
#include <Reference.h>
#include <SDGMemBioSeq.h>
#include <SDGFastaOstream.h>
#include <SDGString.h>

#include "BLRBlast.h"
#include "BLRMatchList.h"
#include "BLRBlasterParameter.h"

/**
 * \class BLRNCBIBlast
 * \brief Wrapper for the BLAST from the NCBI
 */
class BLRNCBIBlast: public BLRBlast
{

 private:
  SDGString auto_blastparam;
  SDGString auto_blastnparam;
  SDGString auto_megablastparam;

 public:
  BLRNCBIBlast(const BLRBlasterParameter& p):BLRBlast(p)
  {init(p);};

	 virtual ~BLRNCBIBlast()
	 {};

	 void init(const BLRBlasterParameter& p)
	 {
		 para=p;
		 bank_name=para.getBankCut();
		 query_name=para.getQuery();

		 switch(para.getSensitivityLvl())
		 {
		 case 0: //NCBI default
		 {
			 auto_blastparam=" -e 0.001 -m 8 -v 0 -b 100000000";
			 auto_blastnparam=" ";
			 break;
		 }
		 case 1: //WU default
		 {
			 auto_blastparam=" -e 0.001 -m 8 -v 0 -b 100000000";
			 auto_blastnparam=" -r 5 -q -4 -G 10 -E 6 -g T";
			 break;
		 }
		 case 2:
		 {
			 auto_blastparam=" -e 0.001 -m 8 -v 0 -b 100000000";
			 auto_blastnparam=" -r 10 -q -10 -G 40 -E 10 -g T -W 8";
			 break;
		 }
		 case 3:
		 {
			 auto_blastparam=" -e 0.001 -m 8 -v 0 -b 100000000";
			 auto_blastnparam=" -r 10 -q -10 -G 30 -E 10 -g T -W 7";
			 break;
		 }
		 case 4:
		 {
			 auto_blastparam=" -e 0.001 -m 8 -v 0 -b 100000000";
			 auto_blastnparam=" -r 10 -q -10 -G 20 -E 10 -g T -W 7";
			 break;
		 }
		 default:
		 {
			 throw SDGException(NULL,"BLRNCBIBlast.h: sensitivity value must be between 0 and 4",-1);
		 }
		 }

		 auto_megablastparam=" -W 12 -D 3";
	 };

	 void blast( int verbose=0 );

	 void pressdb( int verbose=0 );

	 Reference<BLRMatchList> getList()
	 {
		 return new BLRMatchList(alist);
	 };

};

#endif
