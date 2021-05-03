/**
 * \file BLRWUBlast.h
 * \brief Header file for the class BLRWUBlast
 */

#ifndef BLRWUBLAST_H
#define BLRWUBLAST_H

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
 * \class BLRWUBlast
 * \brief Wrapper for the BLAST from WU
 */
class BLRWUBlast: public BLRBlast
{

 private:
  SDGString auto_blastparam;
  SDGString auto_blastnparam;

 public:
	 BLRWUBlast(const BLRBlasterParameter& p):BLRBlast(p)
	 {init(p);};

	 virtual ~BLRWUBlast(){};

	 void init(const BLRBlasterParameter& p)
	 {
		 para=p;
		 bank_name=para.getBankCut();
		 query_name=para.getQuery();
		 switch(para.getSensitivityLvl())
		 {
		 case 0: //NCBI default
		 {
			 auto_blastparam=" warnings gapE2=0.001 V=0 gapW=76 evalues hspmax=0 gspmax=0 mformat=2 B=100000000";
			 auto_blastnparam=" M=1 N=-3 Q=5 R=2 gapall W=11";
			 break;
		 }
		 case 1: //WU default
		 {
			 auto_blastparam=" warnings gapE2=0.001 V=0 gapW=76 evalues hspmax=0 gspmax=0 mformat=2 B=100000000";
			 auto_blastnparam=" ";
			 break;
		 }
		 case 2:
		 {
			 auto_blastparam=" warnings gapE2=0.001 V=0 gapW=76 evalues hspmax=0 gspmax=0 mformat=2 B=100000000";
			 auto_blastnparam=" M=10 N=-12 Q=30 R=5 gapall S2=150 gapS2=220 S=220 X=220 gapX=450 W=8";
			 break;
		 }
		 case 3:
		 {
			 auto_blastparam=" warnings gapE2=0.001 V=0 gapW=76 evalues hspmax=0 gspmax=0 mformat=2 B=100000000";
			 auto_blastnparam=" M=10 N=-12 Q=25 R=5 gapall S2=100 gapS2=220 S=220 X=220 gapX=450 W=7";
			 break;
		 }
		 case 4:
		 {
			 auto_blastparam=" warnings gapE2=0.001 V=0 gapW=76 evalues hspmax=0 gspmax=0 mformat=2 B=100000000";
			 auto_blastnparam=" M=10 N=-12 Q=25 R=5 gapall S2=100 gapS2=200 S=200 X=200 gapX=450 W=6";
			 break;
		 }
		 default:
		 {
			 throw SDGException(NULL,"BLRWUBlast.h:sensitivity value must be between 0 and 4",-1);
		 }

		 }
	 };

	 void blast( int verbose=0 );

	 void pressdb( int verbose=0 );

	 Reference<BLRMatchList> getList()
	 {
		 return new BLRMatchList(alist);
	 };

};

#endif
