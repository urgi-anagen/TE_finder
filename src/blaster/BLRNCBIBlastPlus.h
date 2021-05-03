/**
 * \file BLRNCBIBlastPlus.h
 * \brief Header file for the class BLRNCBIBlastPlus
 */

#ifndef BLRNCBIBLASTPLUS_H
#define BLRNCBIBLASTPLUS_H

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
 * \class BLRNCBIBlastPlus
 * \brief Wrapper for the BLASTPLUS from the NCBI
 */
class BLRNCBIBlastPlus: public BLRBlast
{

 private:
  SDGString auto_blastparam;
  SDGString auto_blastnparam;
  SDGString auto_megablastparam;

 public:
  BLRNCBIBlastPlus(const BLRBlasterParameter& p):BLRBlast(p)
  {init(p);};

	 virtual ~BLRNCBIBlastPlus()
	 {};

	 void init(const BLRBlasterParameter& p)
	 {
		 para=p;
		 bank_name=para.getBankCut();
		 query_name=para.getQuery();

		 switch(para.getSensitivityLvl())
		 {
		 case 0:
		 {
			 auto_blastparam=" -evalue 0.001 -outfmt 6 -max_target_seqs 100000000";
			 auto_blastnparam=" -reward 1 -penalty -3 -gapopen 5 -gapextend 2 -word_size 11";
			 break;
		 }
		 case 1:
		 {
			 auto_blastparam=" -evalue 0.001 -outfmt 6 -max_target_seqs 100000000";
			 auto_blastnparam=" -reward 5 -penalty -4 -gapopen 10 -gapextend 6 -word_size 11";
			 break;
		 }
		 case 2:
		 {
			 auto_blastparam=" -evalue 0.001 -outfmt 6 -max_target_seqs 100000000";
			 auto_blastnparam=" -reward 10 -penalty -10 -gapopen 40 -gapextend 10 -word_size 8";
			 break;
		 }
		 case 3:
		 {
			 auto_blastparam=" -evalue 0.001 -outfmt 6 -max_target_seqs 100000000";
			 auto_blastnparam=" -reward 10 -penalty -10 -gapopen 30 -gapextend 10 -word_size 7";
			 break;
		 }
		 case 4:
		 {
			 auto_blastparam=" -evalue 0.001 -outfmt 6 -max_target_seqs 100000000";
			 auto_blastnparam=" -reward 10 -penalty -10 -gapopen 20 -gapextend 10 -word_size 7";
			 break;
		 }
		 default:
		 {
			 throw SDGException(NULL,"BLRNCBIBlastPlus.h: sensitivity value must be between 0 and 4",-1);
		 }
		 }

		 auto_megablastparam=" -word_size 11 -outfmt 6 -reward 1 -penalty -3 -dust no -xdrop_gap_final 20 -gapopen 5 -gapextend 2 -xdrop_ungap 10"; //TODO: see filter, xdrop off
	 };

	 void blast( int verbose=0 );

	 void pressdb( int verbose=0 );

	 Reference<BLRMatchList> getList()
	 {
		 return new BLRMatchList(alist);
	 };

};

#endif
