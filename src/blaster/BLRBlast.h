/**
 * \file BLRBlast.h
 * \brief Header file for the class BLRBlast
 */

#ifndef BLRBLAST_H
#define BLRBLAST_H

#include <iostream>
#include <stdlib.h>
#include "Reference.h"
#include "BioSeqDB.h"
#include "BioSeq.h"
#include "SDGFastaOstream.h"
#include "FastaOstream.h"
#include "SDGString.h"
#include "BLRMatchList.h"
#include "BLRBlasterParameter.h"

/**
 * \class BLRBlast
 * \brief Interface for BLRNCBIBlast and BLRWUBlast
 */
class BLRBlast
{

 protected:
  SDGString bank_name,query_name;
  BLRBlasterParameter para;
  SDGString query_filename;
  SDGString result_filename;

 public:
  BLRBlast(const BLRBlasterParameter& p) : para(p){};

  virtual ~BLRBlast(){};

  /**
   * \fn virtual void init(const BLRBlasterParameter& p)=0;
   * \brief set the BLAST parameters (banks and scoring scheme)
   */
  virtual void init(const BLRBlasterParameter& p)=0;

  /**
   * \fn void prepblast(const std::list<SDGBioSeq>& lquery,
   * \brief save a list of queries in a fasta file
   * \param lquery list of SDGBioSeq instances
   * \param first_num_seq number of the first sequence
   */
  void prepblast(const std::list<unsigned>& batch_num_seq, const BioSeqDB& query_db,
		  unsigned first_num_seq)
  {
    query_filename=query_name+"_query"+SDGString(first_num_seq)+".fa";
    FastaOstream fastout(query_filename);
    BioSeq s;
    for(std::list<unsigned>::const_iterator i=batch_num_seq.begin();
	i!=batch_num_seq.end(); i++){
        s = query_db[(*i) - 1];
        fastout<<s;
    }

    fastout.close();

    result_filename=query_name+"_"+bank_name.afterlast("/")
      +"_"+SDGString(first_num_seq)+".res";
  };

  /**
   * \fn virtual void blast( int verbose=0 )=0;
   * \brief launch the given BLAST executable on the cut, indexed banks
   */
  virtual void blast( int verbose=0 )=0;

  /**
   * \fn virtual void pressdb( int verbose=0 )=0;
   * \brief check if the cut bank is indexed, otherwise index it
   */
  virtual void pressdb( int verbose=0 )=0;

};

#endif
