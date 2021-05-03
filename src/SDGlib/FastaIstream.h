/*
 *
 *      Definit les fonction permettant la lecture d'un fichier
 *      de sequence au format fasta/pearson
 */


#ifndef FASTAISTREAM_H
#define FASTAISTREAM_H

#include <fstream>
#include <BioSeq.h>


class FastaIstream : public std::ifstream
{ 	
 public:

  FastaIstream(void):std::ifstream() {};
  FastaIstream(const std::string &titre):std::ifstream(titre) {};

  FastaIstream &operator >> (BioSeq &p);
};




#endif   // SDGFASTAISTREAM_H

