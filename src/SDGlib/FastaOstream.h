/*
 *
 *      Definit les fonction permettant la lecture d'un fichier
 *      de sequence au format fasta/pearson
 *
 */


#ifndef FASTAOSTREAM_H
#define FASTAOSTREAM_H

#include <fstream>
#include <string>
#include <BioSeq.h>

class FastaOstream : public std::ofstream
{ 	
 public:

  FastaOstream(void):std::ofstream()         {};
  FastaOstream(std::string titre):std::ofstream(titre) {};

  FastaOstream &operator << (BioSeq p);

};

#endif   // SDGFASTAOSTREAM_H









