/*
 *
 *      Definit les fonction permettant la lecture d'un fichier
 *      de sequence au format fasta/pearson
 *
 *           si l'entete est de la forme
 *      >sp|X56794|ACCTHY ezea fdsf dsfsdfsdf dfssdfsd
 *
 *     la sequence lu est initialise correctement avec le numero d'acces
 *                                                     l'identificateur
 *                                                     le titre
 *     en absence de caractere | dans l'entete alors seul le titre est
 *     initialise
 * 
 */


#ifndef SDGFASTAISTREAM_H
#define SDGFASTAISTREAM_H

#include <fstream>
#include <SDGString.h>
#include <SDGMemBioSeq.h>
#include <SDGReference.h>

class SDGFastaIstream : public std::ifstream
{ 	
 public:

  SDGFastaIstream(void):std::ifstream() {};
  SDGFastaIstream(const SDGString &titre):std::ifstream(titre) {};
  SDGString name() const {return "SDGFastaIstream";};

  SDGFastaIstream &operator >> (SDGBioSeq &p);
};

extern SDGFastaIstream SDGFastaCin;



#endif   // SDGFASTAISTREAM_H

