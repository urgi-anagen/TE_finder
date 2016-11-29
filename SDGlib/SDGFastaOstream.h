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


#ifndef SDGFASTAOSTREAM_H
#define SDGFASTAOSTREAM_H

#include <fstream>

#include <SDGString.h>
#include <SDGMemBioSeq.h>

class SDGFastaOstream : public std::ofstream
{ 	
 public:

  SDGFastaOstream(void):std::ofstream()         {};
  SDGFastaOstream(SDGString titre):std::ofstream(titre) {};
  SDGString name() const {return "SDGFastaOstream";};

  SDGFastaOstream &operator << (SDGBioSeq p);

};

// Version de std::cout ecrivant les sequences au format fasta

//extern SDGFastaOstream SDGFastaCout;
//extern SDGFastaOstream SDGFastaCerr;

#endif   // SDGFASTAOSTREAM_H









