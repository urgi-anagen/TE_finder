/*
 *
 * Definition d'une classe permettant de stocker un sequence biologique
 *
 *
 */

#include <SDGMemBioSeq.h>
#include <stdlib.h>
#include <ctype.h>
#include <iostream>
#include <sstream>
#include <SDGError.h>
#include <string.h>
#include <strings.h>
#include <math.h>



/*
 *
 * Fonctions membres de la classe sequence
 *
 *
 */



/*
 * 
 * Constructeurs de la classe
 *
 */

__SDGMemBioSeq::__SDGMemBioSeq(const __SDGMemBioSeq &ch):
  molecule(ch.molecule),
  forme(ch.forme),
  sequence(ch.sequence), 
  banque(ch.banque),
  access(ch.access),
  identificateur(ch.identificateur),
  definition(ch.definition)
{};


__SDGMemBioSeq::__SDGMemBioSeq(SDGString ch, type_molecule mol):
  forme(1),  
  sequence(ch),
  banque(), access(),
  identificateur(),
  definition()
{
  molecule= (mol) ?  mol:sequenceType(sequence.start());
  switch (molecule)
    {
    case PRT : sequence.setAlphabet(AAalphabet,'-');
               break;
    case NUC : sequence.setAlphabet(std::string("ACGTacgt"),'-');
               break;
    case NDG : sequence.setAlphabet(DNADGalphabet,'-');
               break;
    case NDF : break;
    };
}
	
SDGString __SDGMemBioSeq::toString(unsigned long d, unsigned long lg) const
{
  checkPos(d);
  if (forme)
    {
      if ((d+lg) > length()) {lg=0;};
      if (lg==0)            {lg=length()-d;};

      if (lg == 0) return ""; // Cas d'une sous sequence de longueur nulle
      if (( d == 0) && (lg == length()))
          return sequence; // Cas de la sequence complete
      return sequence.substr(d,lg);        // Sequence Lineaire
    }

  if (lg==0)            {lg=length()-d;};
  SDGString tmp;

  if (lg == 0) return tmp;  // Cas d'une sous sequence de longueur nulle
  if (( d == 0) && (lg == length())) return sequence; // Cas de la sequence complete

  long tmpLg = (d+lg > length()) ? length()-d : lg;
  
  tmp += sequence.substr(d,tmpLg);
  lg -= tmpLg;
  
  while (lg >0)
    {
      tmpLg = (lg > length()) ? length(): lg;
      tmp += sequence.substr(0,tmpLg);
      lg -= tmpLg;	  
    }

  return tmp;
}
     
SDGBioSeq __SDGMemBioSeq::complement() const
{
  std::ostringstream ostr;
  for(std::string::const_reverse_iterator i=sequence.rbegin();
      i!=sequence.rend();i++)
    ostr<<__SDGBioSeq::complement(*i);
  SDGBioSeq resultat=newSDGMemBioSeq(ostr.str());

  std::ostringstream def;
  def<<definition<<" re-oriented";
  resultat.setDE(def.str());
  //resultat.setDE(SDGString(definition) + SDGString(" re-oriented"));
  resultat.setAC(access);
  resultat.setID(identificateur);
  
 return resultat;
}



SDGBioSeq __SDGMemBioSeq::concat(SDGBioSeq seq) const
{
  return concat(seq.toString());
}

SDGBioSeq __SDGMemBioSeq::concat(SDGString seq) const
{
  SDGString new_seq(sequence);
  new_seq=new_seq+seq;
  SDGBioSeq resultat(__SDGMemBioSeq(new_seq,getMO()));

  resultat.setDE(definition);
  resultat.setAC(access);
  resultat.setID(identificateur);
  
  return resultat;
}

SDGBioSeq __SDGMemBioSeq::concat(char base) const
{
  char tmp[2];
  tmp[0]=base;
  tmp[1]=0;

  return concat(SDGString(tmp));
}

void *__SDGMemBioSeq::clone() const
{
  return (void*) new __SDGMemBioSeq(*this);
}










