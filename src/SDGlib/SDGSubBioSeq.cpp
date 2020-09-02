#include <SDGSubBioSeq.h>

__SDGSubBioSeq::__SDGSubBioSeq(SDGBioSeq &seq, unsigned long d, 
		unsigned long lg)
{
  debut.push_back(d);

  if (lg == 0) lg = seq.length() - d;

  longueur.push_back(lg);

  bioseq = seq;
  molecule       = seq.getMO();
  banque         = seq.getBQ();
  definition     = seq.getDE();
  access         = seq.getAC();
  identificateur = seq.getID();
}
__SDGSubBioSeq::__SDGSubBioSeq(const __SDGBioSeq &seq, unsigned long d, 
		unsigned long lg)
{
  debut.push_back(d);

  if (lg == 0) lg = seq.length() - d;

  longueur.push_back(lg);


  bioseq=seq;

  molecule       = seq.getMO();
  banque         = seq.getBQ();
  definition     = seq.getDE();
  access         = seq.getAC();
  identificateur = seq.getID();
}

__SDGSubBioSeq::__SDGSubBioSeq(const __SDGSubBioSeq &seq)
{
  
  bioseq = seq.bioseq;

  molecule = seq.getMO();
  banque = seq.banque;
  access = seq.access;
  identificateur = seq.identificateur;
  definition = seq.definition;

  debut = seq.debut;
  longueur = seq.longueur;

}
 
char __SDGSubBioSeq::charAt(unsigned long indice) const
{
  char tmp = 0;
  checkPos(indice);

  unsigned long lt = 0;

    for (size_t i=0; (i < longueur.size()) && (tmp == 0); i++)
    {
      if ( (indice >=lt) && (indice < lt+longueur[i]) )
	{
	  tmp = bioseq.charAt(debut[i]+indice - lt);
	}

      lt+=longueur[i];      
      
    }

    return tmp;
} 

SDGString __SDGSubBioSeq::toString(unsigned long d, unsigned long lg) const
{
  SDGString sortie;

  for (size_t i=0; i < debut.size(); i++)
         sortie += bioseq.toString(debut[i],longueur[i]);

  return sortie.substr(d,lg);
}



SDGBioSeq __SDGSubBioSeq::subseq(unsigned long d, unsigned long lg)  const
{
  unsigned long lt=0;
  char start=0;

  checkPos(d);              // Verifie la validite de la position de depart

  if (lg == 0)             // Si lg = 0 (parametre par defaut
    lg = length() - d;      //   |--> calcule lg pour le reste de la sequence

  __SDGSubBioSeq tmp(*this);

  tmp.debut.erase(tmp.debut.begin(),tmp.debut.end());
  tmp.longueur.erase(tmp.longueur.begin(),tmp.longueur.end());

  for (size_t i=0; i < longueur.size(); i++)
    {
      lt += longueur[i];
      if (start)
		if (lg >= lt)
		  {
			tmp.debut.push_back(debut[i]);
			tmp.longueur.push_back(longueur[i]);
		  }
		else
		  {
			start=0;
			tmp.debut.push_back(debut[i]);
			tmp.longueur.push_back(lg-lt+longueur[i]);
			break;
		  }
      else
		if (d < lt)
		  {
			start=1;
			tmp.debut.push_back(debut[i]+d);
			tmp.longueur.push_back(longueur[i]-d);
			if (tmp.longueur[0] > lg)
			  {
				tmp.longueur[0] = lg;
				start=0;
				break;
			  }
			lt = tmp.longueur[0];
		  }
    }

  tmp.setDE(definition);
  tmp.setAC(access);
  tmp.setID(identificateur);

  return SDGBioSeq(tmp);
      
}
  
SDGBioSeq __SDGSubBioSeq::complement() const
{
  __SDGSubBioSeq tmp(*this);
  return SDGBioSeq(tmp.complement());
}

SDGBioSeq __SDGSubBioSeq::reverse() const
{
    __SDGSubBioSeq tmp(*this);
    return SDGBioSeq(tmp.reverse());
}

SDGBioSeq __SDGSubBioSeq::concat(SDGBioSeq seq) const
{
  return concat(seq.toString());
}
  
SDGBioSeq __SDGSubBioSeq::concat(SDGString seq) const
{
  SDGString str=toString();
  str+=seq;
  __SDGMemBioSeq resultat(str,getMO());

 resultat.setDE(definition);
 resultat.setAC(access);
 resultat.setID(identificateur);

 return SDGBioSeq(resultat);
}

SDGBioSeq __SDGSubBioSeq::concat(char base) const
{
  char tmp[2];
  tmp[0]=base;
  tmp[1]=0;

  return concat(SDGString(tmp));
}

unsigned long __SDGSubBioSeq::length() const
{
  long tmp=0;

  for (size_t i=0; i < longueur.size(); i++)
    tmp+=longueur[i];

  return tmp;  
}

void *__SDGSubBioSeq::clone() const 
{
  return (void*) new __SDGSubBioSeq(*this);
}

	// retourne des informations sur la sequence
SDGString __SDGSubBioSeq::getBQ(void) const 
{
  return banque;
}
SDGString __SDGSubBioSeq::getAC(void) const
{
  return access;
}

SDGString __SDGSubBioSeq::getID(void) const
{
  return identificateur;
}

SDGString __SDGSubBioSeq::getDE(void) const
{
  return definition;
}

__SDGBioSeq::type_molecule __SDGSubBioSeq::getMO(void) const
{
  return molecule;
}

	// Permet de modifier les champs privees definissant la sequence   
void __SDGSubBioSeq::setBQ(const SDGString new_BQ)
{
  banque = new_BQ;
}

void __SDGSubBioSeq::setAC(const SDGString new_AC)
{
  access = new_AC;
}

void __SDGSubBioSeq::setID(const SDGString new_ID)
{
  identificateur = new_ID;
}

void __SDGSubBioSeq::setDE(const SDGString new_DE)
{
  definition = new_DE;
}

unsigned long __SDGSubBioSeq::posInRef(unsigned long indice) const
{
  long tmp = 0;
  checkPos(indice);
  

  unsigned long lt = 0;

  for (size_t i=0; (i < longueur.size()) && (tmp == 0); i++)
    {
      if ( (indice >=lt) && (indice < lt+longueur[i]) )
		{
		  tmp = debut[i]+indice - lt;
		}
      lt+=longueur[i];      
    }
    return tmp;
}
SDGBioSeq __SDGSubBioSeq::getBioseq() const
{
  return bioseq;
}
unsigned long __SDGSubBioSeq::posInAbsoluteRef(unsigned long indice) const
{
  if (bioseq.name() == "SDGSubBioSeq")
    return bioseq.posInAbsoluteRef(posInRef(indice));
  else
    return posInRef(indice);

}
SDGBioSeq __SDGSubBioSeq::getAbsoluteBioseq() const
{
    if (bioseq.name() == "SDGSubBioSeq")
      return bioseq.getAbsoluteBioseq();
    else
      return bioseq;
}



