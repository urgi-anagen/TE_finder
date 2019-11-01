#include <string.h>
#include <ctype.h>
#include <SDGFastaBioSeq.h>
#include <SDGSubBioSeq.h>
#include <SDGMemBioSeq.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>

SDGString __SDGFastaBioSeq::litSeq(unsigned long from, unsigned long lg) const
{
  if(from>=buffseq_pos && from<=buffseq_pos+buffseq.length()-1 
     && lg<buffseq.length()-(from-buffseq_pos)
     && buffseq.length()!=0)
    return buffseq.substr(from-buffseq_pos,lg);
  else
    {
      ((__SDGFastaBioSeq*)this)->buffseq_pos=from;
      std::ifstream file(fichier.start());
      long decalage = (long)offset_sequence + (from/npl) + from;
      file.seekg(decalage);

      std::ostringstream buffstr;
      buffstr.setf(std::ios_base::uppercase);
      unsigned long count=0;
      unsigned long len=lg<buffseq_size?buffseq_size:lg;
      std::string line;
      while (file )
    	{ 
    	  if(std::getline(file,line))
    	  {
    	   count+=file.gcount()-1;
    	   buffstr<<line;
    	  }
    	  if((count < len) && file.peek()!='>') break;
     	}
      file.close();
      ((__SDGFastaBioSeq*)this)->buffseq=buffstr.str();

      return buffseq.substr(0,lg);
    }
}

char __SDGFastaBioSeq::charAt(unsigned long indice) const
{
  return litSeq(indice,1).charAt(0);
}

void __SDGFastaBioSeq::write_idx(std::ostream& out)
{
  out<<offset_sequence<<"\t"<<getDE()<<"\t"<<longueur<<std::endl;
}

void __SDGFastaBioSeq::read_idx(std::istream& in) 
{
	char buff[1024];

	in.getline(buff,1024,'\t');
	offset_sequence=atoi(buff);

	in.getline(buff,1024,'\t');
	setDE(buff);

	in.getline(buff,1024,'\n');
	longueur=atoi(buff);
}

__SDGFastaBioSeq::__SDGFastaBioSeq(SDGString fich, unsigned long decalage,
			       __SDGBioSeq::type_molecule mol,
				   unsigned buffsize)
{
  fichier=fich;
  std::ifstream file(fichier.start());
  forme =1;
  file.seekg(decalage);

  std::string line;
  SDGString titre;
    while (file && file.peek() != '>')
        file.get();
    if (std::getline(file, line)) {
        titre = SDGString(line).substr(1);
        titre = titre.trimR();
        offset_sequence = file.tellg();
    }
    setDE(titre);
    longueur = npl = 0;
    while (file) {
        if(file.peek() == '>')
            break;
        if (std::getline(file, line) )
            longueur += line.size();
        if (npl == 0)
            npl = longueur;
    }
/*
  buffseq="";
  buffseq_pos=0;
  buffseq_size=buffsize<longueur?buffsize:longueur;

  if (mol == __SDGBioSeq::NDF)
    {
      SDGString tmp = litSeq(0,buffseq_size);
      molecule = sequenceType(tmp.start());
      if(molecule==__SDGBioSeq::NUC)
	molecule=__SDGBioSeq::NDG;
    }
  else
    molecule=mol;*/
}
 
__SDGFastaBioSeq::__SDGFastaBioSeq(const __SDGFastaBioSeq &seq)
{
  buffseq=seq.buffseq;
  buffseq_pos=seq.buffseq_pos;
  buffseq_size=seq.buffseq_size;

  fichier=seq.fichier;
  offset_sequence = seq.offset_sequence;
  npl=seq.npl;
  molecule=seq.molecule;
  longueur=seq.longueur;
  banque = seq.banque;
  access = seq.access;
  identificateur = seq.identificateur;
  definition = seq.definition;
  forme = seq.forme;
}
 

SDGString __SDGFastaBioSeq::toString(unsigned long d, unsigned long lg) const
{
  checkPos(d);

  if (forme) // Sequence Lineaire
    {
      if ((d+lg) > length()) {lg=0;};
      if (lg==0)            {lg=length()-d;};

      return litSeq(d,lg);
    }
  std::cout<<"sequence non lineaire"<<std::endl;
  if (lg==0)            {lg=length()-d;};
  SDGString tmp;

  long tmpLg = (d+lg > length()) ? length()-d : lg;
  
  tmp += litSeq(d,tmpLg);
  lg -= tmpLg;
  
  while (lg >0)
    {
      tmpLg = (lg > length()) ? length(): lg;
      tmp += litSeq(0,tmpLg);
      lg -= tmpLg;	  
    }

  return tmp;
}

SDGBioSeq __SDGFastaBioSeq::subseq(unsigned long d, unsigned long lg)  const
{
  SDGBioSeq seq=newSDGMemBioSeq(toString(d,lg),getMO());
  seq.setDE(getDE());
  return seq   ;
}
  
SDGBioSeq __SDGFastaBioSeq::complement() const
{
  SDGBioSeq seq=newSDGMemBioSeq(toString(),getMO());
  seq.setDE(getDE());
  seq=seq.complement();
  return seq   ;
}

//  SDGBioSeq __SDGFastaBioSeq::subseq(unsigned long d, unsigned long lg) 
//  {
//    return newSDGSubBioSeq(*this,d,lg);
//  }
  
//  SDGBioSeq __SDGFastaBioSeq::complement()  const 
//  {
//    return newSDGSubBioSeq(*this,0,this->length(),0);
//  }
  
SDGBioSeq __SDGFastaBioSeq::concat(SDGBioSeq seq) const
{
  return concat(seq.toString());
}
 
SDGBioSeq __SDGFastaBioSeq::concat(SDGString seq) const
{
  SDGString str=toString();
  str=str+seq;
  __SDGMemBioSeq resultat(str,getMO());

 resultat.setDE(definition);
 resultat.setAC(access);
 resultat.setID(identificateur);

 return  SDGBioSeq(resultat);
}

SDGBioSeq __SDGFastaBioSeq::concat(char base) const
{
  char tmp[2];
  tmp[0]=base;
  tmp[1]=0;

  return concat(SDGString(tmp));
}

unsigned long __SDGFastaBioSeq::length() const
{
  return longueur;
}

void *__SDGFastaBioSeq::clone() const 
{
  return (void*) new __SDGFastaBioSeq(*this);
}

	// retourne des informations sur la sequence
SDGString __SDGFastaBioSeq::getBQ(void)  const 
{
  return banque;
}
SDGString __SDGFastaBioSeq::getAC(void)   const 
{
  return access;
}

SDGString __SDGFastaBioSeq::getID(void)   const 
{
  return identificateur;
}

SDGString __SDGFastaBioSeq::getDE(void)   const 
{
  return definition;
}

__SDGBioSeq::type_molecule  __SDGFastaBioSeq::getMO(void)   const 
{
  return molecule;
}

	// Permet de modifier les champs privees definissant la sequence   
void __SDGFastaBioSeq::setBQ( const SDGString new_BQ)
{
  banque = new_BQ;
}

void __SDGFastaBioSeq::setAC( const SDGString new_AC)
{
  access = new_AC;
}

void __SDGFastaBioSeq::setID( const SDGString new_ID)
{
  identificateur = new_ID;
}

void __SDGFastaBioSeq::setDE( const SDGString new_DE)
{
  definition = new_DE;
}





