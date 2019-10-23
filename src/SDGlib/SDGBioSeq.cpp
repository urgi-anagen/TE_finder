/*
 * \file
 * \brief Definition d'une classe permettant de stocker un sequence biologique
 */

#include <cstdlib>
#include <cctype>
#include <string>

#include <SDGMemBioSeq.h>

const std::string __SDGBioSeq::AAalphabet = "ACDEFGHIKLMNPQRSTVWXY*acdefghiklmnpqrstvwxy";
const std::string __SDGBioSeq::DNADGalphabet = "ACGTMRSVWYHKDBNXacgtmrsvwyhkdbnx";
const std::string __SDGBioSeq::DNAalphabet = "ACGTacgt";

char __SDGBioSeq::complement(char a)
{
	const std::string DNAcomp = "TGCAKYSBWRDMHVNXtgcakysbwrdmhvnx";
	std::string::size_type rep = DNADGalphabet.find(a);
	if(rep==std::string::npos)
	{
		std::string msg="Can't find a complement to nucleotide: ";
		msg+="a";
		msg+=" !!\n";
		throw SDGBioSeq::NotDefinedException(NULL,msg);
	}
	return DNAcomp[rep];
}

__SDGBioSeq::type_molecule __SDGBioSeq::sequenceType(const char *seq)
{
	long A,C,G,T,U,total;

	A=C=G=T=U=total=0;

	while(*seq)
	{total++;
	switch ((*seq) & (255-32))
	{ case 'A' : A++; break;
	case 'C' : C++; break;
	case 'G' : G++; break;
	case 'T' : T++; break;
	case 'U' : U++; break;
	}
	seq++;
	}

	if ((A+C+G+T+U) == total)
		return NUC;
	else
		if ((A+C+G+T+U) > ((total*2)/5))
			return NDG;
		else
			return PRT;
}

void __SDGBioSeq::checkPos(unsigned long &p) const
{
	if( ( p != 0 ) && ( length() > 0 ) && ( p >= length() ) ) //to be adapted for circular sequences
		throw new SDGException(	this, "too large coordinate" );
}

bool __SDGBioSeq::isLinear() const
{
  return 1;
}

void __SDGBioSeq::setLinear(bool)
{
  throw new SDGBioSeq::NotDefinedException(this,"Function isLinear() not implemented");
}

SDGBioSeq  __SDGBioSeq::fivePrime(unsigned long pos, unsigned long lg) const
{
  throw new SDGBioSeq::NotDefinedException(this,"Function fivePrime() not implemented");
}

SDGBioSeq  __SDGBioSeq::threePrime(unsigned long pos, unsigned long lg) const
{
  throw new SDGBioSeq::NotDefinedException(this,"Function threePrime() not implemented");
}

unsigned long __SDGBioSeq::posInRef(unsigned long indice) const
{
  throw new SDGBioSeq::NotDefinedException(this,"Function posInRef() not implemented");
}

SDGBioSeq __SDGBioSeq::getBioseq() const
{
  throw new SDGBioSeq::NotDefinedException(this,"Function getBioseq() not implemented");
}

unsigned long __SDGBioSeq::posInAbsoluteRef(unsigned long indice) const
{
  throw new SDGBioSeq::NotDefinedException(this,"Function absoluteRef() not implemented");
}

SDGBioSeq __SDGBioSeq::getAbsoluteBioseq() const
{
  throw new SDGBioSeq::NotDefinedException(this,"Function getAbsoluteBioseq() not implemented");
}
