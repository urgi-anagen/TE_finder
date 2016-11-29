#ifndef HASHDNASEQ_H
#define HASHDNASEQ_H

#include <string>

struct HashDNASeq // Compute hashing value of a word
{
	unsigned word_size;
	unsigned msk;
	unsigned lb;
	unsigned hb;
	HashDNASeq(unsigned w,unsigned m=0): word_size(w), msk(m)
	{
		unsigned halfword=word_size/2;
		unsigned halfkmsk=msk/2;
		if(msk!=0 && halfkmsk==0)
			halfkmsk=1;
		if(halfword>halfkmsk)
		{
			lb=halfword-halfkmsk;
			hb=halfword+halfkmsk;
		}
		else
		{
			lb=hb=0;
		}

	  if(w>16)
	throw SDGException(NULL,"HashDNASeq: Word size must be <= 16 !!");
	};

	unsigned getEffectiveKmerSize(void){return word_size-(hb-lb);};

	unsigned hash(const char* p)
	{

		unsigned h=0;
		for(unsigned i=0;i<word_size && *p!='\0';i++)
		{
		  if(i>lb && i<=hb) continue;
		  h<<=2;
		  unsigned val_nuc=0;
		  switch(*p)
		  {
			  case 'C':
			  {
				  val_nuc=1;
				  break;
			  }
			  case 'G':
			  {
				  val_nuc=2;
				  break;
			  }
			  case 'T':
			  {
				  val_nuc=3;
				  break;
			  }
			  case 'c':
			  {
				  val_nuc=1;
				  break;
			  }
			  case 'g':
			  {
				  val_nuc=2;
				  break;
			  }
			  case 't':
			  {
				  val_nuc=3;
				  break;
			  }		  	  default:
			  {
				val_nuc=0;
				break;
			  }
		  }

		  h|=val_nuc;
		  p++;
		}
	  return h;
	};

	std::string reverse_hash(unsigned key)
	{
		unsigned maskw=1;
		maskw<<=2;
		maskw--;

		char car;
		std::string kmer;

		unsigned val=0;
		for(unsigned i=0;i<word_size;i++)
		{
		  val=key&maskw;
		  key>>=2;
		  switch(val)
		  {
			  case 1:
			  {
				  car='C';
				  break;
			  }
			  case 2:
			  {
				  car='G';
				  break;
			  }
			  case 3:
			  {
				  car='T';
				  break;
			  }

			  default:
			  {
				car='A';
				break;
			  }
		  }
		  kmer.insert(kmer.begin(),car);
		}

		return kmer;
	};

	unsigned operator()(const char* p)
	   {return hash(p);};

};

#endif






