#ifndef HASHDNASEQ_H
#define HASHDNASEQ_H

#include <iostream>
#include <string>


struct HashDNASeq // Compute hashing value of a word
{
	unsigned kmer_size;
	unsigned hole_period;
	unsigned effectiveKmerSize;

	HashDNASeq(unsigned w,unsigned p=100): kmer_size(w), hole_period(p)
	{
		effectiveKmerSize=0;
		if(hole_period<2) hole_period=w+1;
		for(unsigned i=1;i<=kmer_size;i++)
		{
			if(i%hole_period!=0)
				effectiveKmerSize++;
		}
	  if(effectiveKmerSize>16)
	throw SDGException(NULL,"HashDNASeq: Effective word size must be <= 16 !!");
	};

	unsigned getEffectiveKmerSize(void){return effectiveKmerSize;};

	unsigned hash(const char* p)
	{

		unsigned h=0, val_nuc=0;
		for(unsigned i=0;i<kmer_size && *p!='\0';i++,p++)
		{
			if((i+1)%hole_period==0)
			{
				continue;
			}
			else
			{
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
				  }
				  default:
				  {
					val_nuc=0;
					break;
				  }
			  }
			  h<<=2;
			  h|=val_nuc;
			}
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
		for(unsigned i=0;i<kmer_size;i++)
		{
			if((kmer_size-i)%hole_period==0)
				{
					car='-';
				}
			else
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
			}
		  kmer.insert(kmer.begin(),car);
		}

		return kmer;
	};

	unsigned operator()(const char* p)
	   {return hash(p);};

};

#endif






