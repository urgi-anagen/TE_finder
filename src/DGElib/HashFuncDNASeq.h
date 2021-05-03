#ifndef HASHFUNCDNASEQ_H
#define HASHFUNCDNASEQ_H

#include <iostream>
#include <string>

struct HashFuncDNASeq // Compute hashing value of a word
{
	unsigned kmer_size;
	unsigned hole_period;
    unsigned hole_size;
	unsigned effectiveKmerSize;
	std::string mask;

	HashFuncDNASeq(unsigned w, unsigned p=0, unsigned s=1): kmer_size(w), hole_period(p), hole_size(s)
	{
        if(hole_period<2) {
            hole_period=w+1;
            hole_size=0;
        }
        mask.resize(kmer_size) ;
        effectiveKmerSize=0;
        build_mask_spaced_hole();
        if(effectiveKmerSize>15)
            throw SDGException(NULL,"HashDNASeq: Effective word size must be < 16 !!");
	};

    virtual ~HashFuncDNASeq() {
    }

    void build_mask_spaced_hole(void){
        unsigned i=0;
        while(i<kmer_size)
        {
            if((i+1)%hole_period==0)
            {
                for(unsigned j=0;j<hole_size && i<kmer_size;j++,i++){
                    mask[i]='-';
                }
            }
            else
            {
                mask[i]='+';
                effectiveKmerSize++;
                i++;
            }
        }
    }
    unsigned getEffectiveKmerSize(void){return effectiveKmerSize;};
    std::string getMask(void){return std::string(mask);};

	unsigned hash(const char* p)
	{

		unsigned h=0, val_nuc=0;
		for(unsigned i=0;i<kmer_size && *p!='\0';i++,p++)
		{
			if(mask[i]=='-')
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
		for(unsigned i=1;i<=kmer_size;i++)
		{
			if(mask[kmer_size-i]=='-')
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






