#include <iostream>
#include <stdlib.h>
#include <SDGString.h>
#include <SDGBioSeq.h>
#include <SDGFastaBioSeq.h>
#include <SDGMemBioSeq.h>
#include <SDGSubBioSeq.h>
#include <SDGFastaIstream.h>
#include <SDGFastaOstream.h>
#include <SDGBioSeqDB.h>

int main(int argc, char* argv[])
{
	if (argc!=2)
		{
			std::cerr<<argv[0]<<" <sequence.fa> "<<std::endl;
			std::cerr<<"Pass a fasta sequence as argument!!"<<std::endl;
			exit(-1);
		}
  try{


    std::cout<<"Load SDGFastaBioSeq ..."<<std::endl;
    SDGBioSeq seq1=newSDGFastaBioSeq(argv[1]);
    std::cout<<seq1.getDE()<<std::endl;
    std::cout<<seq1.length()<<std::endl;
    SDGFastaOstream out1("fasta.fa");
    out1<<seq1;
    out1.close();

	std::cout<<"Load SDGMemBioSeq ..."<<std::endl;
	SDGFastaIstream in(argv[1]);
	SDGBioSeq seq2;
	in>>seq2;
	std::cout<<seq2.getDE()<<std::endl;
	std::cout<<seq2.length()<<std::endl;
	SDGFastaOstream out2("mem.fa");
	out2<<seq2;
	out2.close();

    std::cout<<"Load SDGSubBioSeq ..."<<std::endl;
    SDGBioSeq seq3=newSDGSubBioSeq(seq1,10,20);
    std::cout<<seq3.getDE()<<std::endl;
    std::cout<<seq3.length()<<std::endl;
    SDGFastaOstream out3("subfasta.fa");
    out3<<seq3;
    out3.close();
  }
  catch(SDGException e)
    {
      std::cout<<std::endl<<e.message;
    }
};






