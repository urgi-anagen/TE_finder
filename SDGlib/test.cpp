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
  try{
    std::cout<<"Load SDGMemBioSeq ..."<<std::endl;
    SDGFastaIstream in(argv[1]);
    SDGBioSeq seq2;
    in>>seq2;
    std::cout<<seq2.getDE()<<std::endl;
    std::cout<<seq2.length()<<std::endl;
    SDGFastaOstream out2("mem.fa");
    out2<<seq2;
    out2.close();
    //SDGFastaCout<<seq2<<std::flush;

    std::cout<<"Load SDGFastaBioSeq ..."<<std::endl;
    SDGBioSeq seq1=newSDGFastaBioSeq(argv[1]);
    std::cout<<seq1.getDE()<<std::endl;
    std::cout<<seq1.length()<<std::endl;
    SDGFastaOstream out1("fasta.fa");
    out1<<seq1;
    out1.close();
    //   SDGFastaCout<<seq1<<std::flush;
  }
  catch(SDGException e)
    {
      std::cout<<std::endl<<e.message;
    }
};






