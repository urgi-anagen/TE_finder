//
// Created by Hadi Quesneville on 31/05/2022.
//

#include "Test_FastExtAlign.h"
#include <BioSeq.h>

#include "Test_FragAlign.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_FastExtAlign);

//---------------------------------------------------------------------------------------------
void Test_FastExtAlign::test_align( void ) {
    std::ostringstream ostr1;
    ostr1<<"ATATTTATTTTAGCGTTTACGCTATGTGTTGCGTATTGCTAATCGCTATG";
    ostr1<<"TTACGCTATGTGTTATTTTTAGCGTTATTGCTAGCGTTTGCGATATTTAT";
    ostr1<<"ATATTTCGCGCTATGTGTTGCGATAGCGTTTATTATACCTATATCGCTAT";
    BioSeq seq1=BioSeq(ostr1.str());
    std::cout<<"\nSeq1 length:"<<seq1.size()<<std::endl;

    std::ostringstream ostr2;
    ostr2<<"ATATTTCGCGCTATGTGTTGCGATAGCGTTTATTATACCTATATCGCTAT";
    ostr2<<"TTACGCTATGTGTTATTTTTAGCGTTATTGCTAGCGTTTGCGATATTTAT";
    ostr2<<"ATATTTATTTTAGCGTTTACGCTATGTGTTGCGTATTGCTAATCGCTATG";
    BioSeq seq2=BioSeq(ostr2.str());
    std::cout<<"\nSeq2 length:"<<seq2.size()<<std::endl;

    FastExtAlign align;
    align.setSeq1(seq1);
    align.setSeq2(seq2);

    align.setStart(100,100,50);
    unsigned score=align.extend_dir(10);

    std::ostringstream ostr_obs;
    ostr_obs<<"("<<align.getEndSeq1()<<","<<align.getEndSeq2()<<") -> "<<score<<" / ";

    align.reset_align();
    align.setSeq1(seq1);
    align.setSeq2(seq2);

    align.setStart(50,50,50);
    score=align.extend_rev(10);

    ostr_obs<<"("<<align.getEndSeq1()<<","<<align.getEndSeq2()<<") -> "<<score;

    std::ostringstream ostr_exp;
    ostr_exp<<"(106,106) -> 16 / (42,43) -> 16";
    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
