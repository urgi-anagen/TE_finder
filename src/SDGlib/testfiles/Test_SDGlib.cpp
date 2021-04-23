
#include <BioSeq.h>
#include <FastaIstream.h>
#include <FastaOstream.h>
#include <SDGFastaBioSeq.h>
#include <SDGFastaIstream.h>
#include <SDGBioSeqDB.h>
#include <cstring>
#include "Test_SDGlib.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_SDGlib);
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_BioSeq( void )
{

    std::ostringstream ostr_exp;
    ostr_exp<<"ATATTTATTTTAGCGTTTACGCT";
    BioSeq seq(ostr_exp.str(),"test");

    std::ostringstream ostr_obs;
    ostr_obs<<seq;

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_BioSeq_subseq( void )
{
    std::ostringstream ostr;
    ostr<<"ATATTTATTTTAGCGTTTACGCT";
    BioSeq seq(ostr.str(),"test");

    std::ostringstream ostr_exp;
    ostr_exp<<"TTTATTT";

    std::ostringstream ostr_obs;
    ostr_obs<<seq.subseq(3,7);

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_BioSeq_complement( void )
{
    std::ostringstream ostr;
    ostr<<"ATATTTATTTTAGCGTTTACGCT";
    BioSeq seq(ostr.str(),"test");

    std::ostringstream ostr_exp;
    ostr_exp<<"AGCGTAAACGCTAAAATAAATAT";

    std::ostringstream ostr_obs;
    ostr_obs<<seq.complement();

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_BioSeq_reverse( void )
{
    std::ostringstream ostr;
    ostr<<"ATATTTATTTTAGCGTTTACGCT";
    BioSeq seq(ostr.str(),"test");

    std::ostringstream ostr_exp;
    ostr_exp<<"TCGCATTTGCGATTTTATTTATA";

    std::ostringstream ostr_obs;
    ostr_obs<<seq.reverse();

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_FastaIstream( void )
{
    char titre1[]=">test seq1";
    char str11[]="ATATTTATTTTAGCGTTTACGCT";
    char str12[]="CAATCTTGTGCAAATGGCTGTGA";
    char str13[]="AG";

    char titre2[]=">test seq2";
    char str21[]="GCAGTTTTGGATTTAAGTGAATG";
    char str22[]="CTTCCTCTGGTGACTTTAGTAAG";
    char str23[]="CCCCATTTT";

    std::ofstream fout("test_in.fa");
    fout<<titre1<<std::endl
        <<str11<<std::endl
        <<str12<<std::endl
        <<str13<<std::endl;
    fout<<titre2<<std::endl
        <<str21<<std::endl
        <<str22<<std::endl
        <<str23<<std::endl;
    fout.close();

    std::ostringstream ostr_exp1;
    ostr_exp1<<str11<<str12<<str13;

    std::ostringstream ostr_exp2;
    ostr_exp2<<str21<<str22<<str23;

    FastaIstream in("test_in.fa");
    BioSeq seq1,seq2;
    in>>seq1;
    in>>seq2;

    std::ostringstream ostr_obs1,ostr_obs2;
    ostr_obs1<<seq1;
    ostr_obs2<<seq2;

    CPPUNIT_ASSERT_EQUAL(ostr_exp1.str(),ostr_obs1.str());
    CPPUNIT_ASSERT_EQUAL(ostr_exp2.str(),ostr_obs2.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_FastaOstream( void )
{
    char titre1[]="test seq1exp";
    char str1[]="ATATTTATTTTAGCGTTTACGCTCAATCTTGTGCAAATGGCTGTGAAGTTGGATTTAAGTGAATGCTTCTGGTGACTTTATGGATTTAAGTGAATGCTTCTGG";

    char titre2[]="test seq2exp";
    char str2[]="GCAGTTTTGGATTTAAGTGAATGCTTCCTCTGGTGACTTTAGTAAGCCCCATTTTTTAGCGTTTACGCTCAATCTTGTGCAAATGGCTGTGAAGTTGGATTTAAGTGAATGCTTCTG";

    BioSeq seq1exp(str1, titre1);
    BioSeq seq2exp(str2, titre2);

    FastaOstream out("test_out.fa");
    out << seq1exp;
    out << seq2exp;
    out.close();

    BioSeq seq1obs,seq2obs;
    FastaIstream in("test_out.fa");
    in>>seq1obs;
    in>>seq2obs;
    in.close();

    CPPUNIT_ASSERT_EQUAL(seq1exp, seq1obs);
    CPPUNIT_ASSERT_EQUAL(seq2exp, seq2obs);
}


//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_SDGMemBioSeq( void )
{

    std::ostringstream ostr_exp;
    ostr_exp<<"ATATTTATTTTAGCGTTTACGCT";
    SDGBioSeq seq=newSDGMemBioSeq(ostr_exp.str());

    std::ostringstream ostr_obs;
    ostr_obs<<seq.toString();

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_SDGMemBioSeq_subseq( void )
{
    std::ostringstream ostr;
    ostr<<"ATATTTATTTTAGCGTTTACGCT";
    SDGBioSeq seq=newSDGMemBioSeq(ostr.str());

    std::ostringstream ostr_exp;
    ostr_exp<<"TTTATTT";

    std::ostringstream ostr_obs;
    ostr_obs<<seq.subseq(3,7).toString();

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_SDGMemBioSeq_complement( void )
{
    std::ostringstream ostr;
    ostr<<"ATATTTATTTTAGCGTTTACGCT";
    SDGBioSeq seq=newSDGMemBioSeq(ostr.str());

    std::ostringstream ostr_exp;
    ostr_exp<<"AGCGTAAACGCTAAAATAAATAT";

    std::ostringstream ostr_obs;
    ostr_obs<<seq.complement().toString();

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_SDGMemBioSeq_reverse( void )
{
    std::ostringstream ostr;
    ostr<<"ATATTTATTTTAGCGTTTACGCT";
    SDGBioSeq seq=newSDGMemBioSeq(ostr.str());

    std::ostringstream ostr_exp;
    ostr_exp<<"TCGCATTTGCGATTTTATTTATA";

    std::ostringstream ostr_obs;
    ostr_obs<<seq.reverse().toString();

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_SDGFastaBioSeq( void )
{
    char titre[]=">test seq";
    char str[]="ATATTTATTTTAGCGTTTACGCTATG";
    
    std::ofstream fout("test.fa");
    fout<<titre<<std::endl<<str<<std::endl;
    fout.close();
    
    std::ostringstream ostr_exp;
    ostr_exp<<str<<" "<<strlen(str);
    
    SDGBioSeq seq=newSDGFastaBioSeq("test.fa");
    
    std::ostringstream ostr_obs;
    ostr_obs<<seq.toString()<<" "<<seq.length();

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_SDGFastaIstream( void )
{
    char titre1[]=">test seq1";
    char str11[]="ATATTTATTTTAGCGTTTACGCT";
    char str12[]="CAATCTTGTGCAAATGGCTGTGA";
    char str13[]="AG";

    char titre2[]=">test seq2";
    char str21[]="GCAGTTTTGGATTTAAGTGAATG";
    char str22[]="CTTCCTCTGGTGACTTTAGTAAG";
    char str23[]="CCCCATTTT";

    std::ofstream fout("test.fa");
    fout<<titre1<<std::endl
        <<str11<<std::endl
        <<str12<<std::endl
        <<str13<<std::endl;
    fout<<titre2<<std::endl
        <<str21<<std::endl
        <<str22<<std::endl
        <<str23<<std::endl;
    fout.close();

    std::ostringstream ostr_exp1;
    ostr_exp1<<str11<<str12<<str13;

    std::ostringstream ostr_exp2;
    ostr_exp2<<str21<<str22<<str23;

    SDGFastaIstream in("test.fa");
    SDGBioSeq seq1,seq2;
    in>>seq1;
    in>>seq2;

    std::ostringstream ostr_obs1,ostr_obs2;
    ostr_obs1<<seq1.toString();
    ostr_obs2<<seq2.toString();

    CPPUNIT_ASSERT_EQUAL(ostr_exp1.str(),ostr_obs1.str());
    CPPUNIT_ASSERT_EQUAL(ostr_exp2.str(),ostr_obs2.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_SDGBioSeqDB( void )
{

    char titre1[]=">test seq1";
    char str11[]="ATATTTATTTTAGCGTTTACGCT";
    char str12[]="CAATCTTGTGCAAATGGCTGTGA";
    char str13[]="AG";

    char titre2[]=">test seq2";
    char str21[]="GCAGTTTTGGATTTAAGTGAATG";
    char str22[]="CTTCCTCTGGTGACTTTAGTAAG";
    char str23[]="CCCCATTTT";

    std::ofstream fout("test.fa");
    fout<<titre1<<std::endl
        <<str11<<std::endl
        <<str12<<std::endl
        <<str13<<std::endl;
    fout<<titre2<<std::endl
        <<str21<<std::endl
        <<str22<<std::endl
        <<str23<<std::endl;
    fout.close();

    std::ostringstream ostr_exp1;
    ostr_exp1<<str11<<str12<<str13;

    std::ostringstream ostr_exp2;
    ostr_exp2<<str21<<str22<<str23;

    int verbose=1;
    SDGBioSeqDB db("test.fa",verbose);
    std::vector<SDGBioSeq> vSeq;
    for (SDGBioSeqDB::iterator db_it = db.begin(); db_it != db.end(); db_it++) {

        SDGBioSeq seq = (*db_it);
        vSeq.push_back(seq);
    }


    std::ostringstream ostr_obs1,ostr_obs2;
    ostr_obs1<<vSeq[0].toString();
    ostr_obs2<<vSeq[1].toString();

    CPPUNIT_ASSERT_EQUAL(ostr_exp1.str(),ostr_obs1.str());
    CPPUNIT_ASSERT_EQUAL(ostr_exp2.str(),ostr_obs2.str());
}
