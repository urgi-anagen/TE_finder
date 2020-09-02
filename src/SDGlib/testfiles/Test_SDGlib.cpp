#include <SDGMemBioSeq.h>
#include <SDGFastaBioSeq.h>
#include <SDGFastaIstream.h>
#include <SDGBioSeqDB.h>
#include <cstring>
#include "Test_SDGlib.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_SDGlib);
//------------------------------------------------------------------------------------------------------------
void Test_SDGlib::test_SDGMemBioSeq( void )
{
    char str[]="ATATTTATTTTAGCGTTTACGCT";

    unsigned exp_len=std::strlen(str);
    
	SDGBioSeq seq=newSDGMemBioSeq(str);
	
	unsigned obs_len=seq.length();
	
	CPPUNIT_ASSERT_EQUAL(exp_len,obs_len);
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
