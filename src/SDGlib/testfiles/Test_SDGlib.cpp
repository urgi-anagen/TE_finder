#include <SDGMemBioSeq.h>
#include <SDGFastaBioSeq.h>
#include <SDGFastaIstream.h>
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
void Test_SDGlib::test_SDGFastaBioSeq( void )
{
    char titre[]=">test seq";
    char str[]="ATATTTATTTTAGCGTTTACGCT";
    
    std::ofstream fout("test.fa");
    fout<<titre<<std::endl<<str<<std::endl;
    fout.close();
    
    std::ostringstream ostr_exp;
    ostr_exp<<str;
    
    SDGBioSeq seq=newSDGFastaBioSeq("test.fa");
    
    std::ostringstream ostr_obs;
    ostr_obs<<seq.toString();
 
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
/*
    std::cout<<"Load SDGSubBioSeq ..."<<std::endl;
    SDGBioSeq seq3=newSDGSubBioSeq(seq1,10,20);
    std::cout<<seq3.getDE()<<std::endl;
    std::cout<<seq3.length()<<std::endl;
    SDGFastaOstream out3("subfasta.fa");
    out3<<seq3;
    out3.close();
    */
