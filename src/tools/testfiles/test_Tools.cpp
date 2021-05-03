#include <stdio.h>
#include "SDGString.h"
#include "Cutter.h"
#include "test_Tools.h"



CPPUNIT_TEST_SUITE_REGISTRATION(test_Tools);
//----------------------------------------------------------------------------------------------------------
void test_Tools::initTestSeq(void){

    //build fasta file with 1 seq
    const char *text1 = R""""(>test seq
ATATTTATTTTAGCGTTTACGCTATG)"""";
    std::ofstream fout("test1seq.fa");
    fout<<text1<<std::endl;
    fout.close();

    //build fasta file with 2 seq
    const char *text2 = R""""(>test seq1
ATATTTATTTTAGCGTTTACGCT
CAATCTTGTGCAAATGGCTGTGA
AG
>test seq2
GCAGTTTTGGATTTAAGTGAATG
CTTCCTCTGGTGACTTTAGTAAG
CCCCATTTT)"""";
    fout.open("test2seq.fa");
    fout<<text2<<std::endl;
    fout.close();
}
//----------------------------------------------------------------------------------------------------------
void test_Tools::cutterDB(void){

    std::string filename("testCutterDB.fa");
    const char *text_in = R""""(>seq1
AGTTACCATGCCCAGCATTAACCCCCCTCAACAACCACCTCCGCCTATGAAGCCCGCCCG
GTAGGCGACATCAGCAAAGTGCCAACGCTGTATATATATATATTGATCACGAGCTACCAT
GCCAGCATAGCCTCGTCCCCCGCTACCTGAAACTCTGTTGCACCCGATGATGAAATCGGC
ATGAACACACACACACACACATATTCACACATACTGCGTGGTGGCGACGCTCTCGACGTT
GAAATTAACATGAACCTACACACACACATACACACTGCGTAATCGCGACGTCCTCGTCTA
GACTCGCTATCTAGAACTAACGGGATCAAAAGCACTGCTGCTTGCCCGTGCGTATACATT
AAGAATAAAGCTTTCATCATTCTTGATCTTGACACCAAACCGAGCAGTTGATTTATTTAA
AGTGGCAAATATATATAACCTACATATATCATAAGTACACAATAAAGTCATTATTGACTC
ATG
>seq2
AGTTACCATGCCCAGCATTAACCCCCCTCAACAACCACCTCCGCCTATGAAGCCCGCCCG
GTAGGCGACATCAGCAAAGTGCCAACGCTGTATATATATATATTGATCACGAGCTACCAT
GCCAGCATAGCCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGAAACTCTGTTGC
ATGAACACACACACACACACATATTCACACATACTGCGTGGTGGCGACGCTCTCGACGTT
GAAATTAACATGAACCTACACACACACATACACACTGCGTAATCGCGACGTCCTCGTCTA
GACTCGCTATCTAGAACTAACGGGATCAAAAGCACTGCTGCTTGCCCGTGCGTATACATT
AAGAATAAAGCTTTCATCATTCTTGATCTTGACACCAAACCGAGCAGTTGATTTATTTAA
AGTGGCAAATATATATAACCTACATATATCATAAGTACACAATAAAGTCATTATTGACTC
>seq3
>seq4

>seq5
AGTTACCATG
)"""";

    std::ofstream fout(filename);
    fout<<text_in<<std::endl;
    fout.close();

    const char *text_exp = R""""(>1  seq1 {Cut} 1..100
AGTTACCATGCCCAGCATTAACCCCCCTCAACAACCACCTCCGCCTATGAAGCCCGCCCG
GTAGGCGACATCAGCAAAGTGCCAACGCTGTATATATATA
>2  seq1 {Cut} 91..190
TATATATATATATTGATCACGAGCTACCATGCCAGCATAGCCTCGTCCCCCGCTACCTGA
AACTCTGTTGCACCCGATGATGAAATCGGCATGAACACAC
>3  seq1 {Cut} 181..280
ATGAACACACACACACACACATATTCACACATACTGCGTGGTGGCGACGCTCTCGACGTT
GAAATTAACATGAACCTACACACACACATACACACTGCGT
>4  seq1 {Cut} 271..370
CACACTGCGTAATCGCGACGTCCTCGTCTAGACTCGCTATCTAGAACTAACGGGATCAAA
AGCACTGCTGCTTGCCCGTGCGTATACATTAAGAATAAAG
>5  seq1 {Cut} 361..460
AAGAATAAAGCTTTCATCATTCTTGATCTTGACACCAAACCGAGCAGTTGATTTATTTAA
AGTGGCAAATATATATAACCTACATATATCATAAGTACAC
>6  seq1 {Cut} 451..483
ATAAGTACACAATAAAGTCATTATTGACTCATG
>7  seq2 {Cut} 1..100
AGTTACCATGCCCAGCATTAACCCCCCTCAACAACCACCTCCGCCTATGAAGCCCGCCCG
GTAGGCGACATCAGCAAAGTGCCAACGCTGTATATATATA
>8  seq2 {Cut} 91..132
TATATATATATATTGATCACGAGCTACCATGCCAGCATAGCC
>9  seq2 {Cut} 167..266
TGAAACTCTGTTGCATGAACACACACACACACACATATTCACACATACTGCGTGGTGGCG
ACGCTCTCGACGTTGAAATTAACATGAACCTACACACACA
>10  seq2 {Cut} 257..356
TACACACACACATACACACTGCGTAATCGCGACGTCCTCGTCTAGACTCGCTATCTAGAA
CTAACGGGATCAAAAGCACTGCTGCTTGCCCGTGCGTATA
>11  seq2 {Cut} 347..446
CGTGCGTATACATTAAGAATAAAGCTTTCATCATTCTTGATCTTGACACCAAACCGAGCA
GTTGATTTATTTAAAGTGGCAAATATATATAACCTACATA
>12  seq2 {Cut} 437..480
AACCTACATATATCATAAGTACACAATAAAGTCATTATTGACTC)"""";
    fout.open(filename+"_cut_exp");
    fout<<text_exp<<std::endl;
    fout.close();

    Cutter cutter;            //!< parameter of bank cut-out
    int verbose=1;
    cutter.setLength(100);
    cutter.setOver(10);
    cutter.setWord(11);
    if(!cutter.check(filename,verbose))
    {
        SDGString bqName=cutter.cutDB(filename,verbose);
        if(verbose>0)
            std::cout<<"\n\nbank '"<<bqName<<"' created!!\n";
    }
    std::ifstream t_exp(filename+"_cut_exp");
    std::stringstream buffer_exp;
    buffer_exp << t_exp.rdbuf();

    std::ifstream t_obs(filename+"_cut");
    std::stringstream buffer_obs;
    buffer_obs << t_obs.rdbuf();

    CPPUNIT_ASSERT_EQUAL(buffer_exp.str(),buffer_obs.str());
}

//----------------------------------------------------------------------------------------------------------
void test_Tools::test_a_tool(std::string tool_name, std::string parameters)
{
    std::string obsFileName=tool_name+".outobs";
    std::string expFileName=tool_name+".outexp";
    std::string diff_result=tool_name+".diff";
    std::ostringstream cmd;
    cmd<<"../../../cmake-build-debug/src/tools/"<<tool_name;
    cmd<<" "<<parameters<<" > "<<tool_name<<".outobs";
    std::system(cmd.str().c_str());

    std::ostringstream cmd_diff;
    cmd_diff<<"diff --side-by-side --suppress-common-lines "<<obsFileName<<" "<<expFileName<<" > "<<diff_result;
    std::system(cmd_diff.str().c_str());

    std::ostringstream obsStr;
    std::ifstream fin_obs(diff_result);
    char buff[2048];
    while(fin_obs.getline(buff,2047,'\n'))
        obsStr<<buff<<std::endl;

    bool condition=(obsStr.str()=="");
    CPPUNIT_ASSERT_MESSAGE("Files "+obsFileName+" and "+expFileName+" are differents",condition);
    if(condition) {
        std::remove(diff_result.c_str());
        std::remove(obsFileName.c_str());
    }
}
void test_Tools::test_all_tools(void)
{
 /*
    TestTool hsearch "-o out ../DmelChr4.fa ../DmelChr4_refTEs.fa" out
    TestTool hrepeat "../DmelChr4_refTEs.fa" ../DmelChr4_refTEs.fa.hrepeat.set
*/

    test_a_tool("NWalign","seq1.fa seq2.fa");
    test_a_tool("SWalign","seq1.fa seq2.fa");
    test_a_tool("TRsearch","DmelChr4_refTEs.fa");
    test_a_tool("ltrsearch","DmelChr4_refTEs.fa");
    test_a_tool("fastlalign","seq1.fa seq2.fa");
    test_a_tool("rpt_map","seq12.fa 10 -8 4 2");

}