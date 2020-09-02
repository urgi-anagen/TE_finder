#include <SDGMemBioSeq.h>

#include "Test_Hasher.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_Hasher);

//-------------------------------------------------------------------------
void Test_Hasher::test_runAsScript(void ){
    SDGString inputFileNameGenome = "DmelChr4.fa";
    SDGString inputFileNameTE = "DmelChr4_denovoLibTEs.fa";
    SDGString expFileName = "expDmelChr4.fa.hasher.align";

    SDGString prefixFileName = "test_runAsScript";
    SDGString obsFileName = "DmelChr4.fa.final.hasher.align";
    SDGString diff_result = prefixFileName+"result.txt";

    std::ostringstream cmd;
    cmd<<"../../../cmake-build-debug/src/hasher/hasher"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd<<" -w 15 -k 4 -d 5 -S 7 -n 1 -s 50 -c 100 -v 1 "<<inputFileNameGenome<<" "<<inputFileNameTE
            <<" &> "<<inputFileNameGenome<<".log";
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
        remove(diff_result.c_str());
        remove(obsFileName.c_str());

//        SDGString file=inputFileNameGenome+".1.duster.bed";
//        remove(file.c_str());
//        file=inputFileNameGenome+".1.duster.bed.fa";
//        remove(file.c_str());
//        file=inputFileNameTE+".kidx";
//        remove(file.c_str());
    }
}
//-------------------------------------------------------------------------
void Test_Hasher::test_fragMerge(void)
{
    unsigned word_len=10;
    unsigned word_dist=1;
    Hasher hsrch(word_len, word_dist);

    std::list< RangePair > frag;

    frag.push_back(RangePair(1, 10, 50,1,110,150,40,0.0, 0.0));
    frag.push_back(RangePair(1, 100, 150,1,110,150,40,0.0, 0.0));
    frag.push_back(RangePair(1, 110, 140,1,120,140,40,0.0, 0.0));

    hsrch.fragMerge(frag);

    frag.sort(RangePair::greater);
    std::ostringstream ostr_obs;
    for(std::list< RangePair >::iterator rp_it=frag.begin(); rp_it!=frag.end();rp_it++)
    {
        rp_it->writetxt(ostr_obs);
    }

    //std::cout<<"\n"<<ostr_obs.str()<<std::endl;
    std::list< RangePair > frag_exp;

    frag_exp.push_back(RangePair(1, 10, 50,1,110,150,40,0.0, 0.0));
    frag_exp.push_back(RangePair(1, 100, 150,1,110,150,40,0.0, 0.0));

    frag_exp.sort(RangePair::greater);
    std::ostringstream ostr_exp;
    for(std::list< RangePair >::iterator rp_it=frag_exp.begin(); rp_it!=frag_exp.end();rp_it++)
    {
        rp_it->writetxt(ostr_exp);
    }

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}