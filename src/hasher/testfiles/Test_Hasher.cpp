#include <SDGMemBioSeq.h>

#include "Test_Hasher.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_Hasher);


void Test_Hasher::test_runAsScript(void ){
    SDGString inputFileNameGenome = "DmelChr4.fa";
    SDGString inputFileNameTE = "DmelChr4_denovoLibTEs.fa";
    SDGString expFileName = "expDmelChr4.fa.hasher.align";

    SDGString prefixFileName = "test_runAsScript";
    SDGString obsFileName = "DmelChr4.fa.hasher.align";
    SDGString diff_result = prefixFileName+"result.txt";

    std::ostringstream cmd;
    cmd<<"../../../cmake-build-debug/src/hasher/hasher"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd<<" -w 15 -k 4 -d 5 -S 7 -s 50 -c 100 -v 1 "<<inputFileNameGenome<<" "<<inputFileNameTE;
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