
#include "SDGString.h"
#include "FileUtils.h"
#include "Test_matcherThreads.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_matcherThreads );

void Test_matcherThreads::setUp()
{
}

void Test_matcherThreads::tearDown()
{
}
void Test_matcherThreads::test_runAsScript_join_simple( void ){
	SDGString inputFileName = "input.align";
	std::ofstream inputFileStream(inputFileName);
	inputFileStream<<"chunk1\t100\t500\trefTE1\t100\t500\t9.4e-19\t400\t100.00\n";
	inputFileStream<<"chunk1\t600\t1000\trefTE1\t600\t1000\t0\t400\t100.00\n";
	inputFileStream<<"chunk1\t1000\t2000\trefTE1\t1500\t2500\t0\t1000\t100.00\n";
	inputFileStream.close();

	std::ostringstream expStr;
	expStr<<"1\tchunk1\t100\t500\trefTE1\t100\t500\t9.4e-19\t401\t100\n";
	expStr<<"1\tchunk1\t600\t1000\trefTE1\t600\t1000\t0\t401\t100\n";
	expStr<<"1\tchunk1\t1000\t2000\trefTE1\t1500\t2500\t0\t1001\t100\n";

	SDGString obsFileName = "input.align.clean_match.path";
	SDGString obsFileNameBed = "input.align.clean_match.bed";
	SDGString outParamFileName = "input.align.clean_match.param";
    SDGString obsGFF3FileName = "input.align.clean_match.gff3";

    std::ostringstream cmd;
	cmd<<"../matcherThreads"<<std::fixed<<std::setprecision(2)<<VERSION;
	cmd<<" -m "<<inputFileName<<" -j -M -x -v 1";
	std::system(cmd.str().c_str());

	std::ostringstream obsStr;
	std::ifstream fin(obsFileName);
	char buff[1024];
	while(fin.getline(buff,1023,'\n'))
		obsStr<<buff<<std::endl;


	CPPUNIT_ASSERT_EQUAL(expStr.str(), obsStr.str());

	FileUtils::removeFile(inputFileName);
	FileUtils::removeFile(obsFileNameBed);
	FileUtils::removeFile(obsFileName);
	FileUtils::removeFile(outParamFileName);
    FileUtils::removeFile(obsGFF3FileName);
}
void Test_matcherThreads::test_runAsScript_join( void ){
    SDGString inputFileName = "blasterAthaBest.align";
    SDGString expFileName = "blasterAthaBestNoMerge.match.path";

    SDGString prefixFileName = "test_runAsScript_join_threads_";
    SDGString prefixObsFileName = "obs.match";
    SDGString obsFileName = prefixFileName+prefixObsFileName+".path";
    SDGString diff_result = prefixFileName+"result.txt";
    SDGString obsGFF3FileName = prefixFileName+prefixObsFileName+".gff3";

    std::ostringstream cmd;
    cmd<<"../matcherThreads"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd<<" -m "<<inputFileName<<" -j -B "<<prefixFileName<<"obs";
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
        FileUtils::removeFile(diff_result);
        FileUtils::removeFile(obsFileName);
        FileUtils::removeFile(obsGFF3FileName);

        SDGString file=prefixFileName+prefixObsFileName+".param";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixObsFileName+".bed";
        FileUtils::removeFile(file);
    }
}
void Test_matcherThreads::test_runAsScript_join_threads( void ){
    SDGString inputFileName = "blasterAthaBest.align";
    SDGString prefixFileName = "test_runAsScript_join_threads_";
    SDGString prefixObsFileName = "obs.match";
    SDGString prefixExpFileName = "exp.match";
    SDGString obsFileName = prefixFileName+prefixObsFileName+".path";
    SDGString diff_result = prefixFileName+"result.txt";
    SDGString expFileName = prefixFileName+prefixExpFileName+".path";
    SDGString obsGFF3FileName = prefixFileName+prefixObsFileName+".gff3";
    SDGString expGFF3FileName = prefixFileName+prefixExpFileName+".gff3";

    std::ostringstream cmd_exp;
    cmd_exp<<"../matcherThreads"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd_exp<<" -m "<<inputFileName<<" -j -B "<<prefixFileName<<"exp";
    std::system(cmd_exp.str().c_str());

    std::ostringstream cmd_obs;
    cmd_obs<<"../matcherThreads"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd_obs<<" -m "<<inputFileName<<" -t 2 -j -B "<<prefixFileName<<"obs";
    std::system(cmd_obs.str().c_str());

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
    if(condition)
    {
        FileUtils::removeFile(diff_result);
        FileUtils::removeFile(expFileName);
        FileUtils::removeFile(obsFileName);
        FileUtils::removeFile(obsGFF3FileName);
        FileUtils::removeFile(expGFF3FileName);

        SDGString file=prefixFileName+prefixObsFileName+".param";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixObsFileName+".bed";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixExpFileName+".param";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixExpFileName+".bed";
        FileUtils::removeFile(file);
    }
}
void Test_matcherThreads::test_runAsScript_join_clean_threads( void ){
    SDGString inputFileName = "blasterAthaBest.align";
    SDGString prefixFileName = "test_runAsScript_join_clean_threads_";
    SDGString prefixObsFileName = "obs.clean_match";
    SDGString prefixExpFileName = "exp.clean_match";
    SDGString obsFileName = prefixFileName+prefixObsFileName+".path";
    SDGString diff_result = prefixFileName+"result.txt";
    SDGString expFileName = prefixFileName+prefixExpFileName+".path";
    SDGString obsGFF3FileName = prefixFileName+prefixObsFileName+".gff3";
    SDGString expGFF3FileName = prefixFileName+prefixExpFileName+".gff3";

    std::ostringstream cmd_exp;
    cmd_exp<<"../matcherThreads"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd_exp<<" -m "<<inputFileName<<" -j -x -B "<<prefixFileName<<"exp";
    std::system(cmd_exp.str().c_str());

    std::ostringstream cmd_obs;
    cmd_obs<<"../matcherThreads"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd_obs<<" -m "<<inputFileName<<" -t 2 -j -x -B "<<prefixFileName<<"obs";
    std::system(cmd_obs.str().c_str());

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
    if(condition)
    {
        FileUtils::removeFile(diff_result);
        FileUtils::removeFile(expFileName);
        FileUtils::removeFile(obsFileName);
        FileUtils::removeFile(obsGFF3FileName);
        FileUtils::removeFile(expGFF3FileName);

        SDGString file=prefixFileName+prefixObsFileName+".param";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixObsFileName+".bed";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixExpFileName+".param";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixExpFileName+".bed";
        FileUtils::removeFile(file);
    }

}
void Test_matcherThreads::test_runAsScript_join_merge_threads( void ){

    SDGString inputFileName = "blasterAthaBest.align";
    SDGString prefixFileName = "test_runAsScript_join_merge_threads_";
    SDGString prefixObsFileName = "obs.match";
    SDGString prefixExpFileName = "exp.match";
    SDGString obsFileName = prefixFileName+prefixObsFileName+".path";
    SDGString diff_result = prefixFileName+"result.txt";
    SDGString expFileName = prefixFileName+prefixExpFileName+".path";
    SDGString obsGFF3FileName = prefixFileName+prefixObsFileName+".gff3";
    SDGString expGFF3FileName = prefixFileName+prefixExpFileName+".gff3";

    std::ostringstream cmd_exp;
    cmd_exp<<"../matcherThreads"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd_exp<<" -m "<<inputFileName<<" -j -M -B "<<prefixFileName<<"exp";
    std::system(cmd_exp.str().c_str());

    std::ostringstream cmd_obs;
    cmd_obs<<"../matcherThreads"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd_obs<<" -m "<<inputFileName<<" -t 2 -j -M -B "<<prefixFileName<<"obs";
    std::system(cmd_obs.str().c_str());

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
    if(condition)
    {
        FileUtils::removeFile(diff_result);
        FileUtils::removeFile(expFileName);
        FileUtils::removeFile(obsFileName);
        FileUtils::removeFile(obsGFF3FileName);
        FileUtils::removeFile(expGFF3FileName);

        SDGString file=prefixFileName+prefixObsFileName+".param";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixObsFileName+".bed";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixExpFileName+".param";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixExpFileName+".bed";
        FileUtils::removeFile(file);
    }

}
void Test_matcherThreads::test_runAsScript_join_merge_clean_threads( void ){

    SDGString inputFileName = "blasterAthaBest.align";
    SDGString prefixFileName = "test_runAsScript_join_merge_clean_threads_";
    SDGString prefixObsFileName = "obs.clean_match";
    SDGString prefixExpFileName = "exp.clean_match";
    SDGString obsFileName = prefixFileName+prefixObsFileName+".path";
    SDGString diff_result = prefixFileName+"result.txt";
    SDGString expFileName = prefixFileName+prefixExpFileName+".path";
    SDGString obsGFF3FileName = prefixFileName+prefixObsFileName+".gff3";
    SDGString expGFF3FileName = prefixFileName+prefixExpFileName+".gff3";


    std::ostringstream cmd_exp;
    cmd_exp<<"../matcherThreads"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd_exp<<" -m "<<inputFileName<<" -j -M -x -B "<<prefixFileName<<"exp";
    std::system(cmd_exp.str().c_str());



    std::ostringstream cmd_obs;
    cmd_obs<<"../matcherThreads"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd_obs<<" -m "<<inputFileName<<" -t 2 -j -x -M -B "<<prefixFileName<<"obs";
    std::system(cmd_obs.str().c_str());

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
    if(condition)
    {
        FileUtils::removeFile(diff_result);
        FileUtils::removeFile(expFileName);
        FileUtils::removeFile(obsFileName);
        FileUtils::removeFile(obsGFF3FileName);
        FileUtils::removeFile(expGFF3FileName);

        SDGString file=prefixFileName+prefixObsFileName+".param";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixObsFileName+".bed";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixExpFileName+".param";
        FileUtils::removeFile(file);

        file=prefixFileName+prefixExpFileName+".bed";
        FileUtils::removeFile(file);
    }

}