#include "Test_matcherThreads.h"
#include "SDGString.h"
#include "FileUtils.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_matcherThreads );

void Test_matcherThreads::setUp()
{
}

void Test_matcherThreads::tearDown()
{
}
void Test_matcherThreads::test_runAsProcess( void ){


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


	//set parameters
	BLRMatcherThreadsParameter para;
	para.setAlignFile(inputFileName);
	para.setPrefix_filename(inputFileName);
	para.setMerge(true);
	para.setJoin_frag(true);

	SDGString obsFileName = "input.align.match.path";
	SDGString obsFileNameBed = "input.align.match.bed";
	SDGString outParamFileName = "input.align.match.param";
//	SDGString outMapFileName = "input.align.clean_match.map";
//	SDGString obsFileNameAttr = "input.align.clean_match.path.attr";

    BLRMatcherThreads matcherThreads;
	matcherThreads.process(para);

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
//	FileUtils::removeFile(outMapFileName);
//	FileUtils::removeFile(obsFileNameAttr);


}
void Test_matcherThreads::test_runAsScript( void ){
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
//	SDGString outMapFileName = "input.align.clean_match.map";
//	SDGString obsFileNameAttr = "input.align.clean_match.path.attr";

	SDGString cmd = "../matcher.threads"+SDGString(VERSION)+" -m " + inputFileName + " -j -M -x -v 1";
	std::system(cmd);

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
//	FileUtils::removeFile(outMapFileName);
//	FileUtils::removeFile(obsFileNameAttr);


}
