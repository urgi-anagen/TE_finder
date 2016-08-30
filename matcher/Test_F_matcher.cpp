#include "Test_F_matcher.h"
#include "SDGString.h"
#include "FileUtils.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_F_matcher );

void Test_F_matcher::setUp()
{
}

void Test_F_matcher::tearDown()
{
}

void Test_F_matcher::test_runAsScript_bigData( void ){
	SDGString MATCHER_DATA_SUFFIX = "/TE_finder/matcher/";
	SDGString MATCHER_DATA = std::getenv("REPET_DATA") + MATCHER_DATA_SUFFIX;

	SDGString inputFileName = "input.align";
	SDGString inputFileNameWithTestName = "input_test_runAsScript_bigData.align";
	SDGString cmd = "ln -s " + MATCHER_DATA + inputFileName + " " + inputFileNameWithTestName;
	std::system(cmd);
	
  	SDGString expFileName = "expInput.align.clean_match.path";
	SDGString expFileNameWithTestName = "expInput_test_runAsScript_bigData.align.clean_match.path";
	cmd = "ln -s " + MATCHER_DATA + expFileName + " " + expFileNameWithTestName;
	std::system(cmd);

	SDGString obsFileName = "input_test_runAsScript_bigData.align.clean_match.path";
	SDGString outParamFileName = "input_test_runAsScript_bigData.align.clean_match.param";
	SDGString outMapFileName = "input_test_runAsScript_bigData.align.clean_match.map";

	cmd = "./matcher2.27 -m " + inputFileNameWithTestName + " -j -x";
	std::cout<<" "<<std::endl;
	std::cout<<"cmd: "<<cmd<<std::endl;
  	std::system(cmd);

	bool obs = FileUtils::areTwoFilesIdentical(expFileNameWithTestName, obsFileName);
	bool exp = true;

	CPPUNIT_ASSERT_EQUAL(exp, obs);
	/*
	FileUtils::removeFile(inputFileNameWithTestName);
	FileUtils::removeFile(expFileNameWithTestName);
	FileUtils::removeFile(obsFileName);
	FileUtils::removeFile(outParamFileName);
	FileUtils::removeFile(outMapFileName);
	*/
}


