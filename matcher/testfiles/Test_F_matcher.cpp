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

void Test_F_matcher::test_runAsScript( void ){
	SDGString inputFileName = "input.align";
	std::ofstream inputFileStream;
	FileUtils::openFile(inputFileName, inputFileStream);
	inputFileStream<<"chunk1\t100\t500\trefTE1\t100\t500\t9.4e-19\t400\t100.00\n";
	inputFileStream<<"chunk1\t600\t1000\trefTE1\t600\t1000\t0\t400\t100.00\n";
	inputFileStream<<"chunk1\t1000\t2000\trefTE1\t1500\t2500\t0\t1000\t100.00\n";
	inputFileStream.close();

	SDGString expFileName = "expInput.align.clean_match.path";
	std::ofstream expFileStream;
	FileUtils::openFile(expFileName, expFileStream);
	expFileStream<<"1\tchunk1\t100\t500\trefTE1\t100\t500\t9.4e-19\t401\t100\n";
    expFileStream<<"1\tchunk1\t600\t1000\trefTE1\t600\t1000\t0\t401\t100\n";
    expFileStream<<"1\tchunk1\t1000\t2000\trefTE1\t1500\t2500\t0\t1001\t100\n";
	expFileStream.close();

	SDGString obsFileName = "input.align.clean_match.path";
	SDGString outParamFileName = "input.align.clean_match.param";
	SDGString outMapFileName = "input.align.clean_match.map";

	SDGString cmd = "../matcher"+SDGString(VERSION)+" -m " + inputFileName + " -j -x";
	std::system(cmd);

	bool obs = FileUtils::areTwoFilesIdentical(expFileName, obsFileName);
	bool exp = true;

	CPPUNIT_ASSERT_EQUAL(exp, obs);
	FileUtils::removeFile(inputFileName);
	FileUtils::removeFile(expFileName);
	FileUtils::removeFile(obsFileName);
	FileUtils::removeFile(outParamFileName);
	FileUtils::removeFile(outMapFileName);
}
