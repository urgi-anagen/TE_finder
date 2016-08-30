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
	SDGString MATCHER_DATA_SUFFIX = "/TE_finder/matcher/";
	SDGString MATCHER_DATA = std::getenv("REPET_DATA") + MATCHER_DATA_SUFFIX;


	SDGString inputFileName = "input.align";
	SDGString cmd = "ln -s " + MATCHER_DATA + inputFileName;
	std::system(cmd);

	SDGString expFileName = "expInput.align.clean_match.path";
	cmd = "ln -s " + MATCHER_DATA + expFileName;
	std::system(cmd);

	SDGString obsFileName = "input.align.clean_match.path";
	SDGString outParamFileName = "input.align.clean_match.param";
	SDGString outMapFileName = "input.align.clean_match.map";

	cmd = "matcher -m " + inputFileName + " -j -x";
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

void Test_F_matcher::test_runAsScript_small_data( void ){

	SDGString inputFileName = "input.align";
	std::ofstream inputFileStream;
	FileUtils::openFile(inputFileName, inputFileStream);
	inputFileStream<<"chunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768000\n";
	inputFileStream<<"chunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.860400\n";
	inputFileStream<<"chunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.207600\n";
	inputFileStream.close();

	SDGString expFileName = "expInput.align.clean_match.path";
	std::ofstream expFileStream;
	FileUtils::openFile(expFileName, expFileStream);
	expFileStream<<"1\tchunk682\t109951\t109984\trefTE_266\t119\t149\t9.4e-19\t34\t77.768000\n";
	expFileStream<<"1\tchunk682\t110399\t110566\trefTE_266\t532\t686\t9.4e-19\t168\t77.768000\n";
	expFileStream<<"2\tchunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.860400\n";
	expFileStream<<"3\tchunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.207600\n";
	expFileStream.close();

	SDGString obsFileName = "input.align.clean_match.path";
	SDGString outParamFileName = "input.align.clean_match.param";
	SDGString outMapFileName = "input.align.clean_match.map";

	SDGString cmd = "matcher -m " + inputFileName + " -j -x";
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
