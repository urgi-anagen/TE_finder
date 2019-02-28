#include "Test_BLRMatchPath.h"
#include "SDGString.h"
#include "FileUtils.h"
#include "../../matcher/BLRMatcherParameter.h"
#include "BLRMatchPath.h"
#include "Range.h"
#include "RangePair.h"
#include "RangePairSet.h"
#include "Test_BLRMatchMapUtils.h"
#include "Test_BLRMatchMapLoadUtils.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_BLRMatchPath );
//---------------------------------------------------------------------------
void Test_BLRMatchPath::setUp()
{
}
//---------------------------------------------------------------------------
void Test_BLRMatchPath::tearDown()
{
}
//---------------------------------------------------------------------------
void Test_BLRMatchPath::test_load()
{
   BLRMatcherThreadsParameter para = Test_BLRMatchMapUtils::createParameter();
   BLRMatchPath matchPath;
   SDGString path_file = "match.path";

	std::ostringstream inputData;
	inputData<<"1\tchunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.836100\n";
	inputData<<"1\tchunk682\t105091\t105206\trefTE_251\t3351\t3238\t0\t116\t76.106200\n";
	inputData<<"1\tchunk682\t109527\t109818\trefTE_251\t3005\t3316\t0\t292\t82.846700\n";
	inputData<<"1\tchunk682\t109601\t109708\trefTE_251\t3239\t3338\t3e-28\t108\t81.818200\n";
	inputData<<"2\tchunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768000\n";
	inputData<<"3\tchunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.860400\n";
	inputData<<"4\tchunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.207600\n";
	std::ofstream obs(path_file);
	obs<<inputData.str();
	obs.close();

	matchPath.load(para,path_file);

	std::ostringstream sout_obs;
	matchPath.write(sout_obs);
	matchPath.writeAttribute(sout_obs);

	std::ostringstream sout_exp;
	sout_exp<<"1\tchunk682\t109527\t109818\trefTE_251\t3005\t3316\t0\t292\t82.8467\n";
	sout_exp<<"1\tchunk682\t109601\t109708\trefTE_251\t3239\t3338\t3e-28\t108\t81.8182\n";
	sout_exp<<"1\tchunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.8361\n";
	sout_exp<<"1\tchunk682\t105091\t105206\trefTE_251\t3351\t3238\t0\t116\t76.1062\n";
	sout_exp<<"2\tchunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768\n";
    sout_exp<<"3\tchunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.8604\n";
	sout_exp<<"4\tchunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.2076\n";
	sout_exp<<"[1\tchunk682\t105085\t109818\trefTE_251\t3351\t2921\t0\t515\t82.4409]\n";
	sout_exp<<"[2\tchunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768]\n";
    sout_exp<<"[3\tchunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.8604]\n";
	sout_exp<<"[4\tchunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.2076]\n";


	CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
	FileUtils::removeFile(path_file);

	CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
}
//---------------------------------------------------------------------------
void Test_BLRMatchPath::test_read(void)
{
	BLRMatcherThreadsParameter para = Test_BLRMatchMapUtils::createParameter();
	BLRMatchPath matchPath;

	std::ostringstream inputData;
	inputData<<"1\tchunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.836100\n";
	inputData<<"1\tchunk682\t105091\t105206\trefTE_251\t3351\t3238\t0\t116\t76.106200\n";
	inputData<<"1\tchunk682\t109527\t109818\trefTE_251\t3005\t3316\t0\t292\t82.846700\n";
	inputData<<"1\tchunk682\t109601\t109708\trefTE_251\t3239\t3338\t3e-28\t108\t81.818200\n";
	inputData<<"2\tchunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768000\n";
	inputData<<"3\tchunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.860400\n";
	inputData<<"4\tchunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.207600\n";

	std::istringstream inputDataStream(inputData.str());
	matchPath.read(para,inputDataStream, 0);

	std::ostringstream sout_obs;
	matchPath.write(sout_obs);
	matchPath.writeAttribute(sout_obs);

	std::ostringstream sout_exp;
	sout_exp<<"1\tchunk682\t109527\t109818\trefTE_251\t3005\t3316\t0\t292\t82.8467\n";
	sout_exp<<"1\tchunk682\t109601\t109708\trefTE_251\t3239\t3338\t3e-28\t108\t81.8182\n";
	sout_exp<<"1\tchunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.8361\n";
	sout_exp<<"1\tchunk682\t105091\t105206\trefTE_251\t3351\t3238\t0\t116\t76.1062\n";
    sout_exp<<"2\tchunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768\n";
	sout_exp<<"3\tchunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.8604\n";
	sout_exp<<"4\tchunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.2076\n";
	sout_exp<<"[1\tchunk682\t105085\t109818\trefTE_251\t3351\t2921\t0\t515\t82.4409]\n";
 	sout_exp<<"[2\tchunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768]\n";
    sout_exp<<"[3\tchunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.8604]\n";
	sout_exp<<"[4\tchunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.2076]\n";


	CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
}
//---------------------------------------------------------------------------
void Test_BLRMatchPath::test_writeBED(void){

	BLRMatcherThreadsParameter para = Test_BLRMatchMapUtils::createParameter();
	BLRMatchPath matchPath;

	std::ostringstream inputData;
	inputData << "10\tCHR1v01212004\t100\t250\tTNAT1A\t36\t523\t2e-128\t460\t81.56\n";
	inputData << "10\tCHR1v01212004\t800\t1000\tTNAT1A\t36\t523\t2e-128\t460\t81.56\n";
	inputData << "20\tCHR1v01212004\t1050\t2000\tTNAT1A\t580\t480\t7e-78\t87\t79.14\n";

	std::istringstream inputDataStream(inputData.str());
	matchPath.read(para,inputDataStream, 0);



	std::ostringstream obs;
	SDGString color = "255,127,0";

	matchPath.writeBED(obs, color);

	std::ostringstream exp;
	exp << "CHR1v01212004\t100\t1000\tTNAT1A,TNAT1A\t868\t+\t100\t1000\t255,127,0\t2\t151,201,\t0,700,\n";
	exp << "CHR1v01212004\t1050\t2000\tTNAT1A\t87\t-\t1050\t2000\t255,127,0\n";

	CPPUNIT_ASSERT_EQUAL(exp.str(), obs.str());
}




















