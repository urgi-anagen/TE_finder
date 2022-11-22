#include "Test_BLRMatchMapLoad.h"
#include "SDGString.h"
#include "FileUtils.h"
#include "BLRMatcherThreadsParameter.h"
#include "BLRMatchMap.h"
#include "Range.h"
#include "RangePair.h"
#include "RangePairSet.h"
#include "BLRRangePairSetList.h"


CPPUNIT_TEST_SUITE_REGISTRATION( Test_BLRMatchMapLoad );

void Test_BLRMatchMapLoad::setUp()
{
}

void Test_BLRMatchMapLoad::tearDown()
{
}

void Test_BLRMatchMapLoad::test_readAlign()
{
   BLRMatcherThreadsParameter para = createParameter();
   BLRMatchMap matchMap(para);
 
   std::ostringstream inputData;
   inputData<<"chunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.836100\n";
   inputData<<"chunk682\t105091\t105206\trefTE_251\t3351\t3238\t0\t116\t76.106200\n";
   inputData<<"chunk682\t109527\t109818\trefTE_251\t3005\t3316\t0\t292\t82.846700\n";
   inputData<<"chunk682\t109601\t109708\trefTE_251\t3239\t3338\t3e-28\t108\t81.818200\n";
   inputData<<"chunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768000\n";
   inputData<<"chunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.860400\n";
   inputData<<"chunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.207600\n";
 
   std::istringstream inputDataStream(inputData.str());
   matchMap.readAlign(inputDataStream, 0);
   std::ostringstream obs;
   matchMap.writeMapAlign(obs);

   std::ostringstream exp;
   exp<<"chunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.8604\n";
   exp<<"chunk682\t109601\t109708\trefTE_251\t3239\t3338\t3e-28\t108\t81.8182\n";
   exp<<"chunk682\t109527\t109818\trefTE_251\t3005\t3316\t0\t292\t82.8467\n";
   exp<<"chunk682\t105091\t105206\trefTE_251\t3351\t3238\t0\t116\t76.1062\n";
   exp<<"chunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.8361\n";
   exp<<"chunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768\n";
   exp<<"chunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.2076\n";

   CPPUNIT_ASSERT_EQUAL(exp.str(), obs.str());
}

void Test_BLRMatchMapLoad::test_readPath() {
    std::ostringstream inputData;
    inputData << "1\tchunk1\t105091\t105206\tchunk1\t3351\t3238\t0\t116\t76.1062\n";
    inputData << "2\tchunk1\t109601\t109708\tchunk1\t3239\t3338\t3e-28\t108\t81.8182\n";
    inputData << "3\tchunk1\t105085\t105353\tchunk1\t3202\t2921\t0\t269\t84.8361\n";
    inputData << "3\tchunk1\t109985\t110398\tchunk1\t2429\t2033\t0\t414\t77.8604\n";

    BLRMatcherThreadsParameter para = createParameter();
    BLRMatchMap matchMap(para);

    std::istringstream inputDataStream(inputData.str());
    matchMap.readPath(inputDataStream, 0);

    BLRRangePairSetList rpsList = matchMap.getRpsList();
    std::ostringstream sout_obs;
    rpsList.writePath(sout_obs);
    rpsList.writePathAttr(sout_obs);

    std::ostringstream sout_exp;
    sout_exp << "1\tchunk1\t105091\t105206\tchunk1\t3351\t3238\t0\t116\t76.1062\n";
    sout_exp << "2\tchunk1\t109601\t109708\tchunk1\t3239\t3338\t3e-28\t108\t81.8182\n";
    sout_exp << "3\tchunk1\t105085\t105353\tchunk1\t3202\t2921\t0\t269\t84.8361\n";
    sout_exp << "3\tchunk1\t109985\t110398\tchunk1\t2429\t2033\t0\t414\t77.8604\n";
    sout_exp << "[1\tchunk1\t105091\t105206\tchunk1\t3351\t3238\t0\t116\t76.1062]\n";
    sout_exp << "[2\tchunk1\t109601\t109708\tchunk1\t3239\t3338\t3e-28\t108\t81.8182]\n";
    sout_exp << "[3\tchunk1\t105085\t110398\tchunk1\t3202\t2033\t0\t427\t80.6868]\n";


    CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
}

void Test_BLRMatchMapLoad::test_load(void)
{
 	SDGString match_file = "match.align";


  	BLRMatcherThreadsParameter para = createParameter();
   	BLRMatchMap matchMap(para);
 
	std::ostringstream inputData;
	inputData<<"chunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.836100\n";
	inputData<<"chunk682\t105091\t105206\trefTE_251\t3351\t3238\t0\t116\t76.106200\n";
	inputData<<"chunk682\t109527\t109818\trefTE_251\t3005\t3316\t0\t292\t82.846700\n";
	inputData<<"chunk682\t109601\t109708\trefTE_251\t3239\t3338\t3e-28\t108\t81.818200\n";
	inputData<<"chunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768000\n";
	inputData<<"chunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.860400\n";
	inputData<<"chunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.207600\n";

	std::ofstream obs(match_file);
	obs<<inputData.str();
	obs.close();

  	matchMap.loadAlign(match_file,0);

   	std::ostringstream sout_obs;
   	matchMap.writeMapAlign(sout_obs);

	std::ostringstream sout_exp;
	sout_exp<<"chunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.8604\n";
	sout_exp<<"chunk682\t109601\t109708\trefTE_251\t3239\t3338\t3e-28\t108\t81.8182\n";
	sout_exp<<"chunk682\t109527\t109818\trefTE_251\t3005\t3316\t0\t292\t82.8467\n";
	sout_exp<<"chunk682\t105091\t105206\trefTE_251\t3351\t3238\t0\t116\t76.1062\n";
	sout_exp<<"chunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.8361\n";
	sout_exp<<"chunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768\n";
	sout_exp<<"chunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.2076\n";

	CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
	FileUtils::removeFile(match_file);
}

void Test_BLRMatchMapLoad::test_loadPath(void)
{
	SDGString path_file = "match.path";

	std::ostringstream inputData;
    inputData<<"1\tchunk1\t105091\t105206\tchunk1\t3351\t3238\t0\t116\t76.1062\n";
    inputData<<"2\tchunk1\t109601\t109708\tchunk1\t3239\t3338\t3e-28\t108\t81.8182\n";
    inputData<<"3\tchunk1\t105085\t105353\tchunk1\t3202\t2921\t0\t269\t84.8361\n";
    inputData<<"3\tchunk1\t109985\t110398\tchunk1\t2429\t2033\t0\t414\t77.8604\n";
	std::ofstream obs(path_file);
	obs<<inputData.str();
	obs.close();

	BLRMatcherThreadsParameter para = createParameter();
	BLRMatchMap matchMap(para);

	matchMap.loadPath(path_file);

	// sort obs list
	std::list<RangePairSet> rpsList = matchMap.getRawRpsList();
    matchMap.setRawRpsList(rpsList);

	std::ostringstream sout_obs;
	matchMap.writeRpsList(rpsList,sout_obs);
	matchMap.writeRpsListAttribute(rpsList,sout_obs);

	std::ostringstream sout_exp;
    sout_exp<<"1\tchunk1\t105091\t105206\tchunk1\t3351\t3238\t0\t116\t76.1062\n";
    sout_exp<<"2\tchunk1\t109601\t109708\tchunk1\t3239\t3338\t3e-28\t108\t81.8182\n";
    sout_exp<<"3\tchunk1\t105085\t105353\tchunk1\t3202\t2921\t0\t269\t84.8361\n";
    sout_exp<<"3\tchunk1\t109985\t110398\tchunk1\t2429\t2033\t0\t414\t77.8604\n";
    sout_exp<<"[1\tchunk1\t105091\t105206\tchunk1\t3351\t3238\t0\t116\t76.1062]\n";
    sout_exp<<"[2\tchunk1\t109601\t109708\tchunk1\t3239\t3338\t3e-28\t108\t81.8182]\n";
    sout_exp<<"[3\tchunk1\t105085\t110398\tchunk1\t3202\t2033\t0\t427\t80.6868]\n";


	CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
	FileUtils::removeFile(path_file);
}
























