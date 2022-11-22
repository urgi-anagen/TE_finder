#include "Test_BLRMatchAlign.h"
#include "SDGString.h"
#include "FileUtils.h"
#include "BLRMatchAlign.h"
#include "Range.h"
#include "RangePair.h"
#include "RangePairSet.h"
#include "BLRMatcherThreadsParameter.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_BLRMatchAlign );
//---------------------------------------------------------------------------
void Test_BLRMatchAlign::setUp()
{
}
//---------------------------------------------------------------------------
void Test_BLRMatchAlign::tearDown()
{
}
//---------------------------------------------------------------------------
void Test_BLRMatchAlign::test_read()
{
   BLRMatcherThreadsParameter para = createParameter();
   BLRMatchAlign matchAlign;
 
   std::ostringstream inputData;
   inputData<<"chunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.836100\n";
   inputData<<"chunk682\t105091\t105206\trefTE_251\t3351\t3238\t0\t116\t76.106200\n";
   inputData<<"chunk682\t109527\t109818\trefTE_251\t3005\t3316\t0\t292\t82.846700\n";
   inputData<<"chunk682\t109601\t109708\trefTE_251\t3239\t3338\t3e-28\t108\t81.818200\n";
   inputData<<"chunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768000\n";
   inputData<<"chunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.860400\n";
   inputData<<"chunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.207600\n";
 
   std::istringstream inputDataStream(inputData.str());
   matchAlign.read(para,inputDataStream, 0);
   std::ostringstream obs;
   matchAlign.write(obs);

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
//---------------------------------------------------------------------------
void Test_BLRMatchAlign::test_load(void)
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
//---------------------------------------------------------------------------
void Test_BLRMatchAlign::test_set(void)
{

    BLRMatcherThreadsParameter para = createParameter();
    BLRMatchAlign match_align;

    std::list<RangePair> rp_list;
    rp_list.push_back(RangePair("chunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.836100\n"));
    rp_list.push_back(RangePair("chunk682\t105091\t105206\trefTE_251\t3351\t3238\t0\t116\t76.106200\n"));
    rp_list.push_back(RangePair("chunk682\t109527\t109818\trefTE_251\t3005\t3316\t0\t292\t82.846700\n"));
    rp_list.push_back(RangePair("chunk682\t109601\t109708\trefTE_251\t3239\t3338\t3e-28\t108\t81.818200\n"));
    rp_list.push_back(RangePair("chunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768000\n"));
    rp_list.push_back(RangePair("chunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.860400\n"));
    rp_list.push_back(RangePair("chunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.207600\n"));

    match_align.setFromRpsList(para, rp_list, 0);

    std::ostringstream sout_obs;
    match_align.write(sout_obs);

    std::ostringstream sout_exp;
    sout_exp<<"chunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.8604\n";
    sout_exp<<"chunk682\t109601\t109708\trefTE_251\t3239\t3338\t3e-28\t108\t81.8182\n";
    sout_exp<<"chunk682\t109527\t109818\trefTE_251\t3005\t3316\t0\t292\t82.8467\n";
    sout_exp<<"chunk682\t105091\t105206\trefTE_251\t3351\t3238\t0\t116\t76.1062\n";
    sout_exp<<"chunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.8361\n";
    sout_exp<<"chunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768\n";
    sout_exp<<"chunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.2076\n";

    CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
}























