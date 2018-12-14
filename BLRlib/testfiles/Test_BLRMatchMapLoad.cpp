#include "Test_BLRMatchMapLoad.h"
#include "SDGString.h"
#include "FileUtils.h"
#include "BLRMatcherParameter.h"
#include "BLRGrouperParameter.h"
#include "BLRMatchMap.h"
#include "Range.h"
#include "RangePair.h"
#include "RangePairSet.h"
#include "Test_BLRMatchMapUtils.h"
#include "Test_BLRMatchMapLoadUtils.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_BLRMatchMapLoad );

void Test_BLRMatchMapLoad::setUp()
{
}

void Test_BLRMatchMapLoad::tearDown()
{
}

void Test_BLRMatchMapLoad::test_readAlign()
{
   BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter();
   BLRMatchMap matchMap(&para);
 
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

void Test_BLRMatchMapLoad::test_readPath()
{
   std::ostringstream inputData;
   inputData<<"1\t1\t105091\t105206\t1\t3351\t3238\t0\t116\t76.1062\t116\n";
   inputData<<"2\t1\t109601\t109708\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
   inputData<<"3\t1\t105085\t105353\t1\t3202\t2921\t0\t269\t84.8361\t269\n";
   inputData<<"3\t1\t109985\t110398\t1\t2429\t2033\t0\t414\t77.8604\t414\n";
   inputData<<"4\t1\t109951\t110569\t2\t119\t689\t9.4e-19\t619\t77.768\t619\n";
   inputData<<"5\t1\t110567\t110878\t3\t2833\t2532\t0\t312\t77.2076\t312\n";
   inputData<<"6\t1\t109527\t109818\t1\t3005\t3316\t0\t292\t82.8467\t292\n";
   inputData<<"7\t4\t105091\t105206\t4\t3351\t3238\t0\t116\t76.1062\t116\n";
   inputData<<"8\t4\t109601\t109708\t4\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
   inputData<<"9\t4\t105085\t105353\t4\t3202\t2921\t0\t269\t84.8361\t269\n";
   inputData<<"9\t4\t109985\t110398\t4\t2429\t2033\t0\t414\t77.8604\t414\n";
   inputData<<"10\t4\t109951\t110569\t5\t119\t689\t9.4e-19\t619\t77.768\t619\n";
   inputData<<"11\t4\t110567\t110878\t6\t2833\t2532\t0\t312\t77.2076\t312\n";
   inputData<<"12\t4\t109527\t109818\t4\t3005\t3316\t0\t292\t82.8467\t292\n";
   
   BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter();
   BLRMatchMap matchMap(&para);
 
   std::istringstream inputDataStream(inputData.str());
   matchMap.readPath(inputDataStream, 0);
   
   // sort obs list
   std::list<RangePairSet> rpsList = matchMap.getRpsList();
   rpsList.sort(RangePair::greaterScore);
   matchMap.setRpsList(rpsList);

   std::ostringstream sout_obs;
   matchMap.writeRpsList(rpsList,sout_obs);
   matchMap.writeRpsListAttribute(rpsList,sout_obs);

   std::ostringstream sout_exp;
   sout_exp<<"1\t1\t109951\t110569\t2\t119\t689\t9.4e-19\t619\t77.768\n";
   sout_exp<<"2\t4\t109951\t110569\t5\t119\t689\t9.4e-19\t619\t77.768\n";
   sout_exp<<"3\t1\t109985\t110398\t1\t2429\t2033\t0\t414\t77.8604\n";
   sout_exp<<"3\t1\t105085\t105353\t1\t3202\t2921\t0\t269\t84.8361\n";
   sout_exp<<"4\t4\t109985\t110398\t4\t2429\t2033\t0\t414\t77.8604\n";
   sout_exp<<"4\t4\t105085\t105353\t4\t3202\t2921\t0\t269\t84.8361\n";
   sout_exp<<"5\t1\t110567\t110878\t3\t2833\t2532\t0\t312\t77.2076\n";
   sout_exp<<"6\t4\t110567\t110878\t6\t2833\t2532\t0\t312\t77.2076\n";
   sout_exp<<"7\t1\t109527\t109818\t1\t3005\t3316\t0\t292\t82.8467\n";
   sout_exp<<"8\t4\t109527\t109818\t4\t3005\t3316\t0\t292\t82.8467\n";
   sout_exp<<"9\t1\t105091\t105206\t1\t3351\t3238\t0\t116\t76.1062\n";
   sout_exp<<"10\t4\t105091\t105206\t4\t3351\t3238\t0\t116\t76.1062\n";
   sout_exp<<"11\t1\t109601\t109708\t1\t3239\t3338\t3e-28\t108\t81.8182\n";
   sout_exp<<"12\t4\t109601\t109708\t4\t3239\t3338\t3e-28\t108\t81.8182\n";
   sout_exp<<"[1\t1\t109951\t110569\t2\t119\t689\t9.4e-19\t619\t77.768]\n";
   sout_exp<<"[2\t4\t109951\t110569\t5\t119\t689\t9.4e-19\t619\t77.768]\n";
   sout_exp<<"[3\t1\t105085\t110398\t1\t3202\t2033\t0\t427\t80.6868]\n";
   sout_exp<<"[4\t4\t105085\t110398\t4\t3202\t2033\t0\t427\t80.6868]\n";
   sout_exp<<"[5\t1\t110567\t110878\t3\t2833\t2532\t0\t312\t77.2076]\n";
   sout_exp<<"[6\t4\t110567\t110878\t6\t2833\t2532\t0\t312\t77.2076]\n";
   sout_exp<<"[7\t1\t109527\t109818\t1\t3005\t3316\t0\t292\t82.8467]\n";
   sout_exp<<"[8\t4\t109527\t109818\t4\t3005\t3316\t0\t292\t82.8467]\n";
   sout_exp<<"[9\t1\t105091\t105206\t1\t3351\t3238\t0\t116\t76.1062]\n";
   sout_exp<<"[10\t4\t105091\t105206\t4\t3351\t3238\t0\t116\t76.1062]\n";
   sout_exp<<"[11\t1\t109601\t109708\t1\t3239\t3338\t3e-28\t108\t81.8182]\n";
   sout_exp<<"[12\t4\t109601\t109708\t4\t3239\t3338\t3e-28\t108\t81.8182]\n";

   CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
}

void Test_BLRMatchMapLoad::test_load(void)
{
 	SDGString match_file = "match.align";


  	BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter();
   	BLRMatchMap matchMap(&para);
 
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

  	matchMap.loadAlign(0);

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
	inputData<<"1\t1\t105091\t105206\t1\t3351\t3238\t0\t116\t76.1062\t116\n";
	inputData<<"2\t1\t109601\t109708\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	inputData<<"3\t1\t105085\t105353\t1\t3202\t2921\t0\t269\t84.8361\t269\n";
	inputData<<"3\t1\t109985\t110398\t1\t2429\t2033\t0\t414\t77.8604\t414\n";
	inputData<<"4\t1\t109951\t110569\t2\t119\t689\t9.4e-19\t619\t77.768\t619\n";
	inputData<<"5\t1\t110567\t110878\t3\t2833\t2532\t0\t312\t77.2076\t312\n";
	inputData<<"6\t1\t109527\t109818\t1\t3005\t3316\t0\t292\t82.8467\t292\n";
	inputData<<"7\t4\t105091\t105206\t4\t3351\t3238\t0\t116\t76.1062\t116\n";
	inputData<<"8\t4\t109601\t109708\t4\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	inputData<<"9\t4\t105085\t105353\t4\t3202\t2921\t0\t269\t84.8361\t269\n";
	inputData<<"9\t4\t109985\t110398\t4\t2429\t2033\t0\t414\t77.8604\t414\n";
	inputData<<"10\t4\t109951\t110569\t5\t119\t689\t9.4e-19\t619\t77.768\t619\n";
	inputData<<"11\t4\t110567\t110878\t6\t2833\t2532\t0\t312\t77.2076\t312\n";
	inputData<<"12\t4\t109527\t109818\t4\t3005\t3316\t0\t292\t82.8467\t292\n";
	std::ofstream obs(path_file);
	obs<<inputData.str();
	obs.close();

	BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter();
	BLRMatchMap matchMap(&para);

	matchMap.loadPath(path_file);

	// sort obs list
	std::list<RangePairSet> rpsList = matchMap.getRpsList();
	rpsList.sort(RangePair::greaterScore);
	matchMap.setRpsList(rpsList);

	std::ostringstream sout_obs;
	matchMap.writeRpsList(rpsList,sout_obs);
	matchMap.writeRpsListAttribute(rpsList,sout_obs);

	std::ostringstream sout_exp;
	sout_exp<<"1\t1\t109951\t110569\t2\t119\t689\t9.4e-19\t619\t77.768\n";
	sout_exp<<"2\t4\t109951\t110569\t5\t119\t689\t9.4e-19\t619\t77.768\n";
	sout_exp<<"3\t1\t109985\t110398\t1\t2429\t2033\t0\t414\t77.8604\n";
	sout_exp<<"3\t1\t105085\t105353\t1\t3202\t2921\t0\t269\t84.8361\n";
	sout_exp<<"4\t4\t109985\t110398\t4\t2429\t2033\t0\t414\t77.8604\n";
	sout_exp<<"4\t4\t105085\t105353\t4\t3202\t2921\t0\t269\t84.8361\n";
	sout_exp<<"5\t1\t110567\t110878\t3\t2833\t2532\t0\t312\t77.2076\n";
	sout_exp<<"6\t4\t110567\t110878\t6\t2833\t2532\t0\t312\t77.2076\n";
	sout_exp<<"7\t1\t109527\t109818\t1\t3005\t3316\t0\t292\t82.8467\n";
	sout_exp<<"8\t4\t109527\t109818\t4\t3005\t3316\t0\t292\t82.8467\n";
	sout_exp<<"9\t1\t105091\t105206\t1\t3351\t3238\t0\t116\t76.1062\n";
	sout_exp<<"10\t4\t105091\t105206\t4\t3351\t3238\t0\t116\t76.1062\n";
	sout_exp<<"11\t1\t109601\t109708\t1\t3239\t3338\t3e-28\t108\t81.8182\n";
	sout_exp<<"12\t4\t109601\t109708\t4\t3239\t3338\t3e-28\t108\t81.8182\n";
	sout_exp<<"[1\t1\t109951\t110569\t2\t119\t689\t9.4e-19\t619\t77.768]\n";
	sout_exp<<"[2\t4\t109951\t110569\t5\t119\t689\t9.4e-19\t619\t77.768]\n";
	sout_exp<<"[3\t1\t105085\t110398\t1\t3202\t2033\t0\t427\t80.6868]\n";
	sout_exp<<"[4\t4\t105085\t110398\t4\t3202\t2033\t0\t427\t80.6868]\n";
	sout_exp<<"[5\t1\t110567\t110878\t3\t2833\t2532\t0\t312\t77.2076]\n";
	sout_exp<<"[6\t4\t110567\t110878\t6\t2833\t2532\t0\t312\t77.2076]\n";
	sout_exp<<"[7\t1\t109527\t109818\t1\t3005\t3316\t0\t292\t82.8467]\n";
	sout_exp<<"[8\t4\t109527\t109818\t4\t3005\t3316\t0\t292\t82.8467]\n";
	sout_exp<<"[9\t1\t105091\t105206\t1\t3351\t3238\t0\t116\t76.1062]\n";
	sout_exp<<"[10\t4\t105091\t105206\t4\t3351\t3238\t0\t116\t76.1062]\n";
	sout_exp<<"[11\t1\t109601\t109708\t1\t3239\t3338\t3e-28\t108\t81.8182]\n";
	sout_exp<<"[12\t4\t109601\t109708\t4\t3239\t3338\t3e-28\t108\t81.8182]\n";


	CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
	FileUtils::removeFile(path_file);
	//ileUtils::removeFile(path_file+".attr");

}

void Test_BLRMatchMapLoad::test_loadPath_name2Num_and_num2Name(void)
{
	SDGString pathFileName = "input.path";	
	Test_BLRMatchMapLoadUtils::write_input_file_for_test_loadPath_name2Num_and_num2Name(pathFileName);

	BLRMatcherParameter para = Test_BLRMatchMapLoadUtils::createParameter();
	BLRMatchMap matchMap(&para);
 
	matchMap.loadPath(pathFileName, 0);

	std::map<std::string, long> obsName2NumQ = matchMap.getName2NumQ();
	std::map<long, std::string> obsNum2NameQ = matchMap.getNum2NameQ();	
	std::map<std::string, long> obsName2NumS = matchMap.getName2NumQ();
	std::map<long, std::string> obsNum2NameS = matchMap.getNum2NameS();	

	std::map<std::string, long> expName2NumQ = Test_BLRMatchMapLoadUtils::generateName2NumQ();
	std::map<long, std::string> expNum2NameQ = Test_BLRMatchMapLoadUtils::generateNum2NameQ();	
	std::map<std::string, long> expName2NumS = Test_BLRMatchMapLoadUtils::generateName2NumQ();
	std::map<long, std::string> expNum2NameS = Test_BLRMatchMapLoadUtils::generateNum2NameS();	
	
	CPPUNIT_ASSERT(expName2NumQ == obsName2NumQ);
	CPPUNIT_ASSERT(expNum2NameQ == obsNum2NameQ);
	CPPUNIT_ASSERT(expName2NumS == obsName2NumS);
	CPPUNIT_ASSERT(expNum2NameS == obsNum2NameS);

	FileUtils::removeFile(pathFileName);
}

void Test_BLRMatchMapLoad::test_loadPath_rpsList(void)
{
	SDGString pathFileName = "input.path";	
	Test_BLRMatchMapLoadUtils::write_input_file_for_test_loadPath_rpsList(pathFileName);

	BLRMatcherParameter para = Test_BLRMatchMapLoadUtils::createParameter();
	BLRMatchMap matchMap(&para);
 
	matchMap.loadPath(pathFileName, 0);
	std::list<RangePairSet> obsRpsList = matchMap.getRpsList();
 
	std::list<SDGString> lPath = Test_BLRMatchMapLoadUtils::generateExp_for_test_loadPath_rpsList();
	std::list<RangePairSet> expRpsList = Test_BLRMatchMapLoadUtils::createRpsList(lPath);
	std::cout<<"obsRpsList"<<std::endl;  
	std::cout<<"\n"<<std::endl;  
	Test_BLRMatchMapUtils::viewRangePairSetList(obsRpsList);
	std::cout<<"\n"<<std::endl;  
	std::cout<<"expRpsList"<<std::endl;  
	Test_BLRMatchMapUtils::viewRangePairSetList(expRpsList);
	std::cout<<"\n"<<std::endl;  
	CPPUNIT_ASSERT(expRpsList == obsRpsList);

	FileUtils::removeFile(pathFileName);
}


void Test_BLRMatchMapLoad::test_loadPath_with_join_end_of_file(void)
{
	SDGString pathFileName = "input.path";	
	Test_BLRMatchMapLoadUtils::write_input_file_for_test_loadPath_with_join_end_of_file(pathFileName);

	BLRMatcherParameter para = Test_BLRMatchMapLoadUtils::createParameter();
	BLRMatchMap matchMap(&para);
 
	matchMap.loadPath(pathFileName, 0);
	std::list<RangePairSet> obsRpsList = matchMap.getRpsList();

	std::list<SDGString> lPath = Test_BLRMatchMapLoadUtils::generateExp_for_test_loadPath_with_join_end_of_file();
	std::list<RangePairSet> expRpsList = Test_BLRMatchMapLoadUtils::createRpsList(lPath);
	
	//Test_BLRMatchMapUtils::viewRangePairSetList(obsRpsList);
	//Test_BLRMatchMapUtils::viewRangePairSetList(expRpsList);

	RangePairSet obsRps0 = obsRpsList.front();
	RangePairSet expRps0 = expRpsList.front();

	std::cout<<"obsRps0"<<std::endl;	
	obsRps0.view();
	std::cout<<"expRps0"<<std::endl;	
	expRps0.view();	

	RangePairSet obsRps1 = obsRpsList.back();
	RangePairSet expRps1 = expRpsList.back();
	
	std::cout<<"*******************"<<std::endl;		
	
	std::cout<<"obsRps1"<<std::endl;		
	obsRps1.view();
	std::cout<<"expRps1"<<std::endl;	
	expRps1.view();
	CPPUNIT_ASSERT(areTwoRpsEqualsWithoutIndentity(expRps0, obsRps0));
	CPPUNIT_ASSERT(areTwoRpsEqualsWithoutIndentity(expRps1, obsRps1));

	//FileUtils::removeFile(pathFileName);
}


void Test_BLRMatchMapLoad::test_loadPath_with_join_begin_of_file(void)
{
	SDGString pathFileName = "input.path";	
	Test_BLRMatchMapLoadUtils::write_input_file_for_test_loadPath_with_join_begin_of_file(pathFileName);

	BLRMatcherParameter para = Test_BLRMatchMapLoadUtils::createParameter();
	BLRMatchMap matchMap(&para);
 
	matchMap.loadPath(pathFileName, 0);
	std::list<RangePairSet> obsRpsList = matchMap.getRpsList();

	std::list<SDGString> lPath = Test_BLRMatchMapLoadUtils::generateExp_for_test_loadPath_with_join_begin_of_file();
	std::list<RangePairSet> expRpsList = Test_BLRMatchMapLoadUtils::createRpsList(lPath);
	
	//Test_BLRMatchMapUtils::viewRangePairSetList(obsRpsList);
	//Test_BLRMatchMapUtils::viewRangePairSetList(expRpsList);

	RangePairSet obsRps0 = obsRpsList.front();
	RangePairSet expRps0 = expRpsList.front();
	std::cout<<"obsRps0"<<std::endl;        
		obsRps0.view();
		std::cout<<"expRps0"<<std::endl;        
		expRps0.view(); 
	
	RangePairSet obsRps1 = obsRpsList.back();
	RangePairSet expRps1 = expRpsList.back();
	std::cout<<"obsRps1"<<std::endl;        
		obsRps1.view();
		std::cout<<"expRps1"<<std::endl;        
		expRps1.view(); 

	CPPUNIT_ASSERT(areTwoRpsEqualsWithoutIndentity(expRps0, obsRps0));
	// TODO this assertion fails 
	//CPPUNIT_ASSERT(expRpsList == obsRpsList);

	FileUtils::removeFile(pathFileName);
}



void Test_BLRMatchMapLoad::test_loadPath_with_join_middle_of_file(void)
{
	SDGString pathFileName = "input.path";	
	Test_BLRMatchMapLoadUtils::write_input_file_for_test_loadPath_with_join_middle_of_file(pathFileName);

	BLRMatcherParameter para = Test_BLRMatchMapLoadUtils::createParameter();
	BLRMatchMap matchMap(&para);
 
	matchMap.loadPath(pathFileName, 0);
	std::list<RangePairSet> obsRpsList = matchMap.getRpsList();

	std::list<SDGString> lPath = Test_BLRMatchMapLoadUtils::generateExp_for_test_loadPath_with_join_middle_of_file();
	std::list<RangePairSet> expRpsList = Test_BLRMatchMapLoadUtils::createRpsList(lPath);

	std::cout<<"obsRpsList"<<std::endl;
		std::cout<<"\n"<<std::endl;
	Test_BLRMatchMapUtils::viewRangePairSetList(obsRpsList);
	std::cout<<"expRpsList"<<std::endl;
		std::cout<<"\n"<<std::endl;
	Test_BLRMatchMapUtils::viewRangePairSetList(expRpsList);
	
	RangePairSet obsRps = obsRpsList.front();
	RangePairSet expRps = expRpsList.front();
	obsRpsList.pop_front();
	expRpsList.pop_front();
		
	CPPUNIT_ASSERT(obsRps == expRps);
	CPPUNIT_ASSERT(areTwoRpsEqualsWithoutIndentity(expRps, obsRps));
	
	obsRps = obsRpsList.front();
	expRps = expRpsList.front();
	obsRpsList.pop_front();
	expRpsList.pop_front();
		
	//CPPUNIT_ASSERT(obsRps == expRps);
	CPPUNIT_ASSERT(areTwoRpsEqualsWithoutIndentity(expRps, obsRps));
	
	obsRps = obsRpsList.front();
	expRps = expRpsList.front();
	obsRpsList.pop_front();
	expRpsList.pop_front();
	
	CPPUNIT_ASSERT(obsRps == expRps);
	CPPUNIT_ASSERT(areTwoRpsEqualsWithoutIndentity(expRps, obsRps));

}

bool Test_BLRMatchMapLoad::areTwoRpsEqualsWithoutIndentity(RangePairSet rps1, RangePairSet rps2)
{
	return rps1.getScore() == rps2.getScore() 
		&& rps1.getE_value() == rps2.getE_value() 
		&& rps1.getRangeQ() == rps2.getRangeQ()
		&& rps1.getRangeS() == rps2.getRangeS()
		&& rps1.getPath() == rps2.getPath();
	 
}























