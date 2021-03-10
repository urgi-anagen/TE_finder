#include "Test_BLRMatchMap.h"
#include "Test_BLRUtils.h"
#include "SDGString.h"
#include "FileUtils.h"
#include "BLRMatchMap.h"
#include "Range.h"
#include "RangePair.h"
#include "RangePairSet.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_BLRMatchMap );

void Test_BLRMatchMap::setUp()
{
}

void Test_BLRMatchMap::tearDown()
{
}

// clean_conflicts is never called in other submethod in MapPath.
void Test_BLRMatchMap::test_clean_conflicts(void){
	SDGString match_file = "match.align";
	Test_BLRUtils::writeInputFile();
	
	BLRMatcherThreadsParameter para = createParameter();
	BLRMatchMap matchMap(para);
		 
	matchMap.loadAlign(match_file);
	
	BLRMatchMap::MapAlign mapAlignBefore = matchMap.getMapAlign();
	

	/*std::cout<<" clean_conflicts before: "<<std::endl;
	Test_BLRUtils::viewMapAlign(mapAlignBefore);
	std::cout<<" "<<std::endl;*/


	matchMap.clean_conflicts();     
	BLRMatchMap::MapAlign obsMapAlign = matchMap.getMapAlign();  

	/*std::cout<<"clean_conflicts after: "<<std::endl;
	Test_BLRUtils::viewMapAlign(obsMapAlign);*/

	BLRMatchMap::MapAlign expMapAlign = Test_BLRUtils::createExpMapAlign_for_clean_conflicts();
	 
	/*	std::cout<<"expectation: "<<std::endl;
		Test_BLRUtils::viewMapAlign(expMapAlign);*/
 
	bool expComparaison = true;
	bool obsComparaison = (expMapAlign == obsMapAlign);
	
	CPPUNIT_ASSERT_EQUAL(expComparaison, obsComparaison);

	FileUtils::removeFile(match_file);
}

void Test_BLRMatchMap::test_mapPath(void){
		bool joiningParameter = true;
		bool cleanBefore = false;
		bool cleanAfter = true;
		bool merged = true;
		int verboseParameter = 0;
		
		std::ostringstream inputData;
		inputData<<"chunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100.00\n";
		inputData<<"chunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100.00\n";
		inputData<<"chunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100.00\n";
		inputData<<"chunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100.00\n";
		inputData<<"chunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80.00\n";
		inputData<<"chunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90.00\n";

		BLRMatcherThreadsParameter para = createParameter();
		BLRMatchMap matchMap(para);


		std::istringstream inputDataStream(inputData.str());
		matchMap.readAlign(inputDataStream);

		BLRMatchMap::MapAlign mapAlignBefore = matchMap.getMapAlign();

		matchMap.mapPath(joiningParameter, cleanBefore, cleanAfter, merged, verboseParameter);
		
		std::ostringstream obs;
		std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
		matchMap.writeRpsList(rps_list, obs);
		matchMap.writeRpsListAttribute(rps_list, obs);
		
		std::ostringstream exp;
		exp<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t401\t100\n";
		exp<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t401\t100\n";
		exp<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
		exp<<"1\tchunk1\t2001\t2500\trefTE_2\t2001\t2500\t0\t500\t100\n";
		exp<<"2\tchunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t801\t80\n";
		exp<<"2\tchunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t361\t90\n";
		exp<<"[1\tchunk1\t100\t2500\t-1\t0\t0\t0\t2302\t100]\n";
		exp<<"[2\tchunk1\t10000\t11500\trefTE_2\t1\t1500\t0\t1162\t82.8602]\n";

		CPPUNIT_ASSERT_EQUAL(exp.str(), obs.str());
}

void Test_BLRMatchMap::view_add_clean_path_all_S(void){


	BLRMatcherThreadsParameter para = createParameter();
	BLRMatchMap matchMap(para);


	std::list<RangePairSet> rpList = Test_BLRUtils::createRpList_for_test_add_clean_path_all_S();

	for (std::list<RangePairSet>::iterator i = rpList.begin();
		 i != rpList.end(); i++) {
		matchMap.add_clean_path_all_S(rpList, i, 1);
	}

	BLRMatchMap::MapPath obsMapPath = matchMap.getMapPath();
	//std::cout<<"add_clean_path_all_S Output "<<std::endl;
	//Test_BLRUtils::viewMapPath(obsMapPath);
	
}


void Test_BLRMatchMap::test_isOverlapFound_in_add_split_path(void){
		std::list<RangePairSet> rp_list = Test_BLRUtils::createRpListForTest_isOverloapFound_in_add_split_path();
		rp_list.sort(RangePair::lessIdentity);    
		 
		//std::cout<<" "<<std::endl;
		//std::cout<<"RpList isOverlapFound: "<<std::endl;
		//Test_BLRUtils::viewRangePairSetList(rp_list);
		

		BLRMatchMap::MapPath mapPath = Test_BLRUtils::createMapPathForTest_isOverlapFound_in_add_split_path();
						
		//std::cout<<"MapPath isOverlapFound: "<<std::endl;
		//Test_BLRUtils::viewMapPath(mapPath);
		
		double idTolerance = 2.0; 
		unsigned lenFilter = 20;

		bool expOverlap = true;
		bool obsOverlap = false;
		
		 for (std::list<RangePairSet>::iterator lrp_it=rp_list.begin(); lrp_it!=rp_list.end(); lrp_it++){
				obsOverlap = Test_BLRUtils::isOverlapFound_in_add_split_path(rp_list.begin(), mapPath, idTolerance, lenFilter);
				if (obsOverlap){
						break;
				}
		 }
		
		CPPUNIT_ASSERT_EQUAL(expOverlap, obsOverlap);
}   

void Test_BLRMatchMap::test_mapAlign_equality(void){
	BLRMatchMap::MapAlign mapAlign1 = Test_BLRUtils::createMapAlign_instance_for_test_mapAlign_Equality();
	BLRMatchMap::MapAlign mapAlign2 = Test_BLRUtils::createMapAlign_instance_for_test_mapAlign_Equality();
	bool expComparaison = true;
	bool obsComparaison = (mapAlign1 == mapAlign2);
	CPPUNIT_ASSERT_EQUAL(expComparaison, obsComparaison);
}


void Test_BLRMatchMap::test_merge_on_two_queries(void)
{
	 std::ostringstream inputData;
	 inputData<<"1\tchunk1\t10\t400\tTE1\t3351\t3238\t0\t116\t76.1062\n";
	 inputData<<"1\tchunk1\t700\t900\tTE1\t3239\t3338\t3e-28\t108\t81.8182\n";
	 inputData<<"2\tchunk1\t600\t1000\tTE1\t3202\t2921\t0\t269\t84.8361\n";
	 inputData<<"3\tchunk2\t10\t400\tTE1\t3351\t3238\t0\t116\t76.1062\n";
	 inputData<<"3\tchunk2\t700\t900\tTE1\t3239\t3338\t3e-28\t108\t81.8182\n";
	 inputData<<"4\tchunk2\t600\t1000\tTE1\t3202\t2921\t0\t269\t84.8361\n";

	 BLRMatcherThreadsParameter para = createParameter();
	 BLRMatchMap matchMap(para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);

	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRawRpsList();
	 matchMap.computeScoreWithLength(rpsList);
     matchMap.setRawRpsList(rpsList);
	 
	 matchMap.merge(0);

	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	 
	 std::ostringstream exp;
	 exp<<"1\tchunk1\t10\t400\tTE1\t3351\t3238\t0\t298\t76.1062\n";
	 exp<<"1\tchunk1\t600\t699\tTE1\t3202\t3133\t0\t85\t84.8361\n";
	 exp<<"1\tchunk1\t700\t900\tTE1\t3338\t3239\t3e-28\t164\t81.8182\n";
	 exp<<"1\tchunk1\t901\t1000\tTE1\t2990\t2921\t0\t85\t84.8361\n";
	 exp<<"2\tchunk2\t10\t400\tTE1\t3351\t3238\t0\t298\t76.1062\n";
	 exp<<"2\tchunk2\t600\t699\tTE1\t3202\t3133\t0\t85\t84.8361\n";
	 exp<<"2\tchunk2\t700\t900\tTE1\t3338\t3239\t3e-28\t164\t81.8182\n";
	 exp<<"2\tchunk2\t901\t1000\tTE1\t2990\t2921\t0\t85\t84.8361\n";
	 exp<<"[1\tchunk1\t10\t1000\t-1\t0\t0\t0\t632\t81.3452]\n";
	 exp<<"[2\tchunk2\t10\t1000\t-1\t0\t0\t0\t632\t81.3452]\n";

	 
	 CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}

void Test_BLRMatchMap::test_merge_second_fragment_include(void)
{
	 std::ostringstream inputData;
	 inputData<<"1\t1\t10\t400\t1\t3351\t3238\t0\t116\t76.1062\t116\n";
	 inputData<<"1\t1\t700\t900\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	 inputData<<"2\t1\t600\t1000\t1\t3202\t2921\t0\t269\t84.8361\t269\n";

	 BLRMatcherThreadsParameter para = createParameter();
	 BLRMatchMap matchMap(para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);

	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRawRpsList();
	 matchMap.computeScoreWithLength(rpsList);
    matchMap.setRawRpsList(rpsList);
	 
	 matchMap.merge(0);
	 
	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	
	 std::ostringstream exp;
	 exp<<"1\t1\t10\t400\t1\t3351\t3238\t0\t298\t76.1062\n";
	 exp<<"1\t1\t600\t699\t1\t3202\t3133\t0\t85\t84.8361\n";
	 exp<<"1\t1\t700\t900\t1\t3338\t3239\t3e-28\t164\t81.8182\n";
	 exp<<"1\t1\t901\t1000\t1\t2990\t2921\t0\t85\t84.8361\n";
	 exp<<"[1\t1\t10\t1000\t-1\t0\t0\t0\t632\t81.3452]\n";

		
	 //std::ofstream obsFile("obsFile");
	 //obsFile<<obs.str();

	 //std::ofstream expFile("expFile");
	 //expFile<<exp.str();
	 
	 CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}

void Test_BLRMatchMap::test_merge_overlap_right_on_second_fragment(void)
{
	 std::ostringstream inputData;
	 inputData<<"1\t1\t10\t400\t1\t3351\t3238\t0\t116\t76.1062\t116\n";
	 inputData<<"1\t1\t600\t900\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	 inputData<<"2\t1\t700\t1000\t1\t3202\t2921\t0\t269\t84.8361\t269\n";
		
	 BLRMatcherThreadsParameter para = createParameter();
	 BLRMatchMap matchMap(para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);

	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRawRpsList();
	 matchMap.computeScoreWithLength(rpsList);
    matchMap.setRawRpsList(rpsList);
	 
	 matchMap.merge(0);
	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	
	 std::ostringstream exp;
	 exp<<"1\t1\t10\t400\t1\t3351\t3238\t0\t298\t76.1062\n";
	 exp<<"1\t1\t600\t900\t1\t3338\t3239\t3e-28\t246\t81.8182\n";
	 exp<<"1\t1\t901\t1000\t1\t3014\t2921\t0\t85\t84.8361\n";
	 exp<<"[1\t1\t10\t1000\t-1\t0\t0\t0\t629\t79.5199]\n";
	 
 
	 //std::ofstream obsFile("obsFile");
	 //obsFile<<obs.str();

	 //std::ofstream expFile("expFile");
	 //expFile<<exp.str();
	 

	 CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}

void Test_BLRMatchMap::test_merge_overlap_left_on_second_fragment(void)
{   
	 std::ostringstream inputData;
	 inputData<<"1\t1\t10\t400\t1\t3351\t3238\t0\t116\t76.1062\t116\n";
	 inputData<<"1\t1\t600\t900\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	 inputData<<"2\t1\t500\t800\t1\t3202\t2921\t0\t269\t84.8361\t269\n";
	 
	 BLRMatcherThreadsParameter para = createParameter();
	 BLRMatchMap matchMap(para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);

	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRawRpsList();
	 matchMap.computeScoreWithLength(rpsList);
    matchMap.setRawRpsList(rpsList);
	 
	 matchMap.merge(0);
	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	
	 std::ostringstream exp;
	 exp<<"1\t1\t10\t400\t1\t3351\t3238\t0\t298\t76.1062\n";
	 exp<<"1\t1\t500\t599\t1\t3202\t3109\t0\t85\t84.8361\n";
	 exp<<"1\t1\t600\t900\t1\t3338\t3239\t3e-28\t246\t81.8182\n";
	 exp<<"[1\t1\t10\t900\t-1\t0\t0\t0\t629\t79.5199]\n";
	
	 
	 //std::ofstream obsFile("obsFile");
	 //obsFile<<obs.str();
	 
	 //std::ofstream expFile("expFile");
	 //expFile<<exp.str();
	 
	 CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}

void Test_BLRMatchMap::test_merge_overlap_first_and_second_fragment(void)
{
	 std::ostringstream inputData;
	 inputData<<"1\t1\t1\t400\t1\t3351\t3238\t0\t116\t76.1062\n";
	 inputData<<"1\t1\t600\t900\t1\t3239\t3338\t3e-28\t108\t81.8182\n";
	 inputData<<"2\t1\t200\t900\t2\t3202\t2921\t0\t269\t84.8361\n";
		
	 BLRMatcherThreadsParameter para = createParameter();
	 BLRMatchMap matchMap(para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);

	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRawRpsList();

	 matchMap.computeScoreWithLength(rpsList);
    matchMap.setRawRpsList(rpsList);
	 
	 matchMap.merge(0);
	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	
	 std::ostringstream exp;
	 exp<<"1\t1\t1\t199\t1\t3351\t3295\t0\t151\t76.1062\n";
	 exp<<"1\t1\t200\t900\t2\t3202\t2921\t0\t595\t84.8361\n";
	 exp<<"[1\t1\t1\t900\t-1\t0\t0\t0\t746\t83.0691]\n";
	 
	 CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}

void Test_BLRMatchMap::test_merge_on_mapPath_data(void)
{
	 std::ostringstream inputData;
	inputData<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100.00\n";
	inputData<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100.00\n";
	inputData<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100.00\n";
	inputData<<"2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100.00\n";
	inputData<<"3\tchunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80.00\n";
	inputData<<"3\tchunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90.00\n";


	 BLRMatcherThreadsParameter para = createParameter();
	 BLRMatchMap matchMap(para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);
	 
	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRawRpsList();

	 matchMap.computeScoreWithLength(rpsList);
    matchMap.setRawRpsList(rpsList);
	 
	 matchMap.merge(0);

	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	 
	 std::ostringstream exp;
	exp<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t401\t100\n";
	exp<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t401\t100\n";
	exp<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
	exp<<"1\tchunk1\t2001\t2500\trefTE_2\t2001\t2500\t0\t500\t100\n";
	exp<<"2\tchunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t801\t80\n";
	exp<<"2\tchunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t361\t90\n";
	exp<<"[1\tchunk1\t100\t2500\t-1\t0\t0\t0\t2302\t100]\n";
	exp<<"[2\tchunk1\t10000\t11500\trefTE_2\t1\t1500\t0\t1162\t82.8602]\n";

	 //std::ofstream expFile("mergeOnMapPathDataExpFile");
	 //expFile<<exp.str();
	 
	 //std::ofstream obsFile("mergeOnMapPathDataobsFile");
	 //obsFile<<obs.str();
	 
	 CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}

void Test_BLRMatchMap::test_merge_all_included(void)
{
	 std::ostringstream inputData;
	 inputData<<"1\t1\t10\t400\t1\t3351\t3238\t0\t116\t76.1062\t116\n";
	 inputData<<"1\t1\t800\t900\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	 inputData<<"2\t1\t5\t905\t1\t3202\t2921\t0\t269\t84.8361\t269\n";
	 
	 BLRMatcherThreadsParameter para = createParameter();
	 BLRMatchMap matchMap(para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);
	 
	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRawRpsList();
	 matchMap.computeScoreWithLength(rpsList);
    matchMap.setRawRpsList(rpsList);
	 
	 matchMap.merge(0);
	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	 
	 std::ostringstream exp;
	 exp<<"1\t1\t5\t905\t1\t3202\t2921\t0\t764\t84.8361\n";
	 exp<<"[1\t1\t5\t905\t-1\t0\t0\t0\t764\t84.8361]\n";
	 
			
	 //std::ofstream obsFile("obsFile");
	 //obsFile<<obs.str();
		
	 //std::ofstream expFile("expFile");
	 //expFile<<exp.str();
	 
	 CPPUNIT_ASSERT_EQUAL(exp.str(), obs.str());
}

void Test_BLRMatchMap::test_writeBED(void)
{
	std::ostringstream inputData;
	inputData<<"10\tCHR1v01212004\t100\t250\tTNAT1A\t36\t523\t2e-128\t460\t81.56\n";
	inputData<<"10\tCHR1v01212004\t800\t1000\tTNAT1A\t36\t523\t2e-128\t460\t81.56\n";
	inputData<<"20\tCHR1v01212004\t1050\t2000\tTNAT1A\t580\t480\t7e-78\t87\t79.14\n";

		BLRMatcherThreadsParameter para = createParameter();
		BLRMatchMap matchMap(para);
		
		std::istringstream inputDataStream(inputData.str());
		matchMap.readPath(inputDataStream, 0);

	std::ostringstream obs;
		std::list<RangePairSet> rps_list = matchMap.getRawRpsList();
	SDGString color = "255,127,0";	

	matchMap.writeBED(obs, rps_list, 0);

		std::ostringstream exp;
	exp<<"CHR1v01212004\t100\t1000\tTNAT1A,TNAT1A\t868\t+\t100\t1000\t255,127,0\t2\t151,201,\t0,700,\n";
	exp<<"CHR1v01212004\t1050\t2000\tTNAT1A\t87\t-\t1050\t2000\t255,127,0\n";

	/* 		
		std::ofstream obsFile("obsFileBED");
		obsFile<<obs.str();

		std::ofstream expFile("expFileBED");
		expFile<<exp.str();
		*/
	 CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}
