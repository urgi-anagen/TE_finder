#include "Test_BLRMatchMap.h"
#include "SDGString.h"
#include "FileUtils.h"
#include "BLRMatcherParameter.h"
#include "BLRMatchMap.h"
#include "Range.h"
#include "RangePair.h"
#include "RangePairSet.h"
#include "Test_BLRMatchMapUtils.h"

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
	Test_BLRMatchMapUtils::writeInputFile();
	
	BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter();
	BLRMatchMap matchMap(&para);
		 
	matchMap.load(0);     
	
	BLRMatchMap::MapAlign mapAlignBefore = matchMap.getMapAlign();
	

	/*std::cout<<" clean_conflicts before: "<<std::endl;
	Test_BLRMatchMapUtils::viewMapAlign(mapAlignBefore);
	std::cout<<" "<<std::endl;*/


	matchMap.clean_conflicts();     
	BLRMatchMap::MapAlign obsMapAlign = matchMap.getMapAlign();  

	/*std::cout<<"clean_conflicts after: "<<std::endl;
	Test_BLRMatchMapUtils::viewMapAlign(obsMapAlign);*/

	BLRMatchMap::MapAlign expMapAlign = Test_BLRMatchMapUtils::createExpMapAlign_for_clean_conflicts();  
	 
	/*	std::cout<<"expectation: "<<std::endl;
		Test_BLRMatchMapUtils::viewMapAlign(expMapAlign);*/
 
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

		BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
		BLRMatchMap matchMap(&para);

		std::istringstream inputDataStream(inputData.str());
		matchMap.readAlign(inputDataStream);

		BLRMatchMap::MapAlign mapAlignBefore = matchMap.getMapAlign();

		matchMap.mapPath(joiningParameter, cleanBefore, cleanAfter, merged, verboseParameter);
		
		std::ostringstream obs;
		//std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
		std::list<RangePairSet> rps_list = matchMap.getRpsList();
		matchMap.writeRpsList(rps_list, obs);
		matchMap.writeRpsListAttribute(rps_list, obs);
		
		std::ostringstream exp;
		exp<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t401\t100\n";
		exp<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t401\t100\n";
		exp<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
		exp<<"2\tchunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t1001\t80\n";
		exp<<"2\tchunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t401\t90\n";
		exp<<"3\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1001\t100\n";
		exp<<"[1\tchunk1\t100\t2000\trefTE_1\t1\t3000\t0\t1802\t100]\n";
		exp<<"[2\tchunk1\t10000\t11500\trefTE_2\t1\t1500\t0\t1402\t82.8602]\n";
		exp<<"[3\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1001\t100]\n";


		CPPUNIT_ASSERT_EQUAL(exp.str(), obs.str());
}

void Test_BLRMatchMap::view_add_clean_path_all_S(void){
		SDGString match_file = "match.align";
		Test_BLRMatchMapUtils::writeInputFile();
 
		BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
		BLRMatchMap matchMap(&para);
		matchMap.load(0);     
		std::list<RangePairSet> rpList = Test_BLRMatchMapUtils::createRpList_for_test_add_clean_path_all_S();
		
			 
		//std::cout<<"add_clean_path_all_S Input "<<std::endl;
		//for(std::list<RangePairSet>::iterator i=rpList.begin();
		//   i!=rpList.end();i++){
		//    i->view();
		//}
		
	 for(std::list<RangePairSet>::iterator i=rpList.begin();
			 i!=rpList.end();i++){
				matchMap.add_clean_path_all_S(rpList, i, 1);
		}
		
		//BLRMatchMap::MapPath obsMapPath = matchMap.getMapPath(); 
		//std::cout<<"add_clean_path_all_S Output "<<std::endl;
		//Test_BLRMatchMapUtils::viewMapPath(obsMapPath);
			
		FileUtils::removeFile(match_file);
}

void Test_BLRMatchMap::view_split_path(void){
		bool joiningParameter = true;
		bool cleanBefore = false;
		bool cleanAfter = true;
		int verboseParameter = 0;
		
		SDGString match_file = "match.align";
		Test_BLRMatchMapUtils::writeInputFile();
		
		BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
		BLRMatchMap matchMap(&para);
		 
		matchMap.load(0);     
		BLRMatchMap::MapAlign mapAlignBefore = matchMap.getMapAlign();
	
		matchMap.mapPathJoinAndComputeScoreWithLengthOnly(joiningParameter, cleanBefore, cleanAfter, verboseParameter);
		
			 
		//std::cout<<"mapPath Before "<<std::endl;
		//Test_BLRMatchMapUtils::viewMapAlign(mapAlignBefore);
		//std::cout<<" "<<std::endl;
		

		BLRMatchMap::MapPath inMapPath = matchMap.getMapPath(); 
		
		//std::cout<<"in "<<std::endl;
		//Test_BLRMatchMapUtils::viewMapPath(inMapPath);
		
		matchMap.split_path();
					 
		BLRMatchMap::MapPath obsMapPath = matchMap.getMapPath(); 
				
		//std::cout<<"Obs "<<std::endl;
		//Test_BLRMatchMapUtils::viewMapPath(obsMapPath);
		

		// TODO unable ot understand test assertion error: data obs and exp are EQUALS
		// TODO assertion error from DATA structure ?
			
		FileUtils::removeFile(match_file);
		
}

void Test_BLRMatchMap::test_isOverlapFound_in_add_split_path(void){
		std::list<RangePairSet> rp_list = Test_BLRMatchMapUtils::createRpListForTest_isOverloapFound_in_add_split_path();
		rp_list.sort(RangePair::lessIdentity);    
		 
		//std::cout<<" "<<std::endl;
		//std::cout<<"RpList isOverlapFound: "<<std::endl;
		//Test_BLRMatchMapUtils::viewRangePairSetList(rp_list);
		

		BLRMatchMap::MapPath mapPath = Test_BLRMatchMapUtils::createMapPathForTest_isOverlapFound_in_add_split_path();
						
		//std::cout<<"MapPath isOverlapFound: "<<std::endl;
		//Test_BLRMatchMapUtils::viewMapPath(mapPath);
		
		double idTolerance = 2.0; 
		unsigned lenFilter = 20;

		bool expOverlap = true;
		bool obsOverlap = false;
		
		 for (std::list<RangePairSet>::iterator lrp_it=rp_list.begin(); lrp_it!=rp_list.end(); lrp_it++){
				obsOverlap = BLRMatchMap::isOverlapFound_in_add_split_path(rp_list.begin(), mapPath, idTolerance, lenFilter);
				if (obsOverlap){
						break;
				}
		 }
		
		CPPUNIT_ASSERT_EQUAL(expOverlap, obsOverlap);
}   

void Test_BLRMatchMap::test_mapAlign_equality(void){
	BLRMatchMap::MapAlign mapAlign1 = Test_BLRMatchMapUtils::createMapAlign_instance_for_test_mapAlign_Equality(); 
	BLRMatchMap::MapAlign mapAlign2 = Test_BLRMatchMapUtils::createMapAlign_instance_for_test_mapAlign_Equality();
	bool expComparaison = true;
	bool obsComparaison = (mapAlign1 == mapAlign2);
	CPPUNIT_ASSERT_EQUAL(expComparaison, obsComparaison);
}
/*
// mapPath compraison problem. diagnostic: compare uniquely list behind keys
void Test_BLRMatchMap::test_mapPath_comparaison_diagnostic_compare_rangePairSetList_list_comparaison_key_1_2(void){
		BLRMatchMap::MapPath mapPathAfterJoin = Test_BLRMatchMapUtils::createMapPath_afterJoin();
		BLRMatchMap::MapPath expMapPath = Test_BLRMatchMapUtils::createExpMapPathForTest_mapPath_comparaison_diagnostic();
		
		std::list<RangePairSet>& rpsListFromCreateMapPathAfterJoin = mapPathAfterJoin[BLRMatchMap::Key(1,2)];     
		std::list<RangePairSet>& rpsListFromMapPathCreation = expMapPath[BLRMatchMap::Key(1,2)];     
		
		CPPUNIT_ASSERT(rpsListFromCreateMapPathAfterJoin == rpsListFromMapPathCreation);
}

// mapPath compraison problem. diagnostic: compare uniquely list behind keys
void Test_BLRMatchMap::test_mapPath_comparaison_diagnostic_compare_rangePairSetList_list_comparaison_key_1_3(void){
		BLRMatchMap::MapPath mapPathAfterJoin = Test_BLRMatchMapUtils::createMapPath_afterJoin();
		BLRMatchMap::MapPath expMapPath = Test_BLRMatchMapUtils::createExpMapPathForTest_mapPath_comparaison_diagnostic();
		
		std::list<RangePairSet>& rpsListFromCreateMapPathAfterJoin = mapPathAfterJoin[BLRMatchMap::Key(1,3)];     
		std::list<RangePairSet>& rpsListFromMapPathCreation = expMapPath[BLRMatchMap::Key(1,3)];     
		
		CPPUNIT_ASSERT(rpsListFromCreateMapPathAfterJoin == rpsListFromMapPathCreation);
}
*/
void Test_BLRMatchMap::test_reComputeScoreWithLength_on_a_match_with_one_match_part(void){
		std::list<RangePairSet> inRpsList = Test_BLRMatchMapUtils::createInputRpsList_for_test_reComputeScoreWithLength_on_a_match_with_one_match_part();   
		std::list<RangePairSet> expRpsList = Test_BLRMatchMapUtils::createExpRpsList_for_test_reComputeScoreWithLength_on_a_match_with_one_match_part();  	
		
		BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter();
		BLRMatchMap matchMap(&para);
			
		//std::cout<<" "<<std::endl;    
		//std::cout<<"In "<<std::endl;    
		
		//Test_BLRMatchMapUtils::viewRangePairSetList(inRpsList);
		
		matchMap.computeScoreWithLength(inRpsList);

		std::list<RangePairSet> obsRpsList = inRpsList;
			 
		//std::cout<<" "<<std::endl;    
		//std::cout<<"Obs "<<std::endl;    
		//Test_BLRMatchMapUtils::viewRangePairSetList(obsRpsList);
 
		//std::cout<<" "<<std::endl;    
		//std::cout<<"Exp "<<std::endl;    
		//Test_BLRMatchMapUtils::viewRangePairSetList(expRpsList);
			

		CPPUNIT_ASSERT(expRpsList == obsRpsList);
}

void Test_BLRMatchMap::test_reComputeScoreWithLength_on_a_match_with_two_match_part(void){
		std::list<RangePairSet> inRpsList = Test_BLRMatchMapUtils::createInputRpsList_for_test_reComputeScoreWithLength_on_a_match_with_two_match_part();   
		std::list<RangePairSet> expRpsList = Test_BLRMatchMapUtils::createExpRpsList_for_test_reComputeScoreWithLength_on_a_match_with_two_match_part();  	
		 
		//std::cout<<" "<<std::endl;    
		//std::cout<<"In "<<std::endl;    
		
		//Test_BLRMatchMapUtils::viewRangePairSetList(inRpsList);
		
 
		BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter();
		BLRMatchMap matchMap(&para);
	
		matchMap.computeScoreWithLength(inRpsList);
		
		std::list<RangePairSet> obsRpsList = inRpsList;
		
		//std::cout<<" "<<std::endl;    
		//std::cout<<"Obs "<<std::endl;    
		//Test_BLRMatchMapUtils::viewRangePairSetList(obsRpsList);
 
		//std::cout<<" "<<std::endl;    
		//std::cout<<"Exp "<<std::endl;    
		//Test_BLRMatchMapUtils::viewRangePairSetList(expRpsList);
		
		CPPUNIT_ASSERT(expRpsList == obsRpsList);
}



// keep this test ?
void Test_BLRMatchMap::test_merge_on_two_queries(void)
{
	 std::ostringstream inputData;
	 inputData<<"1\t1\t10\t400\t1\t3351\t3238\t0\t116\t76.1062\t116\n";
	 inputData<<"1\t1\t700\t900\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	 inputData<<"2\t1\t600\t1000\t1\t3202\t2921\t0\t269\t84.8361\t269\n";
	 inputData<<"3\t2\t10\t400\t1\t3351\t3238\t0\t116\t76.1062\t116\n";
	 inputData<<"3\t2\t700\t900\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	 inputData<<"4\t2\t600\t1000\t1\t3202\t2921\t0\t269\t84.8361\t269\n";

	 BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
	 BLRMatchMap matchMap(&para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);

	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRpsList();
	 matchMap.computeScoreWithLength(rpsList);
	 matchMap.setRpsList(rpsList);
	 
	 matchMap.merge(0);

	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	 
	 std::ostringstream exp;
	 exp<<"1\t1\t10\t400\t1\t3351\t3238\t0\t391\t76.1062\n";
	 exp<<"1\t1\t600\t699\t1\t3202\t3133\t0\t100\t84.8361\n";
	 exp<<"1\t1\t700\t900\t1\t3338\t3239\t3e-28\t201\t81.8182\n";
	 exp<<"1\t1\t901\t1000\t1\t2990\t2921\t0\t100	84.8361\n";
	 exp<<"2\t2\t10\t400\t1\t3351\t3238\t0\t391\t76.1062\n";
	 exp<<"2\t2\t600\t699\t1\t3202\t3133\t0\t100\t84.8361\n";
	 exp<<"2\t2\t700\t900\t1\t3338\t3239\t3e-28\t201\t81.8182\n";
	 exp<<"2\t2\t901\t1000\t1\t2990\t2921\t0\t100\t84.8361\n";
	 exp<<"[1\t1\t10\t1000\t-1\t0\t0\t0\t792\t0]\n";
	 exp<<"[2\t2\t10\t1000\t-1\t0\t0\t0\t792\t0]\n";

		
	 //std::ofstream obsFile("obsFile");
	 //obsFile<<obs.str();

	 //std::ofstream expFile("expFile");
	 //expFile<<exp.str();
	 
	 CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}

void Test_BLRMatchMap::test_merge_second_fragment_include(void)
{
	 std::ostringstream inputData;
	 inputData<<"1\t1\t10\t400\t1\t3351\t3238\t0\t116\t76.1062\t116\n";
	 inputData<<"1\t1\t700\t900\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	 inputData<<"2\t1\t600\t1000\t1\t3202\t2921\t0\t269\t84.8361\t269\n";

	 BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
	 BLRMatchMap matchMap(&para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);

	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRpsList();
	 matchMap.computeScoreWithLength(rpsList);
	 matchMap.setRpsList(rpsList);
	 
	 matchMap.merge(0);
	 
	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	
	 std::ostringstream exp;
	 exp<<"1\t1\t10\t400\t1\t3351\t3238\t0\t391\t76.1062\n";
	 exp<<"1\t1\t600\t699\t1\t3202\t3133\t0\t100\t84.8361\n";
	 exp<<"1\t1\t700\t900\t1\t3338\t3239\t3e-28\t201\t81.8182\n";
	 exp<<"1\t1\t901\t1000\t1\t2990\t2921\t0\t100	84.8361\n";
	 exp<<"[1\t1\t10\t1000\t-1\t0\t0\t0\t792\t0]\n";

		
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
		
	 BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
	 BLRMatchMap matchMap(&para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);

	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRpsList();
	 matchMap.computeScoreWithLength(rpsList);
	 matchMap.setRpsList(rpsList);
	 
	 matchMap.merge(0);
	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	
	 std::ostringstream exp;
	 exp<<"1\t1\t10\t400\t1\t3351\t3238\t0\t391\t76.1062\n";
	 exp<<"1\t1\t600\t900\t1\t3338\t3239\t3e-28\t301\t81.8182\n";
	 exp<<"1\t1\t901\t1000\t1\t3014\t2921\t0\t100\t84.8361\n";
	 exp<<"[1\t1\t10\t1000\t-1\t0\t0\t0\t792\t0]\n";
	 
 
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
	 
	 BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
	 BLRMatchMap matchMap(&para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);

	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRpsList();
	 matchMap.computeScoreWithLength(rpsList);
	 matchMap.setRpsList(rpsList);
	 
	 matchMap.merge(0);
	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	
	 std::ostringstream exp;
	 exp<<"1\t1\t10\t400\t1\t3351\t3238\t0\t391\t76.1062\n"; 
	 exp<<"1\t1\t500\t599\t1\t3202\t3109\t0\t100\t84.8361\n";
	 exp<<"1\t1\t600\t900\t1\t3338\t3239\t3e-28\t301\t81.8182\n";
	 exp<<"[1\t1\t10\t900\t-1\t0\t0\t0\t792\t0]\n";
	
	 
	 //std::ofstream obsFile("obsFile");
	 //obsFile<<obs.str();
	 
	 //std::ofstream expFile("expFile");
	 //expFile<<exp.str();
	 
	 CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}

void Test_BLRMatchMap::test_merge_overlap_first_and_second_fragment(void)
{
	 std::ostringstream inputData;
	 inputData<<"1\t1\t1\t400\t1\t3351\t3238\t0\t116\t76.1062\t116\n";
	 inputData<<"1\t1\t600\t900\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	 inputData<<"2\t1\t200\t900\t2\t3202\t2921\t0\t269\t84.8361\t269\n";
		
	 BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
	 BLRMatchMap matchMap(&para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);

	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRpsList();

	 matchMap.computeScoreWithLength(rpsList);
	 matchMap.setRpsList(rpsList);
	 
	 matchMap.merge(0);
	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	
	 std::ostringstream exp;
	 exp<<"1\t1\t1\t400\t1\t3351\t3238\t0\t400\t76.1062\n";
	 exp<<"1\t1\t401\t599\t2\t3121\t3042\t0\t199\t84.8361\n";
	 exp<<"1\t1\t600\t900\t1\t3338\t3239\t3e-28\t301\t81.8182\n";
	 exp<<"[1\t1\t1\t900\t-1\t0\t0\t0\t900\t0]\n";
	
	 
	 //std::ofstream obsFile("obsFile");
	 //obsFile<<obs.str();
		
	 //std::ofstream expFile("expFile");
	 //expFile<<exp.str();
	 
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


	 BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
	 BLRMatchMap matchMap(&para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);
	 
	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRpsList();

	 matchMap.computeScoreWithLength(rpsList);
	 matchMap.setRpsList(rpsList);
	 
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
	exp<<"2\tchunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t1001\t80\n";
	exp<<"2\tchunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t401\t90\n";
	exp<<"[1\tchunk1\t100\t2500\t-1\t0\t0\t0\t2302\t0]\n";
	exp<<"[2\tchunk1\t10000\t11500\trefTE_2\t1\t1500\t0\t1402\t82.8602]\n";

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
	 
	 BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
	 BLRMatchMap matchMap(&para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);
	 
	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRpsList();
	 matchMap.computeScoreWithLength(rpsList);
	 matchMap.setRpsList(rpsList);
	 
	 matchMap.merge(0);
	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(rps_list, obs);
	 matchMap.writeRpsListAttribute(rps_list, obs);
	 
	 std::ostringstream exp;
	 exp<<"1\t1\t5\t905\t1\t3202\t2921\t0\t901\t84.8361\n";
	 exp<<"[1\t1\t5\t905\t-1\t0\t0\t0\t901\t0]\n";
	 
			
	 //std::ofstream obsFile("obsFile");
	 //obsFile<<obs.str();
		
	 //std::ofstream expFile("expFile");
	 //expFile<<exp.str();
	 
	 CPPUNIT_ASSERT_EQUAL(exp.str(), obs.str());
}

void Test_BLRMatchMap::test_clusterizeOverlapingRps(void)
{
	 std::ostringstream inputData;
	 inputData<<"10\t1\t10\t400\t1\t3351\t3238\t0\t116\t76.1062\t116\n";
	 inputData<<"10\t1\t700\t1000\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	 inputData<<"20\t1\t20\t100\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	 inputData<<"30\t1\t600\t900\t1\t3202\t2921\t0\t269\t84.8361\t269\n";
	 inputData<<"40\t1\t800\t1200\t1\t3202\t2921\t0\t269\t84.8361\t269\n";

	 BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
	 BLRMatchMap matchMap(&para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);

	 // the score is the length 
	 std::list<RangePairSet> rpsList = matchMap.getRpsList();
	 matchMap.computeScoreWithLength(rpsList);

	 Graph <unsigned long> graph = matchMap.clusterizeOverlapingRps(rpsList,0); 
		 
	 std::ostringstream obs;
	 graph.toStream(obs);
		
	 std::ostringstream exp;
	 exp<<"node:0 ->10 next: 2(30) 3(40)\n";
	 exp<<"node:1 ->20 next:\n";
	 exp<<"node:2 ->30 next: 0(10)\n";
	 exp<<"node:3 ->40 next: 0(10)\n";
 
		
	 //std::ofstream obsFile("obsFile");
	 //obsFile<<obs.str();

	 //std::ofstream expFile("expFile");
	 //expFile<<exp.str();
	 
	 CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}

void Test_BLRMatchMap::test_mergeOnCluster(void)
{
	 std::ostringstream inputData;
	 inputData<<"10\t1\t10\t400\t1\t3351\t3238\t0\t116\t76.1062\t116\n";
	 inputData<<"10\t1\t700\t1000\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	 inputData<<"20\t1\t20\t100\t1\t3239\t3338\t3e-28\t108\t81.8182\t108\n";
	 inputData<<"30\t1\t600\t900\t1\t3202\t2921\t0\t269\t84.8361\t269\n";
	 inputData<<"40\t1\t800\t1200\t1\t3202\t2921\t0\t269\t84.8361\t269\n";

	
	 BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
	 BLRMatchMap matchMap(&para);
		
	 std::istringstream inputDataStream(inputData.str());
	 matchMap.readPath(inputDataStream, 0);

	 // rpsList is, for the test a cluster 
	 std::list<RangePairSet> rpsList = matchMap.getRpsList();
	 // the score is the length 
	 matchMap.computeScoreWithLength(rpsList);

	 Graph<unsigned long> graph;
	 graph.add_edge(10,30);
	 graph.add_edge(10,40);
	 graph.add_node(20);

	 std::list<RangePairSet> mergedRpsList = matchMap.mergeOnCluster(rpsList, graph, 0);
		
	 std::ostringstream obs;
	 std::list<RangePairSet> rps_list = matchMap.getRpsListFromMapPath();
	 matchMap.writeRpsList(mergedRpsList, obs);
	 matchMap.writeRpsListAttribute(mergedRpsList, obs);
	
	 std::ostringstream exp;
	 exp<<"1\t1\t10\t400\t1\t3351\t3238\t0\t391\t76.1062\n";
	 exp<<"1\t1\t600\t699\t1\t3202\t3109\t0\t100\t84.8361\n";
	 exp<<"1\t1\t700\t1000\t1\t3338\t3239\t3e-28\t301\t81.8182\n";
	 exp<<"1\t1\t1001\t1200\t1\t3061\t2921\t0\t200\t84.8361\n";
	 exp<<"2\t1\t20\t100\t1\t3239\t3338\t3e-28\t81\t81.8182\n";
	 exp<<"[1\t1\t10\t1200\t-1\t0\t0\t0\t992\t0]\n";
	 exp<<"[2\t1\t20\t100\t1\t3239\t3338\t3e-28\t81\t81.8182]\n";
	 
	 //std::ofstream obsFile("obsFile");
	 //obsFile<<obs.str();

	 //std::ofstream expFile("expFile");
	 //expFile<<exp.str();
	
	 CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}

void Test_BLRMatchMap::test_writeBED(void)
{
	std::ostringstream inputData;
	inputData<<"10\tCHR1v01212004\t100\t250\tTNAT1A\t36\t523\t2e-128\t460\t81.56\n";
	inputData<<"10\tCHR1v01212004\t800\t1000\tTNAT1A\t36\t523\t2e-128\t460\t81.56\n";
	inputData<<"20\tCHR1v01212004\t1050\t2000\tTNAT1A\t580\t480\t7e-78\t87\t79.14\n";

		BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter(); 
		BLRMatchMap matchMap(&para);
		
		std::istringstream inputDataStream(inputData.str());
		matchMap.readPath(inputDataStream, 0);

	std::ostringstream obs;
		std::list<RangePairSet> rps_list = matchMap.getRpsList();
	SDGString color = "255,127,0";	

	matchMap.writeBED(obs, rps_list, color, 0);

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
