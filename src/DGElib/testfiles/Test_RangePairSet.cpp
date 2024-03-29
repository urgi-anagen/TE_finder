#include "Test_RangePairSet.h"
#include "Test_RangePairSetUtils.h"
#include "RangePairSet.h"


CPPUNIT_TEST_SUITE_REGISTRATION( Test_RangePairSet );

void Test_RangePairSet::setUp()
{
}

void Test_RangePairSet::tearDown()
{
}

void Test_RangePairSet::test_inserted_false(void)
{

    RangePairSet rps1 = Test_RangePairSetUtils::createRangePairSet1ForTest_inserted();
    RangePairSet rps2 = Test_RangePairSetUtils::createRangePairSet2ForTest_inserted();
    /*
    std::cout<<" "<<std::endl;
    std::cout<<"rps 1: "<<std::endl;
    rps1.view();
    
    std::cout<<" "<<std::endl;
    std::cout<<"rps 2: "<<std::endl;
    rps2.view();
    */
    bool expInserted = false;
    bool obsInserted = rps1.inserted(rps2);

    CPPUNIT_ASSERT_EQUAL(expInserted,obsInserted);
}


void Test_RangePairSet::test_inserted_true(void)
{
    RangePairSet rps1 = Test_RangePairSetUtils::createRangePairSet1ForTest_inserted_true();
    RangePairSet rps2 = Test_RangePairSetUtils::createRangePairSet2ForTest_inserted_true();
    /*
    std::cout<<" "<<std::endl;
    std::cout<<"rps 1: "<<std::endl;
    rps1.view();
    
    std::cout<<" "<<std::endl;
    std::cout<<"rps 2: "<<std::endl;
    rps2.view();
    */
    bool expInserted = true;
    bool obsInserted = rps1.inserted(rps2);

    CPPUNIT_ASSERT_EQUAL(expInserted,obsInserted);
}


void Test_RangePairSet::test_overlapQ(void)
{
    RangePairSet rps1 = Test_RangePairSetUtils::createRangePairSet1ForTest_overlapQ();
    RangePairSet rps2 = Test_RangePairSetUtils::createRangePairSet2ForTest_overlapQ();
    
    /*  
    std::cout<<" "<<std::endl;
    std::cout<<"rps 1: "<<std::endl;
    rps1.view();

    std::cout<<" "<<std::endl;
    std::cout<<"rps 2: "<<std::endl;
    rps2.view();
    */
    
    bool expOverlapQ = true;
    bool obsOverlapQ = rps1.overlapQ(rps2);

    CPPUNIT_ASSERT_EQUAL(expOverlapQ,obsOverlapQ);
}


void Test_RangePairSet::test_overlapQ_rps2_include_in_rps1_but_no_overlap(void)
{
	RangePairSet rps1 = Test_RangePairSetUtils::createRangePairSet1ForTest_overlapQ_rps2_include_in_rps1_but_no_overlap();
	RangePairSet rps2 = Test_RangePairSetUtils::createRangePairSet2ForTest_overlapQ_rps2_include_in_rps1_but_no_overlap();
	/*
	std::cout<<" "<<std::endl;
 	std::cout<<"rps 1: "<<std::endl;
    	rps1.view();

    	std::cout<<" "<<std::endl;
    	std::cout<<"rps 2: "<<std::endl;
    	rps2.view();
	*/
	bool expOverlapQ = true;
	bool obsOverlapQ= rps1.overlapQ(rps2);

    	CPPUNIT_ASSERT_EQUAL(expOverlapQ,obsOverlapQ);
}

void Test_RangePairSet::test_overlapQ_length(void)
{
    RangePairSet rps1 = Test_RangePairSetUtils::createRangePairSet1ForTest_overlapQ_length();
    RangePairSet rps2 = Test_RangePairSetUtils::createRangePairSet2ForTest_overlapQ_length();
    
    /*      
    std::cout<<" "<<std::endl;
    std::cout<<"rps 1: "<<std::endl;
    rps1.view();

    std::cout<<" "<<std::endl;
    std::cout<<"rps 2: "<<std::endl;
    rps2.view();
    
    */
    bool expLength = 414;
    bool obsLength = rps1.overlapQ(rps2);

    CPPUNIT_ASSERT_EQUAL(expLength,obsLength);
}

void Test_RangePairSet::test_equality(void)
{
	RangePairSet rps1 = Test_RangePairSetUtils::createRangePairSet1For_test_equality();
	RangePairSet rps2 = Test_RangePairSetUtils::createRangePairSet2For_test_equality();
	CPPUNIT_ASSERT_EQUAL(rps1, rps2);
}

void Test_RangePairSet::test_equality_not_equal(void)
{
	RangePairSet rps1 = Test_RangePairSetUtils::createRangePairSet1For_test_equality_not_equal();
	RangePairSet rps2 = Test_RangePairSetUtils::createRangePairSet2For_test_equality_not_equal();
	CPPUNIT_ASSERT(rps1 != rps2);
}

void Test_RangePairSet::test_not_equality_on_match_part(void)
{
	RangePairSet rps1 = Test_RangePairSetUtils::createRangePairSet1For_test_equality_on_match_part();
	RangePairSet rps2 = Test_RangePairSetUtils::createRangePairSet2For_test_equality_on_match_part();
	std::list<RangePair> rps1Path = rps1.getPath();
	std::list<RangePair> rps2Path = rps2.getPath();
	CPPUNIT_ASSERT(rps1 != rps2);
}





void Test_RangePairSet::test_overlapQ_length_length(void)
{
	RangePairSet rps1 = Test_RangePairSetUtils::createRangePairSet1ForTest_overlapQ_length();
	RangePairSet rps2 = Test_RangePairSetUtils::createRangePairSet2ForTest_overlapQ_length();

	unsigned obsOverlap = rps2.overlapQ_length(rps1);
	unsigned expOverlap = 414;

	CPPUNIT_ASSERT(expOverlap == obsOverlap);		 
}

void Test_RangePairSet::test_overlapQ_length_length_symetric(void)
{
	RangePairSet rps1 = Test_RangePairSetUtils::createRangePairSet1ForTest_overlapQ_length();
	RangePairSet rps2 = Test_RangePairSetUtils::createRangePairSet2ForTest_overlapQ_length();

	unsigned obsOverlap = rps1.overlapQ_length(rps2);
	unsigned expOverlap = 414;

	CPPUNIT_ASSERT(expOverlap == obsOverlap);		 
}

void Test_RangePairSet::view_diffQ_for_test_mergeQ(void)
{
        RangePairSet rps1 = Test_RangePairSetUtils::createRangePairSet1ForTest_diffQ();
        RangePairSet rps2 = Test_RangePairSetUtils::createRangePairSet2ForTest_diffQ();
   	/*	
	std::cout<<"rps first "<<std::endl;
	rps1.view();
	std::cout<<"  "<<std::endl;
	std::cout<<"rps second "<<std::endl;
	rps2.view();
	*/	
	rps1.diffQ(rps2);
	/*
	std::cout<<"obs rps first "<<std::endl;
	rps1.view();
	std::cout<<"  "<<std::endl;
	std::cout<<"obs rps second "<<std::endl;
	rps2.view();
	*/
	
}

void Test_RangePairSet::test_orientSubjects_on_plus_strand(void)
{
	RangePairSet inRps = Test_RangePairSetUtils::generateInputs_for_test_orientSubjects_on_plus_strand();
	RangePairSet expRps = Test_RangePairSetUtils::generateOutputs_for_test_orientSubjects_on_plus_strand();
	inRps.orientSubjects();
	RangePairSet obsRps = inRps;
	CPPUNIT_ASSERT(expRps == obsRps);
}

void Test_RangePairSet::test_orientSubjects_on_minus_strand(void)
{
	RangePairSet inRps = Test_RangePairSetUtils::generateInputs_for_test_orientSubjects_on_minus_strand();
	RangePairSet expRps = Test_RangePairSetUtils::generateOutputs_for_test_orientSubjects_on_minus_strand();
	inRps.orientSubjects();
	RangePairSet obsRps = inRps;
	CPPUNIT_ASSERT(expRps == obsRps);

}

void Test_RangePairSet::test_orientSubjects_path_size_is_1(void)
{
	RangePairSet inRps = Test_RangePairSetUtils::generateInputs_for_test_orientSubjects_path_size_is_1();
	RangePairSet expRps = Test_RangePairSetUtils::generateOutputs_for_test_orientSubjects_path_size_is_1();
	inRps.orientSubjects();
	RangePairSet obsRps = inRps;
	CPPUNIT_ASSERT(expRps == obsRps);
}

void Test_RangePairSet::test_orientSubjects_orient_more_than_1_fragment(void)
{
	RangePairSet inRps = Test_RangePairSetUtils::generateInputs_for_test_orientSubjects_orient_more_than_1_fragment();
	RangePairSet expRps = Test_RangePairSetUtils::generateOutputs_for_test_orientSubjects_orient_more_than_1_fragment();
	inRps.orientSubjects();
	RangePairSet obsRps = inRps;
	CPPUNIT_ASSERT(expRps == obsRps);
}

void Test_RangePairSet::test_orientSubjects_no_fragments_to_orient(void)
{
	RangePairSet inRps = Test_RangePairSetUtils::generateInputs_for_test_orientSubjects_no_fragments_to_orient();
	RangePairSet expRps = Test_RangePairSetUtils::generateOutputs_for_test_orientSubjects_no_fragments_to_orient();
	inRps.orientSubjects();
	RangePairSet obsRps = inRps;
	CPPUNIT_ASSERT(expRps == obsRps);
}
