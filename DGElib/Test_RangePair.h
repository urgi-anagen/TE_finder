/**
 * \file Test_RangePair.h
 * \brief Unitary tests for class RangePair
 */

#ifndef TEST_RANGEPAIR_H
#define TEST_RANGEPAIR_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "RangePair.h"

/**
 * \class Test_RangePair
 * \brief Unitary tests for the class RangePair
 */
class Test_RangePair: public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE(Test_RangePair);
	CPPUNIT_TEST(test_writetxt);
	CPPUNIT_TEST(test_diffQ_rp1_start_modified);
	CPPUNIT_TEST(test_getRangeQ_affectation);
	CPPUNIT_TEST(test_isPlusStrand_reverse_strand);
	//CPPUNIT_TEST(test_mergeQ);
 	CPPUNIT_TEST_SUITE_END();

	public:
     void setUp();
	 void tearDown();

	protected:
	  void test_set_string( void );
	  void test_writetxt( void );
	  void test_diffQ_rp1_start_modified(void);
	  void test_getRangeQ_affectation(void);
	  void test_isPlusStrand_reverse_strand(void);
	  //void test_mergeQ(void);
	private:
	  RangePair test_createRangePairSet1ForTest_mergeQ(void);
	  RangePair test_createRangePairSet2ForTest_mergeQ(void);
};

#endif
