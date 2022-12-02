/**
 * \file Test_RangePairSet.h
 * \brief Unitary tests for class Test_RangePairSet.h
 */

#ifndef TEST_RANGEPAIRSET_H
#define TEST_RANGEPAIRSET_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "RangePairSet.h"

/**
 * \class Test_RangePairSet.h
 * \brief Unitary tests for the class Range
 */
class Test_RangePairSet: public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( Test_RangePairSet );
	CPPUNIT_TEST( test_inserted_false );
	CPPUNIT_TEST( test_inserted_true );
	CPPUNIT_TEST( test_overlapQ );
	CPPUNIT_TEST( test_overlapQ_rps2_include_in_rps1_but_no_overlap );
	CPPUNIT_TEST( test_overlapQ_length );
	CPPUNIT_TEST( test_equality );
	CPPUNIT_TEST( test_equality_not_equal );
	CPPUNIT_TEST( test_not_equality_on_match_part );
        CPPUNIT_TEST( test_overlapQ_length_length );
        CPPUNIT_TEST( test_overlapQ_length_length_symetric );
	CPPUNIT_TEST(view_diffQ_for_test_mergeQ);
	CPPUNIT_TEST(test_orientSubjects_on_plus_strand);
        CPPUNIT_TEST(test_orientSubjects_on_minus_strand);
        CPPUNIT_TEST(test_orientSubjects_path_size_is_1);
        CPPUNIT_TEST(test_orientSubjects_orient_more_than_1_fragment);
        CPPUNIT_TEST(test_orientSubjects_no_fragments_to_orient);
 	CPPUNIT_TEST_SUITE_END();
  
	public:
   void setUp();
	 void tearDown();

	protected:
   void test_inserted_false(void); 
   void test_inserted_true(void); 
   void test_overlapQ(void);
   void test_overlapQ_rps2_include_in_rps1_but_no_overlap(void);
   void test_overlapQ_length(void);
   void test_equality(void);
   void test_equality_not_equal(void);
   void test_not_equality_on_match_part(void);
   void test_overlapQ_length_length(void);
   void test_overlapQ_length_length_symetric(void);
   void view_diffQ_for_test_mergeQ(void);
   void test_orientSubjects_on_plus_strand(void);
   void test_orientSubjects_on_minus_strand(void);
   void test_orientSubjects_path_size_is_1(void);
   void test_orientSubjects_orient_more_than_1_fragment(void);
   void test_orientSubjects_no_fragments_to_orient(void);

};

#endif

