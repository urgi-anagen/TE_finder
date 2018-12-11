/**
 * \file Test_RangeAlignSet.h
 * \brief Unitary tests for class RangeAlignSet
 */

#ifndef TEST_RANGE_H
#define TEST_RANGE_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../RangeAlignSet.h"

/**
 * \class Test_RangeAlignSet
 * \brief Unitary tests for the class RangeAlignSet
 */
class Test_RangeAlignSet: public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( Test_RangeAlignSet );
	CPPUNIT_TEST( test_reset );
	CPPUNIT_TEST( test_operator_equal_true );
	CPPUNIT_TEST( test_operator_equal_false_coordinates );
	CPPUNIT_TEST( test_operator_equal_false_included );
	CPPUNIT_TEST( test_getStart_plusStrand );
	CPPUNIT_TEST( test_getStart_minusStrand );
	CPPUNIT_TEST( test_getEnd_plusStrand );
	CPPUNIT_TEST( test_getEnd_minusStrand );
	CPPUNIT_TEST( test_isPlusStrand_true );
	CPPUNIT_TEST( test_isPlusStrand_false );
	CPPUNIT_TEST( test_sortUsingStrand_plusStrand );
	CPPUNIT_TEST( test_sortUsingStrand_minusStrand );
	CPPUNIT_TEST( test_getIncluded );
	CPPUNIT_TEST( test_getLengthSet );
	CPPUNIT_TEST( test_reverse );
	CPPUNIT_TEST( test_overlap_length );
	CPPUNIT_TEST( test_isStrictlyContained_sameSize_true );
	CPPUNIT_TEST( test_isStrictlyContained_sameSize_false );
	CPPUNIT_TEST( test_isStrictlyContained_differentSizes_true_2 );
	CPPUNIT_TEST( test_isStrictlyContained_differentSizes_true_3 );
	CPPUNIT_TEST( test_isStrictlyContained_differentSizes_false );
	CPPUNIT_TEST( test_isContained );
	CPPUNIT_TEST( test_doesItContain );
	CPPUNIT_TEST( test_merge_set );
	CPPUNIT_TEST( test_merge_checkCoordinates );
	CPPUNIT_TEST( test_merge_checkIncluded_false_false );
	CPPUNIT_TEST( test_merge_checkIncluded_true_true );
	CPPUNIT_TEST( test_merge_checkIncluded_sameCoordinates_sameStrand_true_false );
	CPPUNIT_TEST( test_merge_checkIncluded_1false_2true_2notContained );
	CPPUNIT_TEST( test_merge_checkIncluded_1true_2false_2contained );
	CPPUNIT_TEST( test_merge_checkIncluded_1true_2false_2notContained );
	CPPUNIT_TEST( test_merge_checkIncluded_1false_2true_1contained );
	CPPUNIT_TEST( test_hasSameRanges_true );
	CPPUNIT_TEST( test_hasSameRanges_true_unsorted );
	CPPUNIT_TEST( test_hasSameRanges_true_diffIncluded );
	CPPUNIT_TEST( test_hasSameRanges_falseCoords );
 	CPPUNIT_TEST_SUITE_END();

	public:
		void setUp();
		void tearDown();

	protected:
		void test_reset( void );
		void test_operator_equal_true( void );
		void test_operator_equal_false_coordinates( void );
		void test_operator_equal_false_included( void );
		void test_getStart_plusStrand( void );
		void test_getStart_minusStrand( void );
		void test_getEnd_plusStrand( void );
		void test_getEnd_minusStrand( void );
		void test_isPlusStrand_true( void );
		void test_isPlusStrand_false( void );
		void test_sortUsingStrand_plusStrand( void );
		void test_sortUsingStrand_minusStrand( void );
		void test_getIncluded( void );
		void test_getLengthSet( void );
		void test_reverse( void );
		void test_overlap_length( void );
		void test_isStrictlyContained_sameSize_true( void );
		void test_isStrictlyContained_sameSize_false( void );
		void test_isStrictlyContained_differentSizes_true_2( void );
		void test_isStrictlyContained_differentSizes_true_3( void );
		void test_isStrictlyContained_differentSizes_false( void );
		void test_isContained( void );
		void test_doesItContain( void );
		void test_merge_set( void );
		void test_merge_checkCoordinates( void );
		void test_merge_checkIncluded_false_false( void );
		void test_merge_checkIncluded_true_true( void );
		void test_merge_checkIncluded_sameCoordinates_sameStrand_true_false( void );
		void test_merge_checkIncluded_1false_2true_2notContained( void );
		void test_merge_checkIncluded_1true_2false_2contained( void );
		void test_merge_checkIncluded_1true_2false_2notContained( void );
		void test_merge_checkIncluded_1false_2true_1contained( void );
		void test_hasSameRanges_true( void );
		void test_hasSameRanges_true_unsorted( void );
		void test_hasSameRanges_true_diffIncluded( void );
		void test_hasSameRanges_falseCoords( void );

};

#endif
