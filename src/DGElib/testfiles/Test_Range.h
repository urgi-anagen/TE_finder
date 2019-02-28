/**
 * \file Test_Range.h
 * \brief Unitary tests for class Range
 */

#ifndef TEST_RANGE_H
#define TEST_RANGE_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "Range.h"

/**
 * \class Test_Range
 * \brief Unitary tests for the class Range
 */
class Test_Range: public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( Test_Range );
	CPPUNIT_TEST( test_getStart );
	CPPUNIT_TEST( test_getEnd );
	CPPUNIT_TEST( test_isPlusStrand );
	CPPUNIT_TEST( test_getLength );
	CPPUNIT_TEST( test_getMin );
	CPPUNIT_TEST( test_getMax );
	CPPUNIT_TEST( test_reverse );
	CPPUNIT_TEST( test_overlap );
	CPPUNIT_TEST( test_isContained );
	CPPUNIT_TEST(test_diff_range1_start_changed_after_diff);
 	CPPUNIT_TEST_SUITE_END();
  
	public:
     void setUp();
	 void tearDown();

	protected:
	 void test_getStart( void );
	 void test_getEnd( void );
	 void test_isPlusStrand( void );
	 void test_getLength( void );
	 void test_getMin( void );
	 void test_getMax( void );
	 void test_reverse( void );
	 void test_overlap( void );
	 void test_isContained(void);
	 void test_diff_range1_start_changed_after_diff(void);

};

#endif

