/**
 * \file Test_Range.h
 * \brief Unitary tests for class Range
 */

#ifndef TEST_DIAGCLONE_H
#define TEST_DIAGCLONE_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../DiagClone.h"

/**
 * \class Test_Range
 * \brief Unitary tests for the class Range
 */
class Test_DiagClone: public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( Test_DiagClone );
	CPPUNIT_TEST( test_add_push_back_case ); 	
	CPPUNIT_TEST( test_add_set_end_case ); 	
	CPPUNIT_TEST( test_add_set_start_case ); 	
	CPPUNIT_TEST_SUITE_END();
 
	protected:
	 void test_add_push_back_case( void );
	 void test_add_set_end_case( void );
	 void test_add_set_start_case( void );
};

#endif

