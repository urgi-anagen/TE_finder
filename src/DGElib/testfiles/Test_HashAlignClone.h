/**
 * \file Test_HashAlign.h
 * \brief Unitary tests for class HashAlign
 */

#ifndef TEST_HASHALIGNCLONE_H
#define TEST_HASHALIGNCLONE_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "HashAlignClone.h"
#include "SDGBioSeqDB.h"

class Test_HashAlignClone: public CPPUNIT_NS::TestFixture
{

	CPPUNIT_TEST_SUITE(Test_HashAlignClone);
	CPPUNIT_TEST(test_add_push_back_case);
	CPPUNIT_TEST(test_add_set_end_case);
	CPPUNIT_TEST(test_add_set_start_case);
 	CPPUNIT_TEST_SUITE_END();
	
	protected:
	  	void test_add_push_back_case(void);
	  	void test_add_set_start_case(void);
	  	void test_add_set_end_case(void);
};

#endif
