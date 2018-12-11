/**
 * \file Test_HashAlign.h
 * \brief Unitary tests for class HashAlign
 */

#ifndef TEST_HASHALIGN_H
#define TEST_HASHALIGN_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../HashAlign.h"
#include "SDGBioSeqDB.h"

class Test_HashAlign: public CPPUNIT_NS::TestFixture
{

	CPPUNIT_TEST_SUITE(Test_HashAlign);
	CPPUNIT_TEST(test_search_on_loop);
 	CPPUNIT_TEST_SUITE_END();
	
	protected:
	  	void test_search_on_loop(void);
	private:
	  	SDGBioSeqDB initBioSeqDB(void);

};

#endif
