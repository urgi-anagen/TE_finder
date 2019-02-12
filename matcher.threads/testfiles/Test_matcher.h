#ifndef TEST_F_MATCHER_H
#define TEST_F_MATCHER_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../BLRMatcherThreads.h"

class Test_matcherThreads: public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE(Test_matcherThreads);
	CPPUNIT_TEST(test_runAsProcess);
	CPPUNIT_TEST(test_runAsScript);
 	CPPUNIT_TEST_SUITE_END();

	public:
		void setUp();
		void tearDown();

	protected:
		void test_runAsScript(void);
		void test_runAsProcess(void);
};

#endif
