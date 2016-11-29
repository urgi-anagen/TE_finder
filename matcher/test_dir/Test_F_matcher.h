#ifndef TEST_F_MATCHER_H
#define TEST_F_MATCHER_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

class Test_F_matcher: public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE(Test_F_matcher);
	CPPUNIT_TEST(test_runAsScript);
	CPPUNIT_TEST(test_runAsScript_small_data);
 	CPPUNIT_TEST_SUITE_END();

	public:
		void setUp();
		void tearDown();

	protected:
		void test_runAsScript(void);
		void test_runAsScript_small_data(void);
};

#endif
