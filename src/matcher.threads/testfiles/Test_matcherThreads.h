#ifndef TEST_F_MATCHER_H
#define TEST_F_MATCHER_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../BLRMatcherThreads.h"
#include <iomanip>

class Test_matcherThreads: public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE(Test_matcherThreads);
	//CPPUNIT_TEST(test_runAsProcess);
	CPPUNIT_TEST(test_runAsScript_join_simple);
	CPPUNIT_TEST(test_runAsScript_join);
	CPPUNIT_TEST(test_runAsScript_join_threads);
	CPPUNIT_TEST(test_runAsScript_join_clean_threads);
	CPPUNIT_TEST(test_runAsScript_join_merge_threads);
	CPPUNIT_TEST(test_runAsScript_join_merge_clean_threads);
	CPPUNIT_TEST_SUITE_END();

	public:
		void setUp();
		void tearDown();

	protected:
		void test_runAsScript_join_simple(void);
        void test_runAsScript_join(void);
        void test_runAsScript_join_threads(void);
        void test_runAsScript_join_clean_threads(void);
        void test_runAsScript_join_merge_threads(void);
        void test_runAsScript_join_merge_clean_threads(void);
		//void test_runAsProcess(void);
        std::string path2exec;
};

#endif
