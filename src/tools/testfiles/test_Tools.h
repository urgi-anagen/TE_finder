#ifndef TEST_ORIENTER_H
#define TEST_ORIENTER_H

#include <cstdlib>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>


class test_Tools: public CPPUNIT_NS::TestFixture
{

	CPPUNIT_TEST_SUITE(test_Tools);
        CPPUNIT_TEST(initTestSeq);
        CPPUNIT_TEST(cutterDB);
        CPPUNIT_TEST(test_all_tools);
  	CPPUNIT_TEST_SUITE_END();

public:
    test_Tools(void) {}

    void setUp() {}

    void tearDown(){}

protected:
    void initTestSeq(void);
    void cutterDB(void);
    void test_a_tool(std::string tool_name, std::string parameters);
    void test_all_tools(void);
};
#endif
