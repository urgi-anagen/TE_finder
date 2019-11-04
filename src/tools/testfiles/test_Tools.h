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
  	CPPUNIT_TEST_SUITE_END();

public:
    test_Tools(void) {}

    void setUp() {}

    void tearDown(){}

protected:
    void initTestSeq(void);
    void cutterDB(void);

};
#endif
