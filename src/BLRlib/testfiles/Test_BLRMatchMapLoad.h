#ifndef TEST_BLRMATCHMAPLOAD_H
#define TEST_BLRMATCHMAPLOAD_H

#include <stdlib.h>
#include <map>
#include <list>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "RangePair.h"
#include "../../matcher/BLRMatcherParameter.h"
#include "BLRMatchMap.h"

class Test_BLRMatchMapLoad : public CPPUNIT_NS::TestFixture {
CPPUNIT_TEST_SUITE(Test_BLRMatchMapLoad);
        CPPUNIT_TEST(test_readAlign);
        CPPUNIT_TEST(test_load);
        CPPUNIT_TEST(test_readPath);
        CPPUNIT_TEST(test_loadPath);
    CPPUNIT_TEST_SUITE_END();

    BLRMatcherThreadsParameter createParameter(void) {
        BLRMatcherThreadsParameter para;
        para.setLenFilter(0);
        para.setEvalFilter(10);
        para.setIdFilter(0);
        return para;
    };

public:
    void setUp();

    void tearDown();

protected:
    void test_readAlign(void);

    void test_readPath(void);

    void test_load(void);

    void test_loadPath(void);


};
#endif
