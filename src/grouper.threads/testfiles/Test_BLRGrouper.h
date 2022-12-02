//
// Created by Hadi Quesneville on 2019-02-14.
//

#ifndef TE_FINDER_TEST_BLRGROUPER_H
#define TE_FINDER_TEST_BLRGROUPER_H
#include <cppunit/extensions/HelperMacros.h>
#include "Range.h"
#include "RangeAlign.h"
#include "RangeAlignSet.h"
#include "../BLRGrouperThreads.h"
#include <list>
#include <vector>
#include <iomanip>

class Test_BLRGrouper : public CppUnit::TestFixture {

    BLRGrouperParameter para;
    BLRMatchMap match_map;


    CPPUNIT_TEST_SUITE(Test_BLRGrouper);

    CPPUNIT_TEST(test_mergeGroupsLists);

    CPPUNIT_TEST_SUITE_END();

    public:

        Test_BLRGrouper(void) {}

    void setUp()
    {
    }
    void tearDown()
    {
    }

    protected:

        void test_mergeGroupsLists(void);

    };

#endif //TE_FINDER_TEST_BLRGROUPER_H
