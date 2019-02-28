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

class Test_BLRGrouper : public CppUnit::TestFixture {

    BLRGrouperParameter *para_ptr;
    BLRMatchMap *match_map_ptr;


    CPPUNIT_TEST_SUITE(Test_BLRGrouper);

    CPPUNIT_TEST(test_mergeGroupsLists);

    CPPUNIT_TEST_SUITE_END();

    public:

        Test_BLRGrouper(void) : para_ptr(0), match_map_ptr(0){}

    void setUp()
    {
        para_ptr=new BLRGrouperParameter();
        match_map_ptr=new BLRMatchMap(para_ptr);
    }
    void tearDown()
    {
        delete match_map_ptr;
        delete para_ptr;
    }

    protected:

        void test_mergeGroupsLists(void);

    };

#endif //TE_FINDER_TEST_BLRGROUPER_H
