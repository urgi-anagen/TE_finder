#ifndef TEST_BLRMATCHMAP_H
#define TEST_BLRMATCHMAP_H

#include <stdlib.h>
#include <map>
#include <list>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "RangePair.h"
#include "BLRMatcherParameter.h"
#include "BLRMatchMap.h"
#include "Graph.h"

class Test_BLRMatchMap: public CPPUNIT_NS::TestFixture
{
  CPPUNIT_TEST_SUITE(Test_BLRMatchMap);
  CPPUNIT_TEST(test_mapAlign_equality);
  CPPUNIT_TEST(test_clean_conflicts);
  CPPUNIT_TEST(test_mapPath);
  CPPUNIT_TEST(view_add_clean_path_all_S);
  CPPUNIT_TEST(view_split_path);
  CPPUNIT_TEST(test_isOverlapFound_in_add_split_path);
  // CPPUNIT_TEST(test_mapPath_comparaison_diagnostic_compare_rangePairSetList_list_comparaison_key_1_2);
  // CPPUNIT_TEST(test_mapPath_comparaison_diagnostic_compare_rangePairSetList_list_comparaison_key_1_3);
  CPPUNIT_TEST(test_reComputeScoreWithLength_on_a_match_with_one_match_part);
  CPPUNIT_TEST(test_reComputeScoreWithLength_on_a_match_with_two_match_part);
  CPPUNIT_TEST(test_merge_on_two_queries);
  CPPUNIT_TEST(test_merge_second_fragment_include);
  CPPUNIT_TEST(test_merge_overlap_right_on_second_fragment);
  CPPUNIT_TEST(test_merge_overlap_left_on_second_fragment);
  CPPUNIT_TEST(test_merge_overlap_first_and_second_fragment);
  CPPUNIT_TEST(test_merge_all_included);
  CPPUNIT_TEST(test_merge_overlap_right_on_second_fragment);
  CPPUNIT_TEST(test_merge_overlap_left_on_second_fragment);
  CPPUNIT_TEST(test_merge_on_mapPath_data);
  CPPUNIT_TEST(test_clusterizeOverlapingRps);
  CPPUNIT_TEST(test_mergeOnCluster);
  CPPUNIT_TEST(test_writeBED);
  CPPUNIT_TEST_SUITE_END();
  public:
    void setUp();
    void tearDown();

	protected:
    
    void test_mapAlign_equality(void);
    void test_clean_conflicts(void);
    void test_mapPath(void);
    void view_add_clean_path_all_S(void);
    void view_split_path(void);
    
    void test_extractRangePairSetListFromMapPath(void);
    void test_isOverlapFound_in_add_split_path(void);
    void test_isOverlapFound_in_add_clean_path_all_S(void);
    // void test_mapPath_comparaison_diagnostic_compare_rangePairSetList_list_comparaison_key_1_2(void);
    // void test_mapPath_comparaison_diagnostic_compare_rangePairSetList_list_comparaison_key_1_3(void);
    void test_reComputeScoreWithLength_on_a_match_with_one_match_part(void);
    void test_reComputeScoreWithLength_on_a_match_with_two_match_part(void);
    void test_merge_on_two_queries(void);
    void test_merge_second_fragment_include(void);
    void test_merge_overlap_right_on_second_fragment(void);
    void test_merge_overlap_left_on_second_fragment(void);
    void test_merge_overlap_first_and_second_fragment(void);
    void test_merge_on_mapPath_data(void);
    void test_merge_all_included(void); 
    void test_clusterizeOverlapingRps(void);
    void test_mergeOnCluster(void);
    void test_writeBED(void);
};

#endif
