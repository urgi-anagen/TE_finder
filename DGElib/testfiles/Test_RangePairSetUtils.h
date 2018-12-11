/***
 *
 * Test_RangePairSetUtils.h
 *
 ***/

#ifndef TEST_RANGEPAIRSETUTILS_H
#define TEST_RANGEPAIRSETUTILS_H

#include <iostream>
#include <stdlib.h>
#include <SDGString.h>
#include <utility>
#include <map>
#include "../RangePairSet.h"

class Test_RangePairSetUtils
{
    public:

        static void viewRangePairSetList(std::list<RangePairSet> rpsList);
        static RangePairSet createRangePairSet1ForTest_inserted(void);
        static RangePairSet createRangePairSet2ForTest_inserted(void);
        static RangePairSet createRangePairSet1ForTest_inserted_true(void);
        static RangePairSet createRangePairSet2ForTest_inserted_true(void);
        static RangePairSet createRangePairSet1ForTest_overlapQ(void);
        static RangePairSet createRangePairSet2ForTest_overlapQ(void);
        static RangePairSet createRangePairSet1ForTest_overlapQ_length(void);
        static RangePairSet createRangePairSet2ForTest_overlapQ_length(void);
        static RangePairSet createRangePairSet1ForTest_diffQ(void);
        static RangePairSet createRangePairSet2ForTest_diffQ(void);
        static RangePairSet createExpRangePairSet1ForTest_diffQ(void);
        static RangePairSet createExpRangePairSet2ForTest_diffQ(void);
	static RangePairSet createRangePairSet1For_test_equality(void);
	static RangePairSet createRangePairSet2For_test_equality(void);
	static RangePairSet createRangePairSet1For_test_equality_not_equal(void);
	static RangePairSet createRangePairSet2For_test_equality_not_equal(void);
	static RangePairSet createRangePairSet1For_test_equality_on_match_part(void);
	static RangePairSet createRangePairSet2For_test_equality_on_match_part(void);
        static RangePairSet createRangePairSet1For_test_equality_rps1_rps2_different_length(void);
        static RangePairSet createRangePairSet2For_test_equality_rps1_rps2_different_length(void);
	static RangePairSet createInputRangePairSetFor_test_computeScoreWithLength_one_match_part(void);
	static RangePairSet createExpRangePairSetFor_test_computeScoreWithLength_one_match_part(void);
	static RangePairSet createInputRangePairSetFor_test_computeScoreWithLength_two_match_part(void);
	static RangePairSet createExpRangePairSetFor_test_computeScoreWithLength_two_match_part(void);
	static RangePairSet createRangePairSet1ForTest_overalpQ_length(void);
	static RangePairSet createRangePairSet2ForTest_overalpQ_length(void);
	static RangePairSet createRangePairSet1ForTest_mergeQ_score_rps1_greater_score_than_rps2(void);
	static RangePairSet createRangePairSet2ForTest_mergeQ_score_rps1_greater_score_than_rps2(void);
	static RangePairSet createExpRangePairSetForTest_mergeQ_score_rps1_greater_score_than_rps2(void);
	static RangePairSet createRangePairSet1ForTest_mergeQ_score_rps1_lower_score_than_rps2(void);
	static RangePairSet createRangePairSet2ForTest_mergeQ_score_rps1_lower_score_than_rps2(void);
	static RangePairSet createExpRangePairSetForTest_mergeQ_score_rps1_lower_score_than_rps2(void);
	static RangePairSet createRangePairSet1ForTest_overlapQ_rps2_include_in_rps1_but_no_overlap(void);
	static RangePairSet createRangePairSet2ForTest_overlapQ_rps2_include_in_rps1_but_no_overlap(void);
	static RangePairSet createRangePairSet1ForTest_overlapQ_length_2(void);
	static RangePairSet createRangePairSet2ForTest_overlapQ_length_2(void);
	static RangePairSet createRangePairSet1ForTest_mergeQ_rps2_overlap_on_rps1_frgt(void);
	static RangePairSet createRangePairSet2ForTest_mergeQ_rps2_overlap_on_rps1_frgt(void);
	static RangePairSet createRangePairSet1ForTest_diffQ_rps2_overlap_on_rps1_frgt_case1(void);
	static RangePairSet createRangePairSet2ForTest_diffQ_rps2_overlap_on_rps1_frgt_case1(void);
	static RangePairSet createRangePairSet1ForTest_diffQ_rps2_overlap_on_rps1_frgt_case2(void);
	static RangePairSet createRangePairSet2ForTest_diffQ_rps2_overlap_on_rps1_frgt_case2(void);
	static RangePairSet createRangePairSet1ForTest_overlap_first_and_second_fragment(void);
	static RangePairSet createRangePairSet2ForTest_overlap_first_and_second_fragment(void);
	
        static RangePairSet generateInputs_for_test_orientSubjects_on_plus_strand(void);
	static RangePairSet generateInputs_for_test_orientSubjects_on_minus_strand(void);
	static RangePairSet generateInputs_for_test_orientSubjects_path_size_is_1(void);
        static RangePairSet generateInputs_for_test_orientSubjects_orient_more_than_1_fragment(void);
        static RangePairSet generateInputs_for_test_orientSubjects_no_fragments_to_orient(void);

	static RangePairSet generateOutputs_for_test_orientSubjects_on_plus_strand(void);
	static RangePairSet generateOutputs_for_test_orientSubjects_on_minus_strand(void);
	static RangePairSet generateOutputs_for_test_orientSubjects_path_size_is_1(void);
        static RangePairSet generateOutputs_for_test_orientSubjects_orient_more_than_1_fragment(void);
        static RangePairSet generateOutputs_for_test_orientSubjects_no_fragments_to_orient(void);
	static RangePairSet createRangePairSet(std::list<SDGString> inputList);
};

#endif
