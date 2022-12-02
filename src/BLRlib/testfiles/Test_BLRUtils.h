//
// Created by Hadi Quesneville on 25/03/2020.
//

#ifndef TE_FINDER_TEST_BLRUTILS_H
#define TE_FINDER_TEST_BLRUTILS_H
#include "BLRJoinParameter.h"
#include "BLRMatchMap.h"
#include "RangePairSet.h"

class Test_BLRUtils {
public:

    static void writeInputFile(void);

    static BLRJoinParameter createParameter(void);

    static BLRMatchMap::MapAlign createExpMapAlign_for_clean_conflicts(void);

    static std::list<RangePairSet> createRpList_for_test_add_clean_path_all_S(void);

    static BLRMatchMap::MapAlign createMapAlign_instance_for_test_mapAlign_Equality(void);

    static std::list<RangePairSet> createRpListForTest_isOverloapFound_in_add_split_path(void);

    static BLRMatchMap::MapPath createMapPathForTest_isOverlapFound_in_add_split_path(void);

    static std::list<RangePairSet>
    createInputRpsList_for_test_reComputeScoreWithLength_on_a_match_with_one_match_part(void);

    static std::list<RangePairSet>
    createExpRpsList_for_test_reComputeScoreWithLength_on_a_match_with_one_match_part(void);

    static std::list<RangePairSet>
    createInputRpsList_for_test_reComputeScoreWithLength_on_a_match_with_two_match_part(void);

    static std::list<RangePairSet>
    createExpRpsList_for_test_reComputeScoreWithLength_on_a_match_with_two_match_part(void);

    static bool
    isOverlapFound_in_add_split_path(std::list<RangePairSet>::iterator iter, BLRMatchMap::MapPath mapPath,
                                     double idTolerance,
                                     unsigned lenFilter);
};


#endif //TE_FINDER_TEST_BLRUTILS_H
