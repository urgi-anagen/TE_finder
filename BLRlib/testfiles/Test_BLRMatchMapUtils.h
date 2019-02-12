/***
 *
 * Test_BLRMatchMapUtils.h
 *
 ***/

#ifndef TEST_BLRMATCHMAPUTILS_H
#define TEST_BLRMATCHMAPUTILS_H

#include <iostream>
#include <stdlib.h>
#include "SDGString.h"
#include <utility>
#include <map>
#include <sstream>
#include "BLRMatchMap.h"

class Test_BLRMatchMapUtils {
public:
	static void viewMapAlign(BLRMatchMap::MapAlign matchAlign);

	static void viewMapPath(BLRMatchMap::MapPath mapPath);

	static void viewMapPathWithLabel(BLRMatchMap::MapPath mapPath);

	static void viewRangePairSetList(std::list<RangePairSet>);

	static void viewNum2Name(std::map<long, std::string>);

	static void writeInputFile(void);

	static void writeInputFileWithTwoMatchesToBeJoined(void);

	static BLRMatcherThreadsParameter createParameter(void);

	static BLRMatcherThreadsParameter createParameterWithThreads(void);

	static BLRMatcherThreadsParameter createParameter(SDGString inputFileName);

	static BLRMatchMap::MapPath createExpMapPath_for_mapPath(void);

	static BLRMatchMap::MapAlign createExpMapAlign_for_clean_conflicts(void);

	static std::list<RangePairSet> createRpList_for_test_add_clean_path_all_S(void);

	static BLRMatchMap::MapAlign createMapAlign_instance_for_test_mapAlign_Equality(void);

	static std::list<RangePairSet> createRpListForTest_extractRangePairSetListFromMapPath(void);

	static BLRMatchMap::MapPath createMapPathForTest_extractRangePairSetListFromMapPath(void);

	static std::list<RangePairSet> createRpListForTest_add_split_path(void);

	static std::list<RangePairSet> createRpListForTest_isOverloapFound_in_add_split_path(void);

	static BLRMatchMap::MapPath createMapPathForTest_isOverlapFound_in_add_split_path(void);

	static BLRMatchMap::MapPath createMapPath_afterJoin(void);

	static BLRMatchMap::MapPath createExpMapPathForTest_reComputScoreWithLength(void);

	static BLRMatchMap::MapPath createExpMapPathForTest_mapPath_comparaison_diagnostic(void);

	static std::list<RangePairSet>
	createInputRpsList_for_test_reComputeScoreWithLength_on_a_match_with_one_match_part(void);

	static std::list<RangePairSet>
	createExpRpsList_for_test_reComputeScoreWithLength_on_a_match_with_one_match_part(void);

	static std::list<RangePairSet>
	createInputRpsList_for_test_reComputeScoreWithLength_on_a_match_with_two_match_part(void);

	static std::list<RangePairSet>
	createExpRpsList_for_test_reComputeScoreWithLength_on_a_match_with_two_match_part(void);

	static std::list<RangePairSet> createExpRpsListForOnlyJoinComputeScore(void);

	static bool areTwoRangePairSetEqualsWithoutIdentity(RangePairSet rps1, RangePairSet rps2);

	static std::list<RangePairSet> createInputRpsListFor_test_add_split_path(void);

	static BLRMatchMap::MapPath createMapPath(std::list<SDGString> inputList);

	static bool
	isOverlapFound_in_add_split_path(std::list<RangePairSet>::iterator iter, BLRMatchMap::MapPath mapPath, double idTolerance,
									 unsigned lenFilter);
};

#endif
