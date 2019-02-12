/***
 *
 * Test_BLRGroupUtils.h
 *
 ***/

#ifndef TEST_BLRGROUPUTILS_H
#define TEST_BLRGROUPUTILS_H

#include <iostream>
#include <stdlib.h>
#include <SDGString.h>
#include <utility>
#include <map>
#include "../BLRGrouperParameter.h"
#include "../BLRMemIdxBin.h"

class Test_BLRGroupUtils
{
	public:

		static BLRGrouperParameter createParameter(void); 
		static void write_input_file_for_test_group(SDGString inputFileName);
		static void write_map_file_for_test_group(SDGString mapFileName);
		static void write_input_file_for_test_group_load_path(SDGString pathFileName);
		static void write_map_file_for_test_group_load_path(SDGString mapFileName);
		static void write_input_file_for_test_compare_join_for_debug(SDGString inputFileName);
		static void write_align_file_for_test_group_compare_intermediate_rpsList_between_2_loads (SDGString inputFileName);
		static void write_path_file_for_test_group_compare_intermediate_rpsList_between_2_loads (SDGString inputFileName);

};

#endif
