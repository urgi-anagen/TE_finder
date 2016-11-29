#ifndef TEST_BLRMATCHMAPLOADUTILS_H
#define TEST_BLRMATCHMAPLOADUTILS_H

#include <iostream>
#include <stdlib.h>
#include "SDGString.h"
#include <utility>
#include <map>
#include "BLRMatchMap.h"

class Test_BLRMatchMapLoadUtils
{
	public:
		static BLRMatchMap::MapAlign createExpMapAlign_for_test_load(void);
		static void write_input_file_for_test_loadPath_name2Num_and_num2Name(SDGString fileName);
		static std::map<std::string,long> generateName2NumQ(void);
		static std::map<long, std::string> generateNum2NameQ(void);
		static std::map<std::string,long> generateName2NumS(void);
		static std::map<long, std::string> generateNum2NameS(void);
		static void viewName2Num(std::map<std::string,long> mapToView);
		static void viewNum2Name(std::map<long,std::string> mapToView);
		static BLRMatcherParameter createParameter(void);
		static std::list<RangePairSet> createRpsList(std::list<SDGString> inputList);
		static std::list<SDGString> generateExp_for_test_loadPath_rpsList(void);
		static void write_input_file_for_test_loadPath_rpsList(SDGString fileName);
		static void write_input_file_for_test_loadPath_with_join_end_of_file(SDGString fileName);
		static std::list<SDGString> generateExp_for_test_loadPath_with_join_end_of_file();
		static void write_input_file_for_test_loadPath_with_join_begin_of_file(SDGString fileName);
		static std::list<SDGString> generateExp_for_test_loadPath_with_join_begin_of_file();
		static void write_input_file_for_test_loadPath_with_join_middle_of_file(SDGString fileName);
		static std::list<SDGString> generateExp_for_test_loadPath_with_join_middle_of_file();


};
#endif
