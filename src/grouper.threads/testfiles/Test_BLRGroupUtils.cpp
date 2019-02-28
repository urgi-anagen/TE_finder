#include <stdlib.h>
#include "Test_BLRGroupUtils.h"
#include "../BLRGrouperParameter.h"
#include "FileUtils.h"
#include "Range.h"
#include "RangeAlignSet.h"
#include "Test_RangePairSetUtils.h"

BLRGrouperParameter Test_BLRGroupUtils::createParameter(void)
{
	
    	BLRGrouperParameter para;
	para.setJoin_frag(true);
	para.setSizeFiler(3);	
	return para;
}

void Test_BLRGroupUtils::write_input_file_for_test_group(SDGString inputFileName)
{
	std::ofstream inputFileStream;
	FileUtils::openFile(inputFileName, inputFileStream);
	inputFileStream<<"dmel_chr4\t434924\t434969\tdmel_chr4\t437318\t437273\t1e-12\t152\t95.65\n";
	inputFileStream<<"dmel_chr4\t434924\t434969\tdmel_chr4\t437318\t437273\t1e-12\t152\t95.65\n";
	inputFileStream<<"dmel_chr4\t435717\t435762\tdmel_chr4\t437318\t437273\t1e-12\t152\t95.65\n";
	inputFileStream<<"dmel_chr4\t435717\t435762\tdmel_chr4\t437318\t437273\t1e-12\t152\t95.65\n";
	inputFileStream<<"dmel_chr4\t74585\t74628\tdmel_chr4\t952036\t951992\t2e-11\t72\t95.56\n";
	inputFileStream<<"dmel_chr4\t74585\t74630\tdmel_chr4\t802413\t802367\t1e-12\t76\t95.74\n";
	inputFileStream<<"dmel_chr4\t74585\t74630\tdmel_chr4\t802413\t802367\t1e-12\t76\t95.74\n";
	inputFileStream<<"dmel_chr4\t879407\t879454\tdmel_chr4\t1170969\t1170922\t2e-11	72\t93.75\n";
	inputFileStream<<"dmel_chr4\t953550\t953597\tdmel_chr4\t1170922\t1170969\t2e-11	72\t93.88\n";
	inputFileStream.close();
}

void Test_BLRGroupUtils::write_map_file_for_test_group(SDGString mapFileName)
{
	std::ofstream mapFileStream;
	FileUtils::openFile(mapFileName, mapFileStream);
	mapFileStream<<"MbQ1Gr1Cl0\tdmel_chr4\t434924\t434969\n";
	mapFileStream<<"MbS2Gr1Cl0\tdmel_chr4\t437318\t437273\n";
	mapFileStream<<"MbQ3Gr1Cl0\tdmel_chr4\t435717\t435762\n";
	mapFileStream<<"MbS4Gr2Cl0\tdmel_chr4\t952036\t951992\n";
	mapFileStream<<"MbS5Gr2Cl0\tdmel_chr4\t802413\t802367\n";
	mapFileStream<<"MbQ6Gr2Cl0\tdmel_chr4\t74585\t74630\n";
	mapFileStream<<"MbQ7Gr3Cl0\tdmel_chr4\t953597\t953550\n";
	mapFileStream<<"MbS8Gr3Cl0\tdmel_chr4\t1170969\t1170922\n";
	mapFileStream<<"MbQ9Gr3Cl0\tdmel_chr4\t879407\t879454\n";
	mapFileStream.close();
}

void Test_BLRGroupUtils::write_input_file_for_test_group_load_path(SDGString pathFileName)
{
	std::ofstream pathFileStream;
	FileUtils::openFile(pathFileName, pathFileStream);
	pathFileStream<<"10\tchunk682\t46799\t46839\trefTE_230\t2977\t3022\t2.2e-23\t41\t95\n";
	pathFileStream<<"11\tchunk682\t54450\t54511\trefTE_230\t3046\t2977\t0\t62\t93.4426\n";
	pathFileStream<<"12\tchunk682\t54604\t54879\trefTE_230\t2353\t2639\t0\t276\t87.8327\n";
	pathFileStream<<"13\tchunk682\t55604\t55920\trefTE_230\t2743\t3047\t0\t317\t91.3333\n";
	pathFileStream<<"14\tchunk682\t62251\t62526\trefTE_230\t2353\t2639\t0\t276\t84.3866\n";
	pathFileStream<<"14\tchunk682\t62662\t62891\trefTE_230\t2613\t2844\t0\t230\t90.75\n";
	pathFileStream<<"15\tchunk682\t62608\t62661\trefTE_230\t2538\t2595\t0\t54\t90.942\n";
	pathFileStream.close();
}

void Test_BLRGroupUtils::write_map_file_for_test_group_load_path(SDGString mapFileName)
{
	std::ofstream mapFileStream;
	FileUtils::openFile(mapFileName, mapFileStream);
	mapFileStream<<"MbQ1Gr1Cl0\tchunk682\t46839\t46799\n";
	mapFileStream<<"MbS2Gr1Cl0\tchunk682\t3046\t2977\n";
	mapFileStream<<"MbQ3Gr1Cl0\tchunk682\t54450\t54511\n";
	mapFileStream<<"MbQ4Gr2Cl0\tchunk682\t54604\t54879\n";
	mapFileStream<<"MbQ5Gr2Cl0\tchunk682\t55604\t55920\n";
	mapFileStream<<"MbS6Gr2Cl0\tchunk682\t2353\t3047\n";
	mapFileStream<<"MbQ7Gr2Cl0\tchunk682\t62251\t62891\n";
	mapFileStream.close();
}


void Test_BLRGroupUtils::write_input_file_for_test_compare_join_for_debug(SDGString inputFileName)
{
	std::ofstream inputFileStream;
	FileUtils::openFile(inputFileName, inputFileStream);
	inputFileStream<<"dmel_chr4\t13293\t13592\tdmel_chr4\t74091\t74394\t8e-75\t281\t87.18\n";
	inputFileStream<<"dmel_chr4\t13590\t13726\tdmel_chr4\t74436\t74572\t1e-52\t208\t94.2\n";
	inputFileStream<<"dmel_chr4\t13718\t13762\tdmel_chr4\t1260453\t1260409\t4e-12\t74\t95.56\n";
	inputFileStream<<"dmel_chr4\t13609\t13726\tdmel_chr4\t1260568\t1260459\t2e-17\t92\t85.59\n";
	inputFileStream<<"dmel_chr4\t13305\t13592\tdmel_chr4\t1197346\t11976366\t8e-72\t272\t87.16\n";
	inputFileStream<<"dmel_chr4\t13590\t13746\tdmel_chr4\t1197677\t1197834\t4e-55\t216\t92.41\n";
	inputFileStream.close();
}

void Test_BLRGroupUtils::write_align_file_for_test_group_compare_intermediate_rpsList_between_2_loads (SDGString inputFileName)
{
	std::ofstream inputFileStream;
	FileUtils::openFile(inputFileName, inputFileStream);
	inputFileStream<<"dmel_chr4\t13293\t13592\tdmel_chr4\t74091\t74394\t8e-75\t281\t87.18\n";
	inputFileStream<<"dmel_chr4\t13590\t13726\tdmel_chr4\t74436\t74572\t1e-52\t208\t94.2\n";
	inputFileStream<<"dmel_chr4\t13718\t13762\tdmel_chr4\t1260453\t1260409\t4e-12\t74\t95.56\n";
	inputFileStream<<"dmel_chr4\t13609\t13726\tdmel_chr4\t1260568\t1260459\t2e-17\t92\t85.59\n";
	inputFileStream<<"dmel_chr4\t13305\t13592\tdmel_chr4\t1197346\t1197636\t8e-72\t272\t87.16\n";
	inputFileStream<<"dmel_chr4\t13590\t13746\tdmel_chr4\t1197677\t1197834\t4e-55\t216\t92.41\n";
	
	inputFileStream.close();
}
void Test_BLRGroupUtils::write_path_file_for_test_group_compare_intermediate_rpsList_between_2_loads (SDGString inputFileName)
{
	std::ofstream inputFileStream;
	FileUtils::openFile(inputFileName, inputFileStream);
	inputFileStream<<"1\tdmel_chr4\t13293\t13592\tdmel_chr4\t74091\t74394\t8e-75\t281\t87.18\n";
	inputFileStream<<"1\tdmel_chr4\t13590\t13726\tdmel_chr4\t74436\t74572\t1e-52\t208\t94.2\n";
	inputFileStream<<"2\tdmel_chr4\t13718\t13762\tdmel_chr4\t1260453\t1260409\t4e-12\t74\t95.56\n";
	inputFileStream<<"2\tdmel_chr4\t13609\t13726\tdmel_chr4\t1260568\t1260459\t2e-17\t92\t85.59\n";
	inputFileStream<<"3\tdmel_chr4\t13305\t13592\tdmel_chr4\t1197346\t1197636\t8e-72\t272\t87.16\n";
	inputFileStream<<"4\tdmel_chr4\t13590\t13746\tdmel_chr4\t1197677\t1197834\t4e-55\t216\t92.41\n";
	inputFileStream.close();

}
