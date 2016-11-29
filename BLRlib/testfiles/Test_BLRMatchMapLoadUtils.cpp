#include <stack>
#include "Test_BLRMatchMapLoadUtils.h"
#include "Test_BLRMatchMapUtils.h"
#include "SDGString.h"
#include "FileUtils.h"
#include "BLRMatcherParameter.h"
#include "BLRMatchMap.h"
#include "Range.h"
#include "RangePair.h"
#include "RangePairSet.h"

void Test_BLRMatchMapLoadUtils::viewName2Num(std::map<std::string,long> mapToView)
{
	for (std::map<std::string,long>::iterator it=mapToView.begin(); it!=mapToView.end(); ++it)
    		std::cout << it->first << " => " << it->second << '\n';
}

void Test_BLRMatchMapLoadUtils::viewNum2Name(std::map<long,std::string> mapToView)
{
	for (std::map<long,std::string>::iterator it=mapToView.begin(); it!=mapToView.end(); ++it)
    		std::cout << it->first << " => " << it->second << '\n';
}

BLRMatcherParameter Test_BLRMatchMapLoadUtils::createParameter(void){
  BLRMatcherParameter para;
  para.setLenFilter(0);
  para.setEvalFilter(10);
  para.setIdFilter(0);
  return para;
}

std::list<RangePairSet> Test_BLRMatchMapLoadUtils::createRpsList(std::list<SDGString> inputList)
{
	std::list<RangePairSet> rpsContainer;	
	for(std::list<SDGString>::iterator input_it=inputList.begin(); input_it!=inputList.end();input_it++){
	
		size_t posFirstTab = input_it->find_first_of("\t");
		size_t posLastTab = input_it->find_last_of("\t");
		size_t posEndLine = input_it->find_first_of("\\n");
		size_t lengthSeq = posEndLine - posLastTab;
		// TODO + 1 is inconsistent with other utils ...
		size_t lengthInit = (posLastTab - posFirstTab) + 1;
		
		// isARangePairSetToCreate		
		if (posFirstTab != 0){
			// createARangePairSet 
			RangePairSet currentRps = RangePairSet(input_it->substr(0, lengthInit));
			SDGString sNameSeq = currentRps.getRangeS().getNameSeq();
			SDGString qNameSeq = currentRps.getRangeQ().getNameSeq();
			currentRps.getRangeS().setNumChr(atol(sNameSeq));
			currentRps.getRangeQ().setNumChr(atol(qNameSeq));
			currentRps.getRangeS().setNameSeq("");
			currentRps.getRangeQ().setNameSeq("");
			 
	
			currentRps.setLength(atol(input_it->substr(posLastTab, lengthSeq)));
			// empty path
			std::list<RangePair> path;
			currentRps.setPathDirectly(path);

			// storeRangePairSetCreatedInContainer
			rpsContainer.push_back(currentRps);

		// isARangePairToCreate	
		}else{	
			// createARangePair
			RangePair rp = RangePair(input_it->substr(posFirstTab, lengthInit));
			SDGString sNameSeq = rp.getRangeS().getNameSeq();
			SDGString qNameSeq = rp.getRangeQ().getNameSeq();
			rp.getRangeS().setNumChr(atol(sNameSeq));
			rp.getRangeQ().setNumChr(atol(qNameSeq));
			rp.setLength(atol(input_it->substr(posLastTab, lengthSeq)));
			
			rp.getRangeS().setNameSeq("");
			rp.getRangeQ().setNameSeq("");
			// getRangePairSetFromContainer
			RangePairSet currentRps = rpsContainer.back();			
			rpsContainer.pop_back();		
			
			// updatePath
			std::list<RangePair> path = currentRps.getPath();
			std::list<RangePair> newPath;
			for(std::list<RangePair>::iterator lrp_it=path.begin(); lrp_it!=path.end();lrp_it++){
				RangePair newRp = *lrp_it;	
            			newPath.push_back(newRp); 
    			}
			newPath.push_back(rp);
			currentRps.setPathDirectly(newPath);

			// storeRangePairSetCreatedInContainer
			rpsContainer.push_back(currentRps);
		}
	
	}
	
	return rpsContainer;
}

BLRMatchMap::MapAlign Test_BLRMatchMapLoadUtils::createExpMapAlign_for_test_load(void)
{
  SDGString match_file = "match.align";
  Test_BLRMatchMapUtils::writeInputFile();
  BLRMatcherParameter para = Test_BLRMatchMapUtils::createParameter();
  BLRMatchMap matchMap(&para);
  matchMap.clear();

  SDGString line1 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair1 = RangePair(line1);
  rangePair1.getRangeQ().setNumChr(1);
  rangePair1.getRangeQ().setNameSeq("");
  rangePair1.getRangeS().setNumChr(1);
  rangePair1.getRangeS().setNameSeq("");

  matchMap.insert(rangePair1);
  
  SDGString line2 = " \t109601\t109708\t \t3239\t3338\t3e-28\t108\t81.8182";
  RangePair rangePair2 = RangePair(line2);
  rangePair2.getRangeQ().setNumChr(1);
  rangePair2.getRangeQ().setNameSeq("");
  rangePair2.getRangeS().setNumChr(1);
  rangePair2.getRangeS().setNameSeq("");

  matchMap.insert(rangePair2);
    
  SDGString line3 = " \t109527\t109818\t \t3005\t3316\t0\t292\t82.8467"; 
  RangePair rangePair3 = RangePair(line3);
  rangePair3.getRangeQ().setNumChr(1);
  rangePair3.getRangeQ().setNameSeq("");
  rangePair3.getRangeS().setNumChr(1);
  rangePair3.getRangeS().setNameSeq("");

  matchMap.insert(rangePair3);
  
  SDGString line4 = " \t105091\t105206\t \t3351\t3238\t0\t116\t76.1062"; 
  RangePair rangePair4 = RangePair(line4);
  rangePair4.getRangeQ().setNumChr(1);
  rangePair4.getRangeQ().setNameSeq("");
  rangePair4.getRangeS().setNumChr(1);
  rangePair4.getRangeS().setNameSeq("");
 
  matchMap.insert(rangePair4);

  SDGString line5 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361"; 
  RangePair rangePair5 = RangePair(line5);
  rangePair5.getRangeQ().setNumChr(1);
  rangePair5.getRangeQ().setNameSeq("");
  rangePair5.getRangeS().setNumChr(1);
  rangePair5.getRangeS().setNameSeq("");
 
  matchMap.insert(rangePair5);

  SDGString line6 = " \t109951\t110569\t \t119\t689\t9.4e-19\t619\t77.768"; 
  RangePair rangePair6 = RangePair(line6);
  rangePair6.getRangeQ().setNumChr(1);
  rangePair6.getRangeQ().setNameSeq("");
  rangePair6.getRangeS().setNumChr(2);
  rangePair6.getRangeS().setNameSeq("");
 
  matchMap.insert(rangePair6);

  SDGString line7 = " \t110567\t110878\t \t2833\t2532\t0\t312\t77.2076"; 
  RangePair rangePair7 = RangePair(line7);
  rangePair7.getRangeQ().setNumChr(1);
  rangePair7.getRangeQ().setNameSeq("");
  rangePair7.getRangeS().setNumChr(3);
  rangePair7.getRangeS().setNameSeq("");
 
  matchMap.insert(rangePair7);
  return matchMap.getMapAlign();
}

void Test_BLRMatchMapLoadUtils::write_input_file_for_test_loadPath_name2Num_and_num2Name(SDGString fileName)
{
	std::ofstream inputFileStream;
	FileUtils::openFile(fileName, inputFileStream);
	inputFileStream<<"10\tchunk682\t46799\t46839\trefTE_230\t2977\t3022\t2.2e-23\t41\t95\n";
	inputFileStream<<"11\tchunk682\t54450\t54511\trefTE_230\t3046\t2977\t0\t62\t93.4426\n";
	inputFileStream<<"12\tchunk682\t54604\t54879\trefTE_230\t2353\t2639\t0\t276\t87.8327\n";
	inputFileStream<<"13\tchunk682\t55604\t55920\trefTE_230\t2743\t3047\t0\t317\t91.3333\n";
	inputFileStream<<"14\tchunk682\t62251\t62526\trefTE_230\t2353\t2639\t0\t276\t84.3866\n";
	inputFileStream.close();
}

std::map<std::string,long> Test_BLRMatchMapLoadUtils::generateName2NumQ(void)
{
	std::map<std::string,long> mapQ;
	mapQ["chunk682"] = 1;
	return mapQ;
}

std::map<long, std::string> Test_BLRMatchMapLoadUtils::generateNum2NameQ(void)
{
	std::map<long,std::string> mapQ;
	mapQ[1] = "chunk682";
	return mapQ;
}

std::map<std::string,long> Test_BLRMatchMapLoadUtils::generateName2NumS(void)
{
	std::map<std::string,long> mapS;
	mapS["refTE_230"] = 1;
	return mapS;
}

std::map<long, std::string> Test_BLRMatchMapLoadUtils::generateNum2NameS(void)
{
	std::map<long,std::string> mapS;
	mapS[1] = "refTE_230";
	return mapS;
}

void Test_BLRMatchMapLoadUtils::write_input_file_for_test_loadPath_rpsList(SDGString fileName)
{
	Test_BLRMatchMapLoadUtils::write_input_file_for_test_loadPath_name2Num_and_num2Name(fileName);
}

std::list<SDGString> Test_BLRMatchMapLoadUtils::generateExp_for_test_loadPath_rpsList(void)
{
	SDGString str1 = "1\t46799\t46839\t1\t2977\t3022\t2.2e-23\t41\t95\t45\n";
	SDGString str11 = "\t1\t46799\t46839\t1\t2977\t3022\t2.2e-23\t41\t95\t45\n";
	
	SDGString str2 = "1\t54450\t54511\t1\t3046\t2977\t0\t62\t93.4426\t69\n";
	SDGString str21 = "\t1\t54450\t54511\t1\t3046\t2977\t0\t62\t93.4426\t69\n";

	SDGString str3 = "1\t54604\t54879\t1\t2353\t2639\t0\t276\t87.8327\t317\n";
	SDGString str31 = "\t1\t54604\t54879\t1\t2353\t2639\t0\t276\t87.8327\t317\n";

	SDGString str4 = "1\t55604\t55920\t1\t2743\t3047\t0\t317\t91.3333\t286\n";
	SDGString str41 = "\t1\t55604\t55920\t1\t2743\t3047\t0\t317\t91.3333\t286\n";
	
	SDGString str5 = "1\t62251\t62526\t1\t2353\t2639\t0\t276\t84.3866\t45\n";
	SDGString str51 = "\t1\t62251\t62526\t1\t2353\t2639\t0\t276\t84.3866\t45\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1); 
	lPath.push_back(str11); 
	lPath.push_back(str2); 
	lPath.push_back(str21);
	lPath.push_back(str3); 
	lPath.push_back(str31); 
	lPath.push_back(str4);
	lPath.push_back(str41);
	lPath.push_back(str5);
	lPath.push_back(str51);

	return lPath;
}

void Test_BLRMatchMapLoadUtils::write_input_file_for_test_loadPath_with_join_end_of_file(SDGString fileName)
{
	std::ofstream inputFileStream;
	FileUtils::openFile(fileName, inputFileStream);
	//inputFileStream<<"10\tchunk682\t46799\t46839\trefTE_230\t2977\t3022\t2.2e-23\t41\t95\n";
	//inputFileStream<<"11\tchunk682\t54450\t54511\trefTE_230\t3046\t2977\t0\t62\t93.4426\n";
	//inputFileStream<<"12\tchunk682\t54604\t54879\trefTE_230\t2353\t2639\t0\t276\t87.8327\n";
	//inputFileStream<<"13\tchunk682\t55604\t55920\trefTE_230\t2743\t3047\t0\t317\t91.3333\n";
	//inputFileStream<<"14\tchunk682\t62251\t62526\trefTE_230\t2353\t2639\t0\t276\t84.3866\n";
	//inputFileStream<<"14\tchunk682\t62662\t62891\trefTE_230\t2613\t2844\t0\t230\t90.75\n";

	inputFileStream<<"15\tchunk682\t62608\t62661\trefTE_230\t2538\t2595\t0\t54\t90.942\n";
	inputFileStream<<"16\tchunk682\t67771\t67929\trefTE_230\t2520\t2353\t0\t159\t87.7419\n";
	inputFileStream<<"16\tchunk682\t67212\t67460\trefTE_230\t3046\t2761\t0\t249\t88.9831\n";



	inputFileStream.close();
}


std::list<SDGString> Test_BLRMatchMapLoadUtils::generateExp_for_test_loadPath_with_join_end_of_file()
{
	SDGString str1 = "1\t62608\t62661\t1\t2538\t2595\t0\t54\t90.942\t57\n";
	SDGString str11 = "\t1\t62608\t62661\t1\t2538\t2595\t0\t54\t90.942\t57\n";
	SDGString str2 = "1\t67212\t67929\t1\t3046\t2353\t0\t338\t88.5245\t452\n";
	SDGString str21= "\t1\t67771\t67929\t1\t2520\t2353\t0\t159\t87.7419\t167\n";
	SDGString str22 ="\t1\t67212\t67460\t1\t3046\t2761\t0\t249\t88.9831\t285\n";

	
	std::list<SDGString> lPath;
	lPath.push_back(str1);
	lPath.push_back(str11);
	lPath.push_back(str2);
	lPath.push_back(str21);
	lPath.push_back(str22);
	return lPath;
}

void Test_BLRMatchMapLoadUtils::write_input_file_for_test_loadPath_with_join_begin_of_file(SDGString fileName)
{
	std::ofstream inputFileStream;
	FileUtils::openFile(fileName, inputFileStream);
	inputFileStream<<"14\tchunk682\t62251\t62526\trefTE_230\t2353\t2639\t0\t276\t84.3866\n";
	inputFileStream<<"14\tchunk682\t62662\t62891\trefTE_230\t2613\t2844\t0\t230\t90.75\n";

	inputFileStream<<"15\tchunk682\t62608\t62661\trefTE_230\t2538\t2595\t0\t54\t90.942\n";
	inputFileStream.close();
}

std::list<SDGString> Test_BLRMatchMapLoadUtils::generateExp_for_test_loadPath_with_join_begin_of_file()
{
	SDGString str1 = "1\t62251\t62891\t1\t2353\t2844\t0\t498\t87.2298\t517\n";
	SDGString str11= "\t1\t62251\t62526\t1\t2353\t2639\t0\t276\t84.3866\t286\n";
	SDGString str12 ="\t1\t62662\t62891\t1\t2613\t2844\t0\t230\t90.75\t231\n";

	SDGString str2 = "1\t62608\t62661\t1\t2538\t2595\t0\t54\t90.942\t57\n";
	SDGString str21 = "\t1\t62608\t62661\t1\t2538\t2595\t0\t54\t90.942\t57\n";
	
	std::list<SDGString> lPath;
	lPath.push_back(str1);
	lPath.push_back(str11);
	lPath.push_back(str12);
	lPath.push_back(str2);
	lPath.push_back(str21);
	return lPath;
}

void Test_BLRMatchMapLoadUtils::write_input_file_for_test_loadPath_with_join_middle_of_file(SDGString fileName)
{
	std::ofstream inputFileStream;
	FileUtils::openFile(fileName, inputFileStream);
	inputFileStream<<"13\tchunk682\t55604\t55920\trefTE_230\t2743\t3047\t0\t317\t91.3333\n";
	inputFileStream<<"14\tchunk682\t62251\t62526\trefTE_230\t2353\t2639\t0\t276\t84.3866\n";
	inputFileStream<<"14\tchunk682\t62662\t62891\trefTE_230\t2613\t2844\t0\t230\t90.75\n";
	inputFileStream<<"15\tchunk682\t62608\t62661\trefTE_230\t2538\t2595\t0\t54\t90.942\n";
	inputFileStream.close();
}

std::list<SDGString> Test_BLRMatchMapLoadUtils::generateExp_for_test_loadPath_with_join_middle_of_file()
{
	SDGString str1 = "1\t55604\t55920\t1\t2743\t3047\t0\t317\t91.3333\t317\n";
	SDGString str11 = "\t1\t55604\t55920\t1\t2743\t3047\t0\t317\t91.3333\t317\n";
	SDGString str2 = "1\t62251\t62891\t1\t2353\t2844\t0\t498\t87.2298\t517\n";
	SDGString str21= "\t1\t62251\t62526\t1\t2353\t2639\t0\t276\t84.3866\t286\n";
	SDGString str22 ="\t1\t62662\t62891\t1\t2613\t2844\t0\t230\t90.75\t231\n";
	SDGString str3 = "1\t62608\t62661\t1\t2538\t2595\t0\t54\t90.942\t57\n";
	SDGString str31 = "\t1\t62608\t62661\t1\t2538\t2595\t0\t54\t90.942\t57\n";
	
	std::list<SDGString> lPath;
	lPath.push_back(str1);
	lPath.push_back(str11);
	lPath.push_back(str2);
	lPath.push_back(str21);
	lPath.push_back(str22);
	lPath.push_back(str3);
	lPath.push_back(str31);
	return lPath;
}


