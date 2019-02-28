#include "Test_F_grouper.h"
#include "SDGString.h"
#include "FileUtils.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_F_grouper );

void Test_F_grouper::setUp()
{
}

void Test_F_grouper::tearDown()
{
}

void Test_F_grouper::test_runAsScript_bigData( void ){
	SDGString GROUPER_DATA_SUFFIX = "/TE_finder/grouper/";
	SDGString GROUPER_DATA = std::getenv("REPET_DATA") + GROUPER_DATA_SUFFIX;

	SDGString inputFileName = "DmelChr4.fa-megablast-316416173.align";
	SDGString inputFileNameWithTestName ="DmelChr4.fa-megablast-316416173_bigData.align"; 
	SDGString cmd = "ln -s " + GROUPER_DATA + inputFileName + " " + inputFileNameWithTestName;
	std::system(cmd);
	
	SDGString genomeFileName = "DmelChr4.fa";
	SDGString genomeFileNameWithTestName = "DmelChr4_bigData.fa"; 
	cmd = "ln -s " + GROUPER_DATA + genomeFileName + " " + genomeFileNameWithTestName;
	std::system(cmd);

  	SDGString expFileName = "DmelChr4.fa-megablast-316416173.align.group.c0.95.map";
	SDGString expFileNameWithTestName = "exp_DmelChr4.fa-megablast-316416173_bigData.align.group.c0.95.map";
	cmd = "ln -s " + GROUPER_DATA + expFileName + " " + expFileNameWithTestName;
	std::system(cmd);

	SDGString outMapFileNameWithTestName = "DmelChr4.fa-megablast-316416173_bigData.align.group.c0.95.map";
	SDGString outDotFileNameWithTestName = "DmelChr4.fa-megablast-316416173_bigData.align.group.c0.95.cluster.dot";
	SDGString outSetFileNameWithTestName = "DmelChr4.fa-megablast-316416173_bigData.align.group.c0.95.set";
	SDGString outTxtFileNameWithTestName = "DmelChr4.fa-megablast-316416173_bigData.align.group.c0.95.txt";
	SDGString outParamFileNameWithTestName = "DmelChr4.fa-megablast-316416173_bigData.align.group.c0.95.param";
	SDGString outFaFileNameWithTestName = "DmelChr4.fa-megablast-316416173_bigData.align.group.c0.95.fa";

        cmd= "./grouper2.25 -m "+ inputFileNameWithTestName + " -q " + genomeFileNameWithTestName + " -s " + genomeFileNameWithTestName +" -j -Z 3 -v 2";
	std::cout<<" "<<std::endl;
	std::cout<<cmd<<std::endl;
  	std::system(cmd);

	bool obs = FileUtils::areTwoFilesIdentical(expFileNameWithTestName, outMapFileNameWithTestName);
	bool exp = true;

	CPPUNIT_ASSERT_EQUAL(exp, obs);
	
	FileUtils::removeFile(inputFileNameWithTestName);

	FileUtils::removeFile(outMapFileNameWithTestName);
	FileUtils::removeFile(outDotFileNameWithTestName);
	FileUtils::removeFile(outSetFileNameWithTestName);
	FileUtils::removeFile(outTxtFileNameWithTestName);
	FileUtils::removeFile(outParamFileNameWithTestName);
	FileUtils::removeFile(outFaFileNameWithTestName);
	FileUtils::removeFile(expFileNameWithTestName);
	FileUtils::removeFile(genomeFileNameWithTestName);
}

void Test_F_grouper::test_runAsScript_smallData_filter_length_40_46( void ){
	SDGString GROUPER_DATA_SUFFIX = "/TE_finder/grouper/";
	SDGString GROUPER_DATA = std::getenv("REPET_DATA") + GROUPER_DATA_SUFFIX;

	SDGString inputFileName = "DmelChr4.fa-megablast-316416173-small.align";
	SDGString inputFileNameWithTestName ="DmelChr4.fa-megablast-316416173_smallData.align"; 
	SDGString cmd = "ln -s " + GROUPER_DATA + inputFileName + " " + inputFileNameWithTestName;
	std::system(cmd);
	
	SDGString genomeFileName = "DmelChr4.fa";
	SDGString genomeFileNameWithTestName = "DmelChr4_smallData.fa"; 
	cmd = "ln -s " + GROUPER_DATA + genomeFileName + " " + genomeFileNameWithTestName;
	std::system(cmd);

  	SDGString expFileName = "DmelChr4.fa-megablast-316416173-small.align.group.c0.95.map";
	SDGString expFileNameWithTestName = "exp_DmelChr4.fa-megablast-316416173_smallData.align.group.c0.95.map";
	cmd = "ln -s " + GROUPER_DATA + expFileName + " " + expFileNameWithTestName;
	std::system(cmd);

	SDGString outMapFileNameWithTestName = "DmelChr4.fa-megablast-316416173_smallData.align.group.c0.95.map";
	SDGString outDotFileNameWithTestName = "DmelChr4.fa-megablast-316416173_smallData.align.group.c0.95.cluster.dot";
	SDGString outSetFileNameWithTestName = "DmelChr4.fa-megablast-316416173_smallData.align.group.c0.95.set";
	SDGString outTxtFileNameWithTestName = "DmelChr4.fa-megablast-316416173_smallData.align.group.c0.95.txt";
	SDGString outParamFileNameWithTestName = "DmelChr4.fa-megablast-316416173_smallData.align.group.c0.95.param";
	SDGString outFaFileNameWithTestName = "DmelChr4.fa-megablast-316416173_smallData.align.group.c0.95.fa";

        cmd= "./grouper2.25 -m "+ inputFileNameWithTestName + " -q " + genomeFileNameWithTestName + " -s " + genomeFileNameWithTestName +" -j -Z 3";
	cmd = cmd + " -v 2";
	std::cout<<" "<<std::endl;
	std::cout<<cmd<<std::endl;
  	std::system(cmd);

	bool obs = FileUtils::areTwoFilesIdentical(expFileNameWithTestName, outMapFileNameWithTestName);
	bool exp = true;

	CPPUNIT_ASSERT_EQUAL(exp, obs);
	
	FileUtils::removeFile(inputFileNameWithTestName);

	FileUtils::removeFile(outMapFileNameWithTestName);
	FileUtils::removeFile(outDotFileNameWithTestName);
	FileUtils::removeFile(outSetFileNameWithTestName);
	FileUtils::removeFile(outTxtFileNameWithTestName);
	FileUtils::removeFile(outParamFileNameWithTestName);
	FileUtils::removeFile(outFaFileNameWithTestName);
	FileUtils::removeFile(expFileNameWithTestName);
	FileUtils::removeFile(genomeFileNameWithTestName);
}

void Test_F_grouper::test_runAsScript_loadPath_option_grouper_only( void ){
	SDGString GROUPER_DATA_SUFFIX = "/TE_finder/grouper/";
	SDGString GROUPER_DATA = std::getenv("REPET_DATA") + GROUPER_DATA_SUFFIX;
	
	SDGString alignFileName = "DmelChr4.fa-megablast-316416173.align";
	SDGString cmd = "ln -s " + GROUPER_DATA + alignFileName + " " + alignFileName;
	std::system(cmd);
	
	SDGString genomeFileName = "DmelChr4.fa";
	cmd = "ln -s " + GROUPER_DATA + genomeFileName + " " + genomeFileName;
	std::system(cmd);
	//run grouper with align file
	
	cmd = "./grouper2.26 -v 2 -d 10 -g 1 -c 2 -m "+ alignFileName + " -j -q " + genomeFileName;
	std::system(cmd);

	SDGString mapFileNameGrp = alignFileName + ".group.c0.95.map";
	SDGString dotFileNameGrp = alignFileName + ".group.c0.95.cluster.dot";
	SDGString setFileNameGrp = alignFileName + ".group.c0.95.set";
	SDGString paramFileNameGrp = alignFileName + ".group.c0.95.param";
	SDGString faFileNameGrp = alignFileName + ".group.c0.95.fa";
	SDGString txtFileNameGrp = alignFileName + ".group.c0.95.txt";
	SDGString gPathFileName = alignFileName + ".gpath";	
	SDGString gPathAttrFileName = alignFileName + ".gpath.attr";
	SDGString rpsListNoLoadPath = alignFileName + ".rpsListNoLoadPath";
	SDGString rpsListNoLoadPathAttr = alignFileName + ".rpsListNoLoadPath.attr";

	// run grouper with path file
	cmd = "./grouper2.26 -d 10 -g 1 -p "+ gPathFileName + " -q " + genomeFileName;
	std::system(cmd);

	SDGString mapFileNameMatGrp = gPathFileName + ".group.c0.95.map";
	SDGString dotFileNameMatGrp = gPathFileName + ".group.c0.95.cluster.dot";
	SDGString setFileNameMatGrp = gPathFileName + ".group.c0.95.set";
	SDGString paramFileNameMatGrp = gPathFileName + ".group.c0.95.param";
	SDGString faFileNameMatGrp = gPathFileName + ".group.c0.95.fa";
	SDGString txtFileNameMatGrp = gPathFileName + ".group.c0.95.txt";
	
	bool obs = FileUtils::areTwoFilesIdentical(txtFileNameMatGrp, txtFileNameGrp);
	bool exp = true;
        
	CPPUNIT_ASSERT_EQUAL(exp, obs);


        FileUtils::removeFile(mapFileNameMatGrp);
        FileUtils::removeFile(dotFileNameMatGrp);
        FileUtils::removeFile(setFileNameMatGrp);
        FileUtils::removeFile(paramFileNameMatGrp);
        FileUtils::removeFile(faFileNameMatGrp);
        FileUtils::removeFile(txtFileNameMatGrp);

        FileUtils::removeFile(mapFileNameGrp);
        FileUtils::removeFile(dotFileNameGrp);
        FileUtils::removeFile(setFileNameGrp);
        FileUtils::removeFile(paramFileNameGrp);
        FileUtils::removeFile(faFileNameGrp);
        FileUtils::removeFile(txtFileNameGrp);

        FileUtils::removeFile(alignFileName);
        FileUtils::removeFile(gPathFileName);
	FileUtils::removeFile(gPathAttrFileName);
	FileUtils::removeFile(rpsListNoLoadPath);
	FileUtils::removeFile(rpsListNoLoadPathAttr);
        FileUtils::removeFile(genomeFileName);


}

// compare run matcher (join) + grouper -j
// with run grouper -m
void Test_F_grouper::test_runAsScript_loadPath_option( void ){
	SDGString GROUPER_DATA_SUFFIX = "/TE_finder/grouper/";
	SDGString GROUPER_DATA = std::getenv("REPET_DATA") + GROUPER_DATA_SUFFIX;
	
	SDGString alignFileName = "DmelChr4.fa-megablast-316416173.align";
	SDGString cmd = "ln -s " + GROUPER_DATA + alignFileName + " " + alignFileName;
	std::system(cmd);
	
	SDGString genomeFileName = "DmelChr4.fa";
	cmd = "ln -s " + GROUPER_DATA + genomeFileName + " " + genomeFileName;
	std::system(cmd);

	// run combo matcher + grouper
	cmd = "./../matcher/matcher2.25 -g 1 -d 10 -c 2 -a -j -m " + alignFileName;
	std::system(cmd);
	
	FileUtils::removeFile(alignFileName+".match.param");
	FileUtils::removeFile(alignFileName+".match.map");

	SDGString pathFileName = alignFileName+".match.path";

	cmd = "./grouper2.26 -g 1 -d 10 -j -p "+ pathFileName + " -q " + genomeFileName;
	std::system(cmd);

	SDGString mapFileNameMatGrp = pathFileName + ".group.c0.95.map";
	SDGString dotFileNameMatGrp = pathFileName + ".group.c0.95.cluster.dot";
	SDGString setFileNameMatGrp = pathFileName + ".group.c0.95.set";
	SDGString paramFileNameMatGrp = pathFileName + ".group.c0.95.param";
	SDGString faFileNameMatGrp = pathFileName + ".group.c0.95.fa";
	SDGString txtFileNameMatGrp = pathFileName + ".group.c0.95.txt";


	//run grouper with align file
	cmd = "./grouper2.26 -v 2 -g 1 -d 10 -c 2 -m "+ alignFileName + " -j -q " + genomeFileName;
	std::system(cmd);

	SDGString mapFileNameGrp = alignFileName + ".group.c0.95.map";
	SDGString dotFileNameGrp = alignFileName + ".group.c0.95.cluster.dot";
	SDGString setFileNameGrp = alignFileName + ".group.c0.95.set";
	SDGString paramFileNameGrp = alignFileName + ".group.c0.95.param";
	SDGString faFileNameGrp = alignFileName + ".group.c0.95.fa";
	SDGString txtFileNameGrp = alignFileName + ".group.c0.95.txt";
	SDGString gPathFileName = alignFileName + ".gpath";
        SDGString gPathAttrFileName = alignFileName + ".gpath.attr";
        SDGString rpsListNoLoadPath = alignFileName + ".rpsListNoLoadPath";
        SDGString rpsListNoLoadPathAttr = alignFileName + ".rpsListNoLoadPath.attr";
	SDGString pathFileNameAttr = alignFileName+".match.path.attr";

	bool obs = FileUtils::areTwoFilesIdentical(txtFileNameMatGrp, txtFileNameGrp);
	bool exp = true;

	CPPUNIT_ASSERT_EQUAL(exp, obs);
	
	FileUtils::removeFile(mapFileNameMatGrp);
	FileUtils::removeFile(dotFileNameMatGrp);
	FileUtils::removeFile(setFileNameMatGrp);
	FileUtils::removeFile(paramFileNameMatGrp);
	FileUtils::removeFile(faFileNameMatGrp);
	FileUtils::removeFile(txtFileNameMatGrp);

	FileUtils::removeFile(mapFileNameGrp);
	FileUtils::removeFile(dotFileNameGrp);
	FileUtils::removeFile(setFileNameGrp);
	FileUtils::removeFile(paramFileNameGrp);
	FileUtils::removeFile(faFileNameGrp);
	FileUtils::removeFile(txtFileNameGrp);

	FileUtils::removeFile(alignFileName);
	FileUtils::removeFile(pathFileName);
	FileUtils::removeFile(genomeFileName);
	FileUtils::removeFile(alignFileName);
        FileUtils::removeFile(gPathFileName);
        FileUtils::removeFile(gPathAttrFileName);
        FileUtils::removeFile(rpsListNoLoadPath);
        FileUtils::removeFile(rpsListNoLoadPathAttr);
	FileUtils::removeFile(pathFileNameAttr);

}

// compare run matcher (join) + grouper -j
// with run grouper -m
void Test_F_grouper::test_runAsScript_loadPath_option_small_data( void )
{
	SDGString GROUPER_DATA_SUFFIX = "/TE_finder/grouper/";
	SDGString GROUPER_DATA = std::getenv("REPET_DATA") + GROUPER_DATA_SUFFIX;
	

	SDGString alignFileName = "input_align_loadPath_small_data.align";

	write_input_align_file_for_test_group_loadPath_option_small_data(alignFileName);
	
	SDGString genomeFileName = "DmelChr4.fa";
	SDGString cmd = "ln -s " + GROUPER_DATA + genomeFileName + " " + genomeFileName;
	std::system(cmd);

	// run combo matcher + grouper
	
	cmd = "./../matcher/matcher2.25 -g 1 -d 20 -c 2 -a -j -m " + alignFileName;
	std::system(cmd);
	
	FileUtils::removeFile(alignFileName+".match.param");
	FileUtils::removeFile(alignFileName+".match.map");

	SDGString pathFileName = alignFileName+".match.path";

	//cmd = "./grouper2.25 -j -p "+ pathFileName + " -q " + genomeFileName;
	
	cmd = "./grouper2.25  -p "+ pathFileName + " -q " + genomeFileName;
	std::system(cmd);

	SDGString mapFileNameMatGrp = pathFileName + ".group.c0.95.map";
	SDGString dotFileNameMatGrp = pathFileName + ".group.c0.95.cluster.dot";
	SDGString setFileNameMatGrp = pathFileName + ".group.c0.95.set";
	SDGString paramFileNameMatGrp = pathFileName + ".group.c0.95.param";
	SDGString faFileNameMatGrp = pathFileName + ".group.c0.95.fa";
	SDGString txtFileNameMatGrp = pathFileName + ".group.c0.95.txt";
	

	//run grouper with align file
	cmd = "./grouper2.25  -g 1 -d 20 -c 2 -m "+ alignFileName + " -j -q " + genomeFileName;
	std::system(cmd);

	SDGString mapFileNameGrp = alignFileName + ".group.c0.95.map";
	SDGString dotFileNameGrp = alignFileName + ".group.c0.95.cluster.dot";
	SDGString setFileNameGrp = alignFileName + ".group.c0.95.set";
	SDGString paramFileNameGrp = alignFileName + ".group.c0.95.param";
	SDGString faFileNameGrp = alignFileName + ".group.c0.95.fa";
	SDGString txtFileNameGrp = alignFileName + ".group.c0.95.txt";

	
	bool obs = FileUtils::areTwoFilesIdentical(txtFileNameMatGrp, txtFileNameGrp);
	bool exp = true;

	CPPUNIT_ASSERT_EQUAL(exp, obs);

	
	FileUtils::removeFile(mapFileNameMatGrp);
	FileUtils::removeFile(dotFileNameMatGrp);
	FileUtils::removeFile(setFileNameMatGrp);
	FileUtils::removeFile(paramFileNameMatGrp);
	FileUtils::removeFile(faFileNameMatGrp);
	FileUtils::removeFile(txtFileNameMatGrp);

	FileUtils::removeFile(mapFileNameGrp);
	FileUtils::removeFile(dotFileNameGrp);
	FileUtils::removeFile(setFileNameGrp);
	FileUtils::removeFile(paramFileNameGrp);
	FileUtils::removeFile(faFileNameGrp);
	FileUtils::removeFile(txtFileNameGrp);

	FileUtils::removeFile(alignFileName);
	FileUtils::removeFile(pathFileName);
	FileUtils::removeFile(genomeFileName);


}



void Test_F_grouper::write_input_align_file_for_test_group_loadPath_option_small_data(SDGString alignFileName)
{
	std::ofstream alignFileStream;
	FileUtils::openFile(alignFileName, alignFileStream);
	alignFileStream<<"dmel_chr4\t13293\t13592\tdmel_chr4\t74091\t74394\t8e-75\t281\t87.18\n";
	alignFileStream<<"dmel_chr4\t13590\t13726\tdmel_chr4\t74436\t74572\t1e-52\t208\t94.2\n";
	alignFileStream<<"dmel_chr4\t13718\t13762\tdmel_chr4\t1260453\t1260409\t4e-12\t74\t95.56\n";
	alignFileStream<<"dmel_chr4\t13609\t13726\tdmel_chr4\t1260568\t1260459\t2e-17\t92\t85.59\n";
	alignFileStream<<"dmel_chr4\t13305\t13592\tdmel_chr4\t1197346\t11976366\t8e-72\t272\t87.16\n";
	alignFileStream<<"dmel_chr4\t13590\t13746\tdmel_chr4\t1197677\t1197834\t4e-55\t216\t92.41\n";
	alignFileStream.close();
}

/*
void Test_F_grouper::write_input_path_file_for_test_group_loadPath_option_small_data(SDGstring pathFileName)
{
	std::ofstream pathFileStream;
	FileUtils::openFile(pathFileName, pathFileStream);
pathFileStream<<"1\tdmel_chr4\t13293\t13592\tdmel_chr4\t74091\t74394\t8e-75\t281\t87.18\n";
pathFileStream<<"1\tdmel_chr4\t13590\t13726\tdmel_chr4\t74436\t74572\t1e-52\t208\t94.2\n";
pathFileStream<<"2\tdmel_chr4\t13718 \t13762\tdmel_chr4\t1260453\t1260409\t4e-12\t74\t95.56\n";
2       dmel_chr4       13609   13726   dmel_chr4       1260568 1260459 2e-17   92      85.59
3       dmel_chr4       13305   13592   dmel_chr4       1197346 11976366        8e-72   272     87.16
4       dmel_chr4       13590   13746   dmel_chr4       1197677 1197834 4e-55   216     92.41

	
	
}*/

void Test_F_grouper::write_input_file_for_test_group_load_path(SDGString pathFileName)
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

void Test_F_grouper::write_map_file_for_test_group_load_path(SDGString mapFileName)
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
