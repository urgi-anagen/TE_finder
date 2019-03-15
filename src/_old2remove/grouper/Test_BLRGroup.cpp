#include "Test_BLRGroup.h"
#include "Test_BLRGroupUtils.h"
#include "BLRGroup.h"
#include "BLRMemIdx.h"
#include "../../grouper.threads/BLRGrouperParameter.h"
#include "SDGString.h"
#include "FileUtils.h"

#include "BLRMatchMap.h"
#include "Test_BLRMatchMapUtils.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_BLRGroup );

void Test_BLRGroup::setUp()
{
}

void Test_BLRGroup::tearDown()
{
}

void Test_BLRGroup::test_group(void)
{
	BLRGrouperParameter para = Test_BLRGroupUtils::createParameter();	
	
	SDGString inputFileName = "input.align";	
	SDGString prefixFileName = "test_group";
	Test_BLRGroupUtils::write_input_file_for_test_group(inputFileName);
	para.setMatch_filename(inputFileName);
	para.setPrefix_filename(prefixFileName);

    	BLRGroup groups(&para);
    	groups.group(2);
	groups.save();	

	SDGString expFileName = "exp.group-test-group.c0.95.map";
	SDGString obsFileName = prefixFileName+".group.c0.95.map";
	
	Test_BLRGroupUtils::write_map_file_for_test_group(expFileName);
	bool obs = FileUtils::areTwoFilesIdentical(expFileName, obsFileName);
	bool exp = true;

	CPPUNIT_ASSERT_EQUAL(exp, obs);
	FileUtils::removeFile(inputFileName);
	FileUtils::removeFile(expFileName);
	FileUtils::removeFile(prefixFileName+".group.c0.95.set");
	FileUtils::removeFile(prefixFileName+".group.c0.95.param");
	FileUtils::removeFile(prefixFileName+".group.c0.95.map");
	FileUtils::removeFile(prefixFileName+".group.c0.95.fa");
	
}
 
void Test_BLRGroup::test_group_load_path(void)
{
	BLRGrouperParameter para = Test_BLRGroupUtils::createParameter();	

	SDGString inputFileName = "input.path";	
	SDGString prefixFileName = "test_group_load_path";
	Test_BLRGroupUtils::write_input_file_for_test_group_load_path(inputFileName);

	para.setPath_filename(inputFileName);
	para.setPrefix_filename(prefixFileName);

	para.setLoad_path(true);
	
 	BLRGroup groups(&para);
    	groups.group(2);
	groups.save();	

	SDGString expFileName = "exp.group-test-group-load-path.c0.95.map";
	SDGString obsFileName = prefixFileName+".group.c0.95.map";
	
	Test_BLRGroupUtils::write_map_file_for_test_group_load_path(expFileName);
	bool obs = FileUtils::areTwoFilesIdentical(expFileName, obsFileName);
	bool exp = true;

	CPPUNIT_ASSERT_EQUAL(exp, obs);
	FileUtils::removeFile(inputFileName);
	FileUtils::removeFile(expFileName);
	FileUtils::removeFile(prefixFileName+".group.c0.95.set");
	FileUtils::removeFile(prefixFileName+".group.c0.95.param");
	FileUtils::removeFile(prefixFileName+".group.c0.95.map");
	FileUtils::removeFile(prefixFileName+".group.c0.95.fa");
}


// before clustering, compare join between mapPath method and matcher -j -a -m invocation
// assertion is done manually
void Test_BLRGroup::test_compare_join_for_debug(void)
{
	BLRGrouperParameter para = Test_BLRGroupUtils::createParameter();	
	
	SDGString inputFileName = "input_test_group_for_debug.align";	
	SDGString prefixFileName = "test_group";
	Test_BLRGroupUtils::write_input_file_for_test_compare_join_for_debug(inputFileName);
	
	para.setMatch_filename(inputFileName);
	para.setPrefix_filename(prefixFileName);

	BLRMatchMap matchMap = BLRMatchMap(&para);
	matchMap.load();
	matchMap.mapPath(true, false, false, false, 0);

	std::list<RangePairSet> rp_list;
      	for( BLRMatchMap::MapPath::iterator m=matchMap.path_begin();
      	m!=matchMap.path_end(); m++ )
      	{
    	  	while(!m->second.empty())
    	  	{
    		  	RangePairSet rp=m->second.back();
    		  	m->second.pop_back();
    		  	rp_list.push_back(rp);
    	  	}
      	}
	
	
        for (std::list<RangePairSet>::iterator it = rp_list.begin();
	it != rp_list.end(); it++)
	{
		it->view();
	}

	SDGString cmd = "./../matcher/matcher2.25 -j -a -m " + inputFileName;
	std::system(cmd);
}

// compare intermediate rpsList between load() and loadPath()
void Test_BLRGroup::test_group_compare_intermediate_rpsList_between_2_loads(void)
{
	SDGString alignFilename = "compare_rpList_align.align";
	SDGString alignPrefixFilename = "compare_rpList_align";
	
	Test_BLRGroupUtils::write_align_file_for_test_group_compare_intermediate_rpsList_between_2_loads(alignFilename);
	
	
	BLRGrouperParameter paraForLoadAlignRun = Test_BLRGroupUtils::createParameter();	

	paraForLoadAlignRun.setMatch_filename(alignFilename);
	paraForLoadAlignRun.setPrefix_filename(alignPrefixFilename);
	
	// load	
 	BLRGroup groupsAlign(&paraForLoadAlignRun);
    	std::list<RangePairSet> rpsListAlign = groupsAlign.getRpsListAfterLoad(2);
	std::cout<<" "<<std::endl;	
	std::cout<<"rpsList Align "<<rpsListAlign.size()<<std::endl;	
    	for(std::list<RangePairSet>::iterator lrp_it=rpsListAlign.begin(); lrp_it!=rpsListAlign.end();lrp_it++){
            lrp_it->view(); 
    	}

	
	SDGString pathFilename = "compare_rpList_path.path";
	SDGString pathPrefixFilename = "compare_rpList_path.path";

	Test_BLRGroupUtils::write_path_file_for_test_group_compare_intermediate_rpsList_between_2_loads(pathFilename);
	BLRGrouperParameter paraForLoadPathRun = Test_BLRGroupUtils::createParameter();	
			
	paraForLoadPathRun.setPath_filename(pathFilename);
	paraForLoadPathRun.setPrefix_filename(pathPrefixFilename);
	paraForLoadPathRun.setLoad_path(true);

	// load path
	BLRGroup groupsPath(&paraForLoadPathRun);
    	std::list<RangePairSet> rpsListPath = groupsPath.getRpsListAfterLoad(2);
	std::cout<<" "<<std::endl;	
	std::cout<<"rpsList Path "<<rpsListPath.size()<<std::endl;	
    					
 	for(std::list<RangePairSet>::iterator lrp_it=rpsListPath.begin(); lrp_it!=rpsListPath.end();lrp_it++){
            lrp_it->view(); 
    	}


}
