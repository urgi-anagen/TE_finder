#include <stack>
#include <sstream>
#include "Test_BLRMatchMapUtils.h"
#include "SDGString.h"
#include "FileUtils.h"
#include "../../matcher/BLRMatcherParameter.h"
#include "BLRMatchMap.h"
#include "Range.h"
#include "RangePair.h"
#include "RangePairSet.h"

void Test_BLRMatchMapUtils::viewMapAlign(BLRMatchMap::MapAlign mapAlign){
  for(BLRMatchMap::MapAlign::iterator m=mapAlign.begin(); m!=mapAlign.end();m++)
  {
      while (!m->second.empty())
      {
        BLRMatchMap::Key key = m->first;
        SDGString numQuery = SDGString(key.first);
        SDGString numSubject  = SDGString(key.second);
        //std::cout<<" "<<std::endl;
        //std::cout<<"Num query: "<<numQuery<<" Num subject: "<<numSubject<<std::endl;
        RangePair range_pair = m->second.back();
        SDGString numChr = SDGString(range_pair.getRangeQ().getNumChr());
        SDGString nameSeq = range_pair.getRangeQ().getNameSeq();
        //std::cout<<"Query Name Seq: "<<nameSeq<<"Query Num chr: "<<numChr<<std::endl; 
	      m->second.pop_back();
        range_pair.view();
      }
  };
}

void Test_BLRMatchMapUtils::viewMapPath(BLRMatchMap::MapPath mapPath){
  for(BLRMatchMap::MapPath::iterator m=mapPath.begin(); m!=mapPath.end();m++)
  {
      while (!m->second.empty())
      {
        BLRMatchMap::Key key = m->first;
        SDGString keyFirst = SDGString(key.first);
        SDGString keySecond  = SDGString(key.second);
        std::cout<<" "<<std::endl;
        std::cout<<"("<<keyFirst<<","<<keySecond<<")"<<std::endl;
        std::cout<<" "<<std::endl;
        RangePairSet range_pair = m->second.back();
        SDGString numChr = SDGString(range_pair.getRangeQ().getNumChr());
        SDGString nameSeq = range_pair.getRangeQ().getNameSeq();
	      m->second.pop_back();
        range_pair.view();
      }
  };
}

void Test_BLRMatchMapUtils::viewMapPathWithLabel(BLRMatchMap::MapPath mapPath){
  for(BLRMatchMap::MapPath::iterator m=mapPath.begin(); m!=mapPath.end();m++)
  {
      while (!m->second.empty())
      {
        BLRMatchMap::Key key = m->first;
        SDGString keyFirst = SDGString(key.first);
        SDGString keySecond  = SDGString(key.second);
        std::cout<<" "<<std::endl;
        std::cout<<"("<<keyFirst<<","<<keySecond<<")"<<std::endl;
        std::cout<<" "<<std::endl;
        RangePairSet range_pair = m->second.back();
        SDGString numChr = SDGString(range_pair.getRangeQ().getNumChr());
        SDGString nameSeq = range_pair.getRangeQ().getNameSeq();
	      m->second.pop_back();
        range_pair.viewWithLabel();
      }
  };
}
void Test_BLRMatchMapUtils::viewNum2Name(std::map<long,std::string> num2Name){
	 for (std::map<long,std::string>::iterator it=num2Name.begin(); it!=num2Name.end(); ++it)
    		std::cout<<it->first<<" => "<< it->second <<std::endl;
  
}

// TODO duplication with Test_BLRRangePairSetUtils
void Test_BLRMatchMapUtils::viewRangePairSetList(std::list<RangePairSet> rpsList){
    for(std::list<RangePairSet>::iterator lrp_it=rpsList.begin(); lrp_it!=rpsList.end();lrp_it++){
            lrp_it->view(); 
    }
}


void Test_BLRMatchMapUtils::writeInputFile(void)
{
  SDGString inputFileName = "match.align";
  std::ofstream inputFileStream;
  FileUtils::openFile(inputFileName, inputFileStream);
  inputFileStream << "chunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.836100\n";
  inputFileStream << "chunk682\t105091\t105206\trefTE_251\t3351\t3238\t0\t116\t76.106200\n";
  inputFileStream << "chunk682\t109527\t109818\trefTE_251\t3005\t3316\t0\t292\t82.846700\n";
  inputFileStream << "chunk682\t109601\t109708\trefTE_251\t3239\t3338\t3e-28\t108\t81.818200\n";
  inputFileStream << "chunk682\t109951\t110569\trefTE_266\t119\t689\t9.4e-19\t619\t77.768000\n";
  inputFileStream << "chunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.860400\n";
  inputFileStream << "chunk682\t110567\t110878\trefTE_230\t2833\t2532\t0\t312\t77.207600\n";
  inputFileStream.close();
}

void Test_BLRMatchMapUtils::writeInputFileWithTwoMatchesToBeJoined(void)
{
  	SDGString inputFileName = "match.align";
	std::ofstream inputFileStream;
	FileUtils::openFile(inputFileName, inputFileStream);
  	inputFileStream<<"chunk682\t105085\t105353\trefTE_251\t3202\t2921\t0\t269\t84.836100\n";
  	inputFileStream<<"chunk682\t109985\t110398\trefTE_251\t2429\t2033\t0\t414\t77.860400\n";
	inputFileStream.close();
}

BLRMatcherThreadsParameter Test_BLRMatchMapUtils::createParameter(void){
  return createParameter("match.param");
}

BLRMatcherThreadsParameter Test_BLRMatchMapUtils::createParameter(SDGString inputFileName){
  BLRMatcherThreadsParameter para;
  para.setMatch_filename(inputFileName);    
  para.setJoin_frag(true);
  // Enable '-x' option
  para.setCleanBefore(false);
  para.setCleanAfter(true); 
  return para;
}

BLRMatcherThreadsParameter Test_BLRMatchMapUtils::createParameterWithThreads(){
  BLRMatcherThreadsParameter para;
  para.setMatch_filename("match.param");
  para.setJoin_frag(true);
  // Enable '-x' option
  para.setCleanBefore(false);
  para.setCleanAfter(true);
  para.setNbThread(2);
  return para;
}

BLRMatchMap::MapPath Test_BLRMatchMapUtils::createExpMapPath_for_mapPath(void){
  BLRMatchMap::MapPath mapPath;
  
  SDGString line1 = " \t105085\t110569\t \t0\t0\t0\t888\t0";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeS().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNameSeq("");
  //rangePairSet1.setLength(888);

  SDGString line11 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeS().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNameSeq("");
  //rangePair11.setLength(269);


  //SDGString line12 = " \t109951\t109984\t \t119\t149\t9.4e-19\t34\t77.768";
  SDGString line12 = " \t109951\t109984\t \t149\t119\t9.4e-19\t34\t77.768";
  RangePair rangePair12 = RangePair(line12);
  rangePair12.getRangeQ().setNumChr(1);
  rangePair12.getRangeS().setNumChr(1);
  rangePair12.getRangeQ().setNameSeq("");
  rangePair12.getRangeS().setNameSeq("");
  //rangePair12.setLength(34);

  SDGString line13 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair13 = RangePair(line13);
  rangePair13.getRangeQ().setNumChr(1);
  rangePair13.getRangeS().setNumChr(1);
  rangePair13.getRangeQ().setNameSeq("");
  rangePair13.getRangeS().setNameSeq("");
  //rangePair13.setLength(414);

  //SDGString line14 = " \t110399\t110569\t \t532\t689\t9.4e-19\t171\t77.768";
  SDGString line14 = " \t110399\t110569\t \t689\t532\t9.4e-19\t171\t77.768";
  RangePair rangePair14 = RangePair(line14);
  rangePair14.getRangeQ().setNumChr(1);
  rangePair14.getRangeS().setNumChr(1);
  rangePair14.getRangeQ().setNameSeq("");
  rangePair14.getRangeS().setNameSeq("");
  //rangePair14.setLength(171);


  // push range to range pair
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);
  rpList1.push_back(rangePair12);
  rpList1.push_back(rangePair13);
  rpList1.push_back(rangePair14);
  

  // set match part to match
  rangePairSet1.setPathDirectly(rpList1);   
  
  
  // insert in mapPath
  std::list<RangePairSet>& rpsList = mapPath[BLRMatchMap::Key(1,1)];
  std::list<RangePairSet>::iterator r=std::lower_bound(rpsList.begin(), rpsList.end(),rangePairSet1);
  rpsList.insert(r, rangePairSet1);   
 

 
  SDGString line2 = " \t109527\t109818\t \t3005\t3316\t0\t292\t82.8467";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNameSeq("");
 //rangePairSet2.setLength(292);

  SDGString line21 = " \t109527\t109818\t \t3005\t3316\t0\t292\t82.8467";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNameSeq("");
  //rangePair21.setLength(292);

  // push range to range pair
  std::list<RangePair> rpList2;
  rpList2.push_back(rangePair21);
  
  rangePairSet2.setPathDirectly(rpList2);   
  
  // insert in mapPath
  std::list<RangePairSet>& rpsList2 = mapPath[BLRMatchMap::Key(1,1)];
  std::list<RangePairSet>::iterator r2=std::lower_bound(rpsList2.begin(), rpsList2.end(),rangePairSet2);
  rpsList2.insert(r2, rangePairSet2);   


  SDGString line3 = " \t105091\t105206\t \t3351\t3238\t0\t116\t76.1062";
  RangePairSet rangePairSet3 = RangePairSet(line3);
  rangePairSet3.getRangeQ().setNumChr(1);
  rangePairSet3.getRangeS().setNumChr(1);
  rangePairSet3.getRangeQ().setNameSeq("");
  rangePairSet3.getRangeS().setNameSeq("");
  //rangePairSet3.setLength(116);
  
  SDGString line31 = " \t105091\t105206\t \t3351\t3238\t0\t116\t76.1062";
  RangePair rangePair31 = RangePair(line31);
  rangePair31.getRangeQ().setNumChr(1);
  rangePair31.getRangeS().setNumChr(1);
  rangePair31.getRangeQ().setNameSeq("");
  rangePair31.getRangeS().setNameSeq("");
  //rangePair31.setLength(116);
  
  // push range to range pair
  std::list<RangePair> rpList3;
  rpList3.push_back(rangePair31);
  
  rangePairSet3.setPathDirectly(rpList3);   
 
  // insert in mapPath
  std::list<RangePairSet>& rpsList3 = mapPath[BLRMatchMap::Key(1,1)];
  std::list<RangePairSet>::iterator r3=std::lower_bound(rpsList3.begin(), rpsList3.end(),rangePairSet3);
  rpsList3.insert(r3, rangePairSet3);   
  
  SDGString line4 = " \t110567\t110878\t \t2833\t2532\t0\t312\t77.2076";
  RangePairSet rangePairSet4 = RangePairSet(line4);
  rangePairSet4.getRangeQ().setNumChr(1);
  rangePairSet4.getRangeS().setNumChr(3);
  rangePairSet4.getRangeQ().setNameSeq("");
  rangePairSet4.getRangeS().setNameSeq("");
  //rangePairSet4.setLength(116);
   
  SDGString line41 = " \t110567\t110878\t \t2833\t2532\t0\t312\t77.2076";
  RangePair rangePair41 = RangePair(line41);
  rangePair41.getRangeQ().setNumChr(1);
  rangePair41.getRangeS().setNumChr(3);
  rangePair41.getRangeQ().setNameSeq("");
  rangePair41.getRangeS().setNameSeq("");
  //rangePair41.setLength(116);
   
   // push range to range pair
  std::list<RangePair> rpList4;
  rpList4.push_back(rangePair41);
  
  rangePairSet4.setPathDirectly(rpList4);   
   
  // insert in mapPath
  std::list<RangePairSet>& rpsList4 = mapPath[BLRMatchMap::Key(1,3)];
  std::list<RangePairSet>::iterator r4=std::lower_bound(rpsList4.begin(), rpsList4.end(),rangePairSet4);
  rpsList3.insert(r4, rangePairSet4);   

 
  return mapPath; 
}

BLRMatchMap::MapAlign Test_BLRMatchMapUtils::createExpMapAlign_for_clean_conflicts(void){
  SDGString match_file = "match.align";
  writeInputFile();
  BLRMatcherThreadsParameter para = createParameter();
  BLRMatchMap matchMap(&para);
  matchMap.clear();

  SDGString line1 = " \t109601\t109708\t \t3239\t3338\t3e-28\t108\t81.8182";
  RangePair rangePair1 = RangePair(line1);
  rangePair1.getRangeQ().setNumChr(1);
  rangePair1.getRangeQ().setNameSeq("chunk682");
  rangePair1.getRangeS().setNumChr(1);
  rangePair1.getRangeS().setNameSeq("refTE_251");

  matchMap.insert(rangePair1);
 
  SDGString line2 = " \t109527\t109818\t \t3005\t3316\t0\t292\t82.8467";
  RangePair rangePair2 = RangePair(line2);
  rangePair2.getRangeQ().setNumChr(1);
  rangePair2.getRangeQ().setNameSeq("chunk682");
  rangePair2.getRangeS().setNumChr(1);
  rangePair2.getRangeS().setNameSeq("refTE_251");

  matchMap.insert(rangePair2);
 
  SDGString line3 = " \t105091\t105206\t \t3351\t3238\t0\t116\t76.1062";
  RangePair rangePair3 = RangePair(line3);
  rangePair3.getRangeQ().setNumChr(1);
  rangePair3.getRangeQ().setNameSeq("chunk682");
  rangePair3.getRangeS().setNumChr(1);
  rangePair3.getRangeS().setNameSeq("refTE_251");

  matchMap.insert(rangePair3);

  SDGString line4 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair4 = RangePair(line4);
  rangePair4.getRangeQ().setNumChr(1);
  rangePair4.getRangeQ().setNameSeq("chunk682");
  rangePair4.getRangeS().setNumChr(1);
  rangePair4.getRangeS().setNameSeq("refTE_251");

  matchMap.insert(rangePair4);

  SDGString line5 = " \t109951\t110569\t \t119\t689\t9.4e-19\t619\t77.768";
  RangePair rangePair5 = RangePair(line5);
  rangePair5.getRangeQ().setNumChr(1);
  rangePair5.getRangeQ().setNameSeq("chunk682");
  rangePair5.getRangeS().setNumChr(2);
  rangePair5.getRangeS().setNameSeq("refTE_266");

  matchMap.insert(rangePair5);

  SDGString line6 = " \t110570\t110878\t \t2830\t2532\t0\t309\t77.2076";
  RangePair rangePair6 = RangePair(line6);
  rangePair6.getRangeQ().setNumChr(1);
  rangePair6.getRangeQ().setNameSeq("chunk682");
  rangePair6.getRangeS().setNumChr(3);
  rangePair6.getRangeS().setNameSeq("refTE_230");
  matchMap.insert(rangePair6);

  return matchMap.getMapAlign();
     
}

std::list<RangePairSet> Test_BLRMatchMapUtils::createRpList_for_test_add_clean_path_all_S(void){
    BLRMatcherThreadsParameter para = createParameter(); 
    BLRMatchMap matchMap(&para);
    bool joiningParameter = true;
    bool cleanBefore = false;
    bool cleanAfter = true;
    int verboseParameter = 0;
    SDGString match_file = "match.align";
    writeInputFile();
    matchMap.loadAlign(match_file,0);
    BLRMatchMap::MapAlign mapAlignBefore = matchMap.getMapAlign();
    matchMap.mapPathJoinOnlyForTest(joiningParameter, cleanBefore, cleanAfter, verboseParameter);
    BLRMatchMap::MapPath map_path = matchMap.getMapPath(); 
    std::list<RangePairSet> rp_list;
    for(BLRMatchMap::MapPath::iterator m=map_path.begin(); m!=map_path.end();m++)
    {
      while(!m->second.empty())
	{
	  RangePairSet rp=m->second.back();
	  m->second.pop_back();
	  rp.computeScoreWithDynaProg();
	  rp_list.push_back(rp);
	}
    }
    return rp_list;
}

BLRMatchMap::MapAlign Test_BLRMatchMapUtils::createMapAlign_instance_for_test_mapAlign_Equality(void){
  SDGString match_file = "match.align";
  writeInputFile();
  BLRMatcherThreadsParameter para = createParameter();
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

  BLRMatchMap::MapAlign mapAlign = matchMap.getMapAlign();
  return mapAlign; 
}

std::list<RangePairSet> Test_BLRMatchMapUtils::createRpListForTest_extractRangePairSetListFromMapPath(void){
    // match  
    SDGString line1 = " \t105091\t105206\t \t3351\t3238\t0\t116\t76.1062";
    RangePairSet rangePairSet1 = RangePairSet(line1);
    rangePairSet1.getRangeQ().setNumChr(1);
    rangePairSet1.getRangeQ().setNameSeq("");
    rangePairSet1.getRangeS().setNumChr(1);
    rangePairSet1.getRangeS().setNameSeq("");
    rangePairSet1.getRangeS().setNameSeq("");
    rangePairSet1.setLength(116);
    
    // match part
    SDGString line11 = " \t105091\t105206\t \t3351\t3238\t0\t116\t76.1062";
    RangePair rangePair11 = RangePair(line11);
    rangePair11.getRangeQ().setNumChr(1);
    rangePair11.getRangeQ().setNameSeq("");
    rangePair11.getRangeS().setNumChr(1);
    rangePair11.getRangeS().setNameSeq("");
    rangePair11.setLength(116);

    // set match part to match
    std::list<RangePair> rpList1;
    rpList1.push_back(rangePair11);
    rangePairSet1.setPathDirectly(rpList1);   

    // match    
    SDGString line2 = " \t109601\t109708\t \t3239\t3338\t3e-28\t108\t81.8182";
    RangePairSet rangePairSet2 = RangePairSet(line2);
    rangePairSet2.getRangeQ().setNumChr(1);
    rangePairSet2.getRangeQ().setNameSeq("");
    rangePairSet2.getRangeS().setNumChr(1);
    rangePairSet2.getRangeS().setNameSeq("");
    rangePairSet2.setLength(108);
    // match part
    SDGString line22 = " \t109601\t109708\t \t3239\t3338\t3e-28\t108\t81.8182";
    RangePair rangePair22 = RangePair(line22);
    rangePair22.getRangeQ().setNumChr(1);
    rangePair22.getRangeQ().setNameSeq("");
    rangePair22.getRangeS().setNumChr(1);
    rangePair22.getRangeS().setNameSeq("");
    rangePair22.setLength(108);
 
   // set match part to match
    std::list<RangePair> rpList2;
    rpList2.push_back(rangePair22);
    rangePairSet2.setPathDirectly(rpList2);
 
    std::list<RangePairSet> expRpsList;
    expRpsList.push_back(rangePairSet1);
    expRpsList.push_back(rangePairSet2);

    return expRpsList;
}


BLRMatchMap::MapPath Test_BLRMatchMapUtils::createMapPathForTest_extractRangePairSetListFromMapPath(void){
  BLRMatchMap::MapPath mapPath;

  // match  
  SDGString line1 = " \t105091\t105206\t \t3351\t3238\t0\t116\t76.1062";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(1);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(116);

  // match part
  SDGString line11 = " \t105091\t105206\t \t3351\t3238\t0\t116\t76.1062";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(1);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(116);
  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);
  rangePairSet1.setPathDirectly(rpList1);   
  
  // TODO insert in a private method
  std::list<RangePairSet>& rpsList = mapPath[BLRMatchMap::Key(rangePairSet1.getRangeQ().getNumChr(),rangePairSet1.getRangeS().getNumChr())];
  std::list<RangePairSet>::iterator r=std::lower_bound(rpsList.begin(), rpsList.end(),rangePairSet1);
  rpsList.insert(r, rangePairSet1);   
    
  // match    
  SDGString line2 = " \t109601\t109708\t \t3239\t3338\t3e-28\t108\t81.8182";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(108);
  // match part
  SDGString line22 = " \t109601\t109708\t \t3239\t3338\t3e-28\t108\t81.8182";
  RangePair rangePair22 = RangePair(line22);
  rangePair22.getRangeQ().setNumChr(1);
  rangePair22.getRangeQ().setNameSeq("");
  rangePair22.getRangeS().setNumChr(1);
  rangePair22.getRangeS().setNameSeq("");
  rangePair22.setLength(108);
  // set match part to match
  std::list<RangePair> rpList2;
  rpList2.push_back(rangePair22);
  rangePairSet2.setPathDirectly(rpList2);
  // insert
  rpsList = mapPath[BLRMatchMap::Key(rangePairSet2.getRangeQ().getNumChr(),rangePairSet2.getRangeS().getNumChr())];
  // TODO use rangePairSet1 as arg in lower_bound without understanding why ...
  r=std::lower_bound(rpsList.begin(), rpsList.end(),rangePairSet1);
  rpsList.insert(r, rangePairSet2);   
  return mapPath; 
}

std::list<RangePairSet> Test_BLRMatchMapUtils::createRpListForTest_add_split_path(void){
  return createRpList_for_test_add_clean_path_all_S(); 
}

std::list<RangePairSet> Test_BLRMatchMapUtils::createRpListForTest_isOverloapFound_in_add_split_path(void){
  BLRMatchMap::MapPath mapPath;
  
  SDGString line2 = " \t105085\t110398\t \t3202\t2033\t0\t695\t82.8467";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(695);
// match part   
  SDGString line21 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeS().setNameSeq("");
  rangePair21.setLength(414);
  // set match part to match
  std::list<RangePair> rpList2;
  rpList2.push_back(rangePair21);
  rangePairSet2.setPathDirectly(rpList2);

  // match part   
  SDGString line22 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair22 = RangePair(line22);
  rangePair22.getRangeQ().setNumChr(1);
  rangePair22.getRangeQ().setNameSeq("");
  rangePair22.getRangeS().setNumChr(1);
  rangePair22.getRangeS().setNameSeq("");
  rangePair22.setLength(281);
  // set match part to match
  rpList2.push_back(rangePair22);
  rangePairSet2.setPathDirectly(rpList2);
  // insert 
  std::list<RangePairSet>& rpsList = mapPath[BLRMatchMap::Key(rangePairSet2.getRangeQ().getNumChr(),rangePairSet2.getRangeS().getNumChr())];
  std::list<RangePairSet>::iterator r=std::lower_bound(rpsList.begin(), rpsList.end(),rangePairSet2);
  rpsList.insert(r, rangePairSet2);   

  return rpsList;

}

BLRMatchMap::MapPath Test_BLRMatchMapUtils::createMapPathForTest_isOverlapFound_in_add_split_path(void){
  BLRMatchMap::MapPath mapPath;
 
  // match 
  SDGString line3 = " \t109951\t110566\t \t119\t686\t9.4e-19\t202\t77.768";
  RangePairSet rangePairSet3 = RangePairSet(line3);
  rangePairSet3.getRangeQ().setNumChr(1);
  rangePairSet3.getRangeQ().setNameSeq("");
  rangePairSet3.getRangeS().setNumChr(2);
  rangePairSet3.getRangeS().setNameSeq("");
  rangePairSet3.setLength(638);
  // match part
  SDGString line31 = " \t109951\t109984\t \t119\t149\t9.4e-19\t34\t77.768";
  RangePair rangePair31 = RangePair(line31);
  rangePair31.getRangeQ().setNumChr(1);
  rangePair31.getRangeQ().setNameSeq("");
  rangePair31.getRangeS().setNumChr(2);
  rangePair31.getRangeS().setNameSeq("");
  rangePair31.setLength(34);
  // match part
  SDGString line32 = " \t110399\t110566\t \t532\t686\t9.4e-19\t168\t77.768";
  RangePair rangePair32 = RangePair(line32);
  rangePair32.getRangeQ().setNumChr(1);
  rangePair32.getRangeQ().setNameSeq("");
  rangePair32.getRangeS().setNumChr(2);
  rangePair32.getRangeS().setNameSeq("");
  rangePair32.setLength(604);
  // set match part to match
  std::list<RangePair> rpList3;
  rpList3.push_back(rangePair31);
  rpList3.push_back(rangePair32);
  rangePairSet3.setPathDirectly(rpList3);
  // insert
  std::list<RangePairSet>& rpsListForSecondKey = mapPath[BLRMatchMap::Key(rangePairSet3.getRangeQ().getNumChr(),rangePairSet3.getRangeS().getNumChr())];
  std::list<RangePairSet>::iterator rForSecondKey=std::lower_bound(rpsListForSecondKey.begin(), rpsListForSecondKey.end(),rangePairSet3);
  rpsListForSecondKey.insert(rForSecondKey, rangePairSet3);   
  return mapPath; 
}

BLRMatchMap::MapPath Test_BLRMatchMapUtils::createMapPath_afterJoin(void){
    bool joiningParameter = true;
    bool cleanBefore = false;
    bool cleanAfter = true;
    int verboseParameter = 0;
    
    SDGString match_file = "match.align";
    Test_BLRMatchMapUtils::writeInputFile();
    
    BLRMatcherThreadsParameter para = Test_BLRMatchMapUtils::createParameter(); 
    BLRMatchMap matchMap(&para);
     
    matchMap.loadAlign(0);
    BLRMatchMap::MapAlign mapAlignBefore = matchMap.getMapAlign();
    matchMap.mapPathJoinOnlyForTest(joiningParameter, cleanBefore, cleanAfter, verboseParameter);
    FileUtils::removeFile(match_file);
    return matchMap.getMapPath(); 
}

BLRMatchMap::MapPath Test_BLRMatchMapUtils::createExpMapPathForTest_mapPath_comparaison_diagnostic(void){
  BLRMatchMap::MapPath mapPath;

  // match  
  SDGString line1 = " \t105091\t105206\t \t3351\t3238\t0\t116\t76.1062";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(1);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(116);

  // match part
  SDGString line11 = " \t105091\t105206\t \t3351\t3238\t0\t116\t76.1062";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(1);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(116);
  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);
  rangePairSet1.setPathDirectly(rpList1);   
  
  // TODO insert in a private method
  std::list<RangePairSet>& rpsList = mapPath[BLRMatchMap::Key(rangePairSet1.getRangeQ().getNumChr(),rangePairSet1.getRangeS().getNumChr())];
  std::list<RangePairSet>::iterator r=std::lower_bound(rpsList.begin(), rpsList.end(),rangePairSet1);
  //rpsList.insert(r, rangePairSet1);   
  rpsList.push_back(rangePairSet1);  
  // match    
  SDGString line2 = " \t109601\t109708\t \t3239\t3338\t3e-28\t108\t81.8182";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(108);
  // match part
  SDGString line22 = " \t109601\t109708\t \t3239\t3338\t3e-28\t108\t81.8182";
  RangePair rangePair22 = RangePair(line22);
  rangePair22.getRangeQ().setNumChr(1);
  rangePair22.getRangeQ().setNameSeq("");
  rangePair22.getRangeS().setNumChr(1);
  rangePair22.getRangeS().setNameSeq("");
  rangePair22.setLength(108);
  // set match part to match
  std::list<RangePair> rpList2;
  rpList2.push_back(rangePair22);
  rangePairSet2.setPathDirectly(rpList2);
  // insert
  rpsList = mapPath[BLRMatchMap::Key(rangePairSet2.getRangeQ().getNumChr(),rangePairSet2.getRangeS().getNumChr())];
  // TODO use rangePairSet1 as arg in lower_bound without understanding why ...
  r=std::lower_bound(rpsList.begin(), rpsList.end(),rangePairSet2);
  //rpsList.insert(r, rangePairSet2);   
  rpsList.push_back(rangePairSet2);   

  // match 
  
  SDGString line3 = " \t105085\t110398\t \t3202\t2033\t0\t359\t80.6808";
  RangePairSet rangePairSet3 = RangePairSet(line3);
  rangePairSet3.getRangeQ().setNumChr(1);
  rangePairSet3.getRangeQ().setNameSeq("");
  rangePairSet3.getRangeS().setNumChr(1);
  rangePairSet3.getRangeS().setNameSeq("");
  rangePairSet3.setLength(695);
  // match part
  SDGString line31 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair31 = RangePair(line31);
  rangePair31.getRangeQ().setNumChr(1);
  rangePair31.getRangeQ().setNameSeq("");
  rangePair31.getRangeS().setNumChr(1);
  rangePair31.getRangeS().setNameSeq("");
  rangePair31.setLength(414);
  // match part
  SDGString line32 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair32 = RangePair(line32);
  rangePair32.getRangeQ().setNumChr(1);
  rangePair32.getRangeQ().setNameSeq("");
  rangePair32.getRangeS().setNumChr(1);
  rangePair32.getRangeS().setNameSeq("");
  rangePair32.setLength(281);
  // set match part to match
  std::list<RangePair> rpList3;
  rpList3.push_back(rangePair31);
  rpList3.push_back(rangePair32);
  rangePairSet3.setPathDirectly(rpList3);
  // insert
  rpsList = mapPath[BLRMatchMap::Key(rangePairSet3.getRangeQ().getNumChr(),rangePairSet3.getRangeS().getNumChr())];
  r=std::lower_bound(rpsList.begin(), rpsList.end(),rangePairSet3);
  //rpsList.insert(r, rangePairSet3);   
  rpsList.push_back(rangePairSet3);   
  
  
  // match    
  SDGString line4 = " \t109527\t109818\t \t3005\t3316\t0\t292\t82.8467";
  RangePairSet rangePairSet4 = RangePairSet(line4);
  rangePairSet4.getRangeQ().setNumChr(1);
  rangePairSet4.getRangeQ().setNameSeq("");
  rangePairSet4.getRangeS().setNumChr(1);
  rangePairSet4.getRangeS().setNameSeq("");
  rangePairSet4.setLength(311);
  // match part
  SDGString line44 = " \t109527\t109818\t \t3005\t3316\t0\t292\t82.8467";
  RangePair rangePair44 = RangePair(line44);
  rangePair44.getRangeQ().setNumChr(1);
  rangePair44.getRangeQ().setNameSeq("");
  rangePair44.getRangeS().setNumChr(1);
  rangePair44.getRangeS().setNameSeq("");
  rangePair44.setLength(311);
  // set match part to match
  std::list<RangePair> rpList4;
  rpList4.push_back(rangePair44);
  rangePairSet4.setPathDirectly(rpList4);
  // insert
  rpsList = mapPath[BLRMatchMap::Key(rangePairSet4.getRangeQ().getNumChr(),rangePairSet4.getRangeS().getNumChr())];
  // TODO use rangePairSet3 as arg in lower_bound without understanding why ...
  r=std::lower_bound(rpsList.begin(), rpsList.end(),rangePairSet3);
  //rpsList.insert(r, rangePairSet4);   
  rpsList.push_back(rangePairSet4);   
 
  // match
  SDGString line5 = " \t109951\t110569\t \t119\t689\t9.4e-19\t619\t77.768";
  RangePairSet rangePairSet5 = RangePairSet(line5);
  rangePairSet5.getRangeQ().setNumChr(1);
  rangePairSet5.getRangeQ().setNameSeq("");
  rangePairSet5.getRangeS().setNumChr(2);
  rangePairSet5.getRangeS().setNameSeq("");
  rangePairSet5.setLength(619);
  // match part
  SDGString line51 = " \t109951\t110569\t \t119\t689\t9.4e-19\t619\t77.768";
  RangePairSet rangePair51 = RangePairSet(line51);
  rangePair51.getRangeQ().setNumChr(1);
  rangePair51.getRangeQ().setNameSeq("");
  rangePair51.getRangeS().setNumChr(2);
  rangePair51.getRangeS().setNameSeq("");
  rangePair51.setLength(619);
  // set match part to match
  std::list<RangePair> rpList5;
  rpList5.push_back(rangePair51);
  rangePairSet5.setPathDirectly(rpList5);
  // insert
  std::list<RangePairSet>& rpsListForSecondKey = mapPath[BLRMatchMap::Key(rangePairSet5.getRangeQ().getNumChr(),rangePairSet5.getRangeS().getNumChr())];
  std::list<RangePairSet>::iterator rForSecondKey=std::lower_bound(rpsListForSecondKey.begin(), rpsListForSecondKey.end(),rangePairSet5);
  rpsListForSecondKey.insert(rForSecondKey, rangePairSet5);   
   
  // match
  SDGString line6 = " \t110567\t110878\t \t2833\t2532\t0\t312\t77.2076";
  RangePairSet rangePairSet6 = RangePairSet(line6);
  rangePairSet6.getRangeQ().setNumChr(1);
  rangePairSet6.getRangeQ().setNameSeq("");
  rangePairSet6.getRangeS().setNumChr(3);
  rangePairSet6.getRangeS().setNameSeq("");
  rangePairSet6.setLength(312);
  // match part
  SDGString line66 = " \t110567\t110878\t \t2833\t2532\t0\t312\t77.2076";
  RangePairSet rangePair66 = RangePairSet(line66);
  rangePair66.getRangeQ().setNumChr(1);
  rangePair66.getRangeQ().setNameSeq("");
  rangePair66.getRangeS().setNumChr(3);
  rangePair66.getRangeS().setNameSeq("");
  rangePair66.setLength(312);
  // set match part to match
  std::list<RangePair> rpList6;
  rpList6.push_back(rangePair66);
  rangePairSet6.setPathDirectly(rpList6);
  // insert
  std::list<RangePairSet>& rpsListForThirdKey = mapPath[BLRMatchMap::Key(rangePairSet6.getRangeQ().getNumChr(),rangePairSet6.getRangeS().getNumChr())];
  std::list<RangePairSet>::iterator rForThirdKey=std::lower_bound(rpsListForThirdKey.begin(), rpsListForThirdKey.end(),rangePairSet6);
  rpsListForThirdKey.insert(rForThirdKey, rangePairSet6);   
  return mapPath; 

}

std::list<RangePairSet> Test_BLRMatchMapUtils::createInputRpsList_for_test_reComputeScoreWithLength_on_a_match_with_one_match_part(void){
  // match  
  SDGString line1 = " \t105091\t105206\t \t3351\t3238\t0\t11\t76.1062";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(1);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(116);
  // match part
  SDGString line11 = " \t105091\t105206\t \t3351\t3238\t0\t11\t76.1062";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(1);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(116);
  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);
  rangePairSet1.setPathDirectly(rpList1);   
  std::list<RangePairSet> rpsList;
  rpsList.push_back(rangePairSet1);
  return rpsList; 
 
}

std::list<RangePairSet> Test_BLRMatchMapUtils::createExpRpsList_for_test_reComputeScoreWithLength_on_a_match_with_one_match_part(void){
  // match  
  SDGString line1 = " \t105091\t105206\t \t3351\t3238\t0\t116\t76.1062";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(1);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(116);
  // match part
  SDGString line11 = " \t105091\t105206\t \t3351\t3238\t0\t116\t76.1062";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(1);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(116);
  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);
  rangePairSet1.setPathDirectly(rpList1);   
  std::list<RangePairSet> rpsList;
  rpsList.push_back(rangePairSet1);
  return rpsList; 
 
}

std::list<RangePairSet> Test_BLRMatchMapUtils::createInputRpsList_for_test_reComputeScoreWithLength_on_a_match_with_two_match_part(void){
  // match 
  SDGString line1 = " \t105085\t110398\t \t3202\t2033\t0\t35\t80.6808";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(1);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(695);
  // match part
  SDGString line11 = " \t109985\t110398\t \t2429\t2033\t0\t41\t77.8604";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(1);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(414);
  // match part
  SDGString line12 = " \t105085\t105353\t \t3202\t2921\t0\t26\t84.8361";
  RangePair rangePair12 = RangePair(line12);
  rangePair12.getRangeQ().setNumChr(1);
  rangePair12.getRangeQ().setNameSeq("");
  rangePair12.getRangeS().setNumChr(1);
  rangePair12.getRangeS().setNameSeq("");
  rangePair12.setLength(281);
  // set match part to match
  std::list<RangePair> rpList;
  rpList.push_back(rangePair11);
  rpList.push_back(rangePair12);
  rangePairSet1.setPathDirectly(rpList);
  // insert
  std::list<RangePairSet> rpsList;
  rpsList.push_back(rangePairSet1);   
  return rpsList;
  
} 

std::list<RangePairSet> Test_BLRMatchMapUtils::createExpRpsList_for_test_reComputeScoreWithLength_on_a_match_with_two_match_part(void){
  // match 
  SDGString line1 = " \t105085\t110398\t \t3202\t2033\t0\t683\t80.6808";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(1);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(683);
  // match part
  SDGString line11 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(1);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(414);
  // match part
  SDGString line12 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair12 = RangePair(line12);
  rangePair12.getRangeQ().setNumChr(1);
  rangePair12.getRangeQ().setNameSeq("");
  rangePair12.getRangeS().setNumChr(1);
  rangePair12.getRangeS().setNameSeq("");
  rangePair12.setLength(269);
  // set match part to match
  std::list<RangePair> rpList;
  rpList.push_back(rangePair11);
  rpList.push_back(rangePair12);
  rangePairSet1.setPathDirectly(rpList);
  // insert
  std::list<RangePairSet> rpsList;
  rpsList.push_back(rangePairSet1);   
  return rpsList;
}


std::list<RangePairSet> Test_BLRMatchMapUtils::createExpRpsListForOnlyJoinComputeScore(void)
{
  BLRMatchMap::MapPath mapPath;

  // match 
  SDGString line3 = " \t105085\t110398\t \t3202\t2033\t0\t683\t80.6808";
  RangePairSet rangePairSet3 = RangePairSet(line3);
  rangePairSet3.getRangeQ().setNumChr(1);
  rangePairSet3.getRangeQ().setNameSeq("");
  rangePairSet3.getRangeS().setNumChr(1);
  rangePairSet3.getRangeS().setNameSeq("");
  rangePairSet3.setLength(683);
  // match part
  SDGString line31 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair31 = RangePair(line31);
  rangePair31.getRangeQ().setNumChr(1);
  rangePair31.getRangeQ().setNameSeq("");
  rangePair31.getRangeS().setNumChr(1);
  rangePair31.getRangeS().setNameSeq("");
  rangePair31.setLength(414);
  // match part
  SDGString line32 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair32 = RangePair(line32);
  rangePair32.getRangeQ().setNumChr(1);
  rangePair32.getRangeQ().setNameSeq("");
  rangePair32.getRangeS().setNumChr(1);
  rangePair32.getRangeS().setNameSeq("");
  rangePair32.setLength(269);
  // set match part to match
  std::list<RangePair> rpList3;
  rpList3.push_back(rangePair31);
  rpList3.push_back(rangePair32);
  rangePairSet3.setPathDirectly(rpList3);
  // insert
  std::list<RangePairSet> rpsList;
  rpsList.push_back(rangePairSet3);   
  return rpsList;

}

bool Test_BLRMatchMapUtils::areTwoRangePairSetEqualsWithoutIdentity(RangePairSet rps1, RangePairSet rps2)
{
	return( rps1.first == rps2.first
			&& rps1.second == rps2.second
			&& rps1.getE_value() == rps2.getE_value()
			&& rps1.getScore() == rps2.getScore()
			&& rps1.getPath() == rps2.getPath());

}

std::list<RangePairSet> Test_BLRMatchMapUtils::createInputRpsListFor_test_add_split_path(void)
{
    bool joiningParameter = true;
    bool cleanBefore = false;
    bool cleanAfter = true;
    int verboseParameter = 0;
 
    SDGString match_file = "match.align";
    Test_BLRMatchMapUtils::writeInputFile();
    
    BLRMatcherThreadsParameter para = Test_BLRMatchMapUtils::createParameter(); 
    BLRMatchMap matchMap(&para);
     
    matchMap.loadAlign(0);
    BLRMatchMap::MapAlign mapAlignBefore = matchMap.getMapAlign();
    matchMap.mapPathJoinAndComputeScoreWithLengthOnly(joiningParameter, cleanBefore, cleanAfter, verboseParameter);
    return matchMap.getRpsList();
}

// TODO: Refactor according comments ...
BLRMatchMap::MapPath Test_BLRMatchMapUtils::createMapPath(std::list<SDGString> inputList)
{
  	BLRMatchMap::MapPath mapPath;
	std::list<RangePairSet> rpsContainer;	
	for(std::list<SDGString>::iterator input_it=inputList.begin(); input_it!=inputList.end();input_it++){
	
		size_t posFirstTab = input_it->find_first_of("\t");
		size_t posLastTab = input_it->find_last_of("\t");
		size_t posEndLine = input_it->find_first_of("\\n");
		size_t lengthSeq = posEndLine - posLastTab;
		size_t lengthInit = posLastTab - posFirstTab;
		
		// isARangePairSetToCreate		
		if (posFirstTab != 0){
			// createARangePairSet 
			RangePairSet currentRps = RangePairSet(input_it->substr(0, lengthInit));
			SDGString sNameSeq = currentRps.getRangeS().getNameSeq();
			SDGString qNameSeq = currentRps.getRangeQ().getNameSeq();
			currentRps.getRangeS().setNumChr(atol(sNameSeq));
			currentRps.getRangeQ().setNumChr(atol(qNameSeq));
		        
	
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
	
	// storeContainerInMapPath
	for(std::list<RangePairSet>::iterator lrp_it=rpsContainer.begin(); lrp_it!=rpsContainer.end();lrp_it++){
		std::list<RangePairSet> mapRpsList = mapPath[BLRMatchMap::Key(lrp_it->getRangeQ().getNumChr(), lrp_it->getRangeS().getNumChr())];
		RangePairSet rps = *lrp_it;
		mapRpsList.push_back(rps); 
		mapPath[BLRMatchMap::Key(lrp_it->getRangeQ().getNumChr(), lrp_it->getRangeS().getNumChr())] = mapRpsList;
    	}

	return mapPath;
}
//----------------------------------------------------------------------------
bool Test_BLRMatchMapUtils::isOverlapFound_in_add_split_path(std::list<RangePairSet>::iterator iter, BLRMatchMap::MapPath mapPath,
                                                   double idTolerance, unsigned lenFilter) {

  bool found_over = false;
  std::list<RangePairSet> lrp;

  lrp.push_back(*iter);
  for (BLRMatchMap::MapPath::iterator m = mapPath.begin(); m != mapPath.end(); m++) {

    if (m->first.first == iter->getRangeQ().getNumChr()) {

      for (std::list<RangePairSet>::iterator lrp_it = lrp.begin();
           lrp_it != lrp.end();
           lrp_it++)
        for (std::list<RangePairSet>::iterator iter_list = m->second.begin();
             iter_list != m->second.end(); iter_list++)
          if (lrp_it->overlapQ(*iter_list)) {
            if (lrp_it->getLength() >= 100 &&
                lrp_it->inserted(*iter_list) &&
                fabs(iter_list->getIdentity()
                     - lrp_it->getIdentity()) <= idTolerance)
              continue;
            std::list<RangePairSet> lrp2;

            if (lrp_it->split(*iter_list, lrp2)) {
              found_over = true;
              for (std::list<RangePairSet>::iterator lrp2_it
                      = lrp2.begin(); lrp2_it != lrp2.end(); lrp2_it++)
                if (lrp2_it->getRangeQ().getLength()
                    > lenFilter)
                  lrp.push_back(*lrp2_it);
            }

          }
    }
  }
  return found_over;
}