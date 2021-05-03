//
// Created by Hadi Quesneville on 25/03/2020.
//

#include "Test_BLRUtils.h"
#include <stack>
#include <sstream>
#include "SDGString.h"
#include "FileUtils.h"
#include "BLRMatchMap.h"
#include "RangePair.h"
#include "RangePairSet.h"

void Test_BLRUtils::writeInputFile(void)
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

BLRJoinParameter Test_BLRUtils::createParameter(void) {
    BLRJoinParameter para;
    para.setLenFilter(0);
    para.setEvalFilter(10);
    para.setIdFilter(0);
    return para;
};

BLRMatchMap::MapAlign Test_BLRUtils::createExpMapAlign_for_clean_conflicts(void){
    SDGString match_file = "match.align";
    writeInputFile();
    BLRJoinParameter para = createParameter();
    BLRMatchMap matchMap(para);
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

std::list<RangePairSet> Test_BLRUtils::createRpList_for_test_add_clean_path_all_S(void){
    BLRJoinParameter para = createParameter();
    BLRMatchMap matchMap(para);
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

BLRMatchMap::MapAlign Test_BLRUtils::createMapAlign_instance_for_test_mapAlign_Equality(void){
    SDGString match_file = "match.align";
    writeInputFile();
    BLRJoinParameter para = createParameter();
    BLRMatchMap matchMap(para);
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

std::list<RangePairSet> Test_BLRUtils::createRpListForTest_isOverloapFound_in_add_split_path(void){
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

BLRMatchMap::MapPath Test_BLRUtils::createMapPathForTest_isOverlapFound_in_add_split_path(void){
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

std::list<RangePairSet> Test_BLRUtils::createInputRpsList_for_test_reComputeScoreWithLength_on_a_match_with_one_match_part(void){
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

std::list<RangePairSet> Test_BLRUtils::createExpRpsList_for_test_reComputeScoreWithLength_on_a_match_with_one_match_part(void){
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

std::list<RangePairSet> Test_BLRUtils::createInputRpsList_for_test_reComputeScoreWithLength_on_a_match_with_two_match_part(void){
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

std::list<RangePairSet> Test_BLRUtils::createExpRpsList_for_test_reComputeScoreWithLength_on_a_match_with_two_match_part(void){
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


//----------------------------------------------------------------------------
bool Test_BLRUtils::isOverlapFound_in_add_split_path(std::list<RangePairSet>::iterator iter, BLRMatchMap::MapPath mapPath,
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