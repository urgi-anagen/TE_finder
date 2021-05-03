#include "Test_RangePairSetUtils.h"
#include "SDGString.h"
#include "RangePairSet.h"

void Test_RangePairSetUtils::viewRangePairSetList(std::list<RangePairSet> rpsList){
    for(std::list<RangePairSet>::iterator lrp_it=rpsList.begin(); lrp_it!=rpsList.end();lrp_it++){
            lrp_it->view(); 
    }
}

RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_inserted(void)
{

  std::list<RangePair> rpList2;
  SDGString line2 = " \t105085\t110398\t \t3202\t2033\t0\t695\t82.8467";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(695);
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

  // match part   
  SDGString line21 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeS().setNameSeq("");
  rangePair21.setLength(414);
  rpList2.push_back(rangePair21);

  rangePairSet2.setPathDirectly(rpList2);

  return rangePairSet2;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_inserted(void)
{
// match 
  SDGString line1 = " \t109951\t110566\t \t119\t686\t9.4e-19\t202\t77.768";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(2);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(638);
  // match part
  SDGString line11 = " \t109951\t109984\t \t119\t149\t9.4e-19\t34\t77.768";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(2);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(34);

  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);

  // match part
  SDGString line12 = " \t110399\t110566\t \t532\t686\t9.4e-19\t168\t77.768";
  RangePair rangePair12 = RangePair(line12);
  rangePair12.getRangeQ().setNumChr(1);
  rangePair12.getRangeQ().setNameSeq("");
  rangePair12.getRangeS().setNumChr(2);
  rangePair12.getRangeS().setNameSeq("");
  rangePair12.setLength(604);

  rpList1.push_back(rangePair12);

  rangePairSet1.setPathDirectly(rpList1);

  return rangePairSet1;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_inserted_true(void)
{
// match 
  SDGString line1 = " \t100\t1000\t \t119\t686\t9.4e-19\t202\t77.768";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(2);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(638);
  // match part
  SDGString line11 = " \t100\t150\t \t119\t149\t9.4e-19\t34\t77.768";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(2);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(34);

  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);

  // match part
  SDGString line12 = " \t900\t1000\t \t532\t686\t9.4e-19\t168\t77.768";
  RangePair rangePair12 = RangePair(line12);
  rangePair12.getRangeQ().setNumChr(1);
  rangePair12.getRangeQ().setNameSeq("");
  rangePair12.getRangeS().setNumChr(2);
  rangePair12.getRangeS().setNameSeq("");
  rangePair12.setLength(604);

  rpList1.push_back(rangePair12);

  rangePairSet1.setPathDirectly(rpList1);


  return rangePairSet1;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_inserted_true(void)
{

  std::list<RangePair> rpList2;
  SDGString line2 = " \t300\t500\t \t3202\t2033\t0\t695\t82.8467";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(695);
  
  // match part   
  SDGString line22 = " \t300\t350\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair22 = RangePair(line22);
  rangePair22.getRangeQ().setNumChr(1);
  rangePair22.getRangeQ().setNameSeq("");
  rangePair22.getRangeS().setNumChr(1);
  rangePair22.getRangeS().setNameSeq("");
  rangePair22.setLength(100);
  // set match part to match
  rpList2.push_back(rangePair22);


// match part   
  SDGString line21 = " \t400\t500\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeS().setNameSeq("");
  rangePair21.setLength(100);
  rpList2.push_back(rangePair21);

  rangePairSet2.setPathDirectly(rpList2);

  return rangePairSet2;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_overlapQ(void)
{
    // match 
  SDGString line1 = " \t109951\t110569\t \t119\t686\t9.4e-19\t202\t77.768";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(2);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(638);
  // match part
  SDGString line11 = " \t109951\t110569\t \t119\t686\t9.4e-19\t202\t77.768";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(2);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(638);

  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);

  rangePairSet1.setPathDirectly(rpList1);

  return rangePairSet1;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_overlapQ(void)
{
  std::list<RangePair> rpList2;
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
  rpList2.push_back(rangePair21);
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

  return rangePairSet2;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_overlapQ_length(void)
{
    return createRangePairSet1ForTest_overlapQ();
}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_overlapQ_length(void)
{
    return createRangePairSet2ForTest_overlapQ();
}

RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_overlapQ_length_2(void)
{
	
	SDGString str1 = "1\t1\t700\t1\t3351\t3238\t0\t431\t76.1062\t431\n";
	SDGString str11 = "\t1\t1\t100\t1\t3351\t3238\t0\t100\t76.1062\t100\n";
	SDGString str12 = "\t1\t370\t700\t1\t3351\t3238\t0\t231\t76.1062\t231\n";
	
	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);		
	lPath.push_back(str12);		

	return  createRangePairSet(lPath);	

}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_overlapQ_length_2(void)
{
	SDGString str2 = "1\t150\t380\t2\t119\t689\t9.4e-19\t231\t77.768\t231\n";
	SDGString str21 = "\t1\t150\t380\t2\t119\t689\t9.4e-19\t231\t77.768\t231\n";
	
	std::list<SDGString> lPath;
	
	lPath.push_back(str2);		
	lPath.push_back(str21);		

	return  createRangePairSet(lPath);	

}
RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_diffQ(void)
{
     // match 
  SDGString line1 = " \t109951\t110569\t \t119\t689\t9.4e-19\t619\t77.768";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(2);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(619);
  // match part
  SDGString line11 = " \t109951\t110569\t \t119\t689\t9.4e-19\t619\t77.768";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(2);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(619);

  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);

  rangePairSet1.setPathDirectly(rpList1);

  return rangePairSet1;
}   

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_diffQ(void)
{
  std::list<RangePair> rpList2;
  SDGString line2 = " \t105085\t110398\t \t3202\t2033\t0\t683\t80.6808";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(683);
  // match part   
  SDGString line21 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeS().setNameSeq("");
  rangePair21.setLength(414);
  rpList2.push_back(rangePair21);
  // match part   
  SDGString line22 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair22 = RangePair(line22);
  rangePair22.getRangeQ().setNumChr(1);
  rangePair22.getRangeQ().setNameSeq("");
  rangePair22.getRangeS().setNumChr(1);
  rangePair22.getRangeS().setNameSeq("");
  rangePair22.setLength(269);
  // set match part to match
  rpList2.push_back(rangePair22);

  rangePairSet2.setPathDirectly(rpList2);

  return rangePairSet2;
}

RangePairSet Test_RangePairSetUtils::createExpRangePairSet1ForTest_diffQ(void)
{
// match 
  SDGString line1 = " \t109951\t110569\t \t119\t689\t9.4e-19\t205\t77.768";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(2);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(619);
  // match part
  SDGString line11 = " \t109951\t109984\t \t119\t149\t9.4e-19\t34\t77.768";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(2);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(34);

  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);

  // match part
  SDGString line12 = " \t110399\t110569\t \t532\t689\t9.4e-19\t171\t77.768";
  RangePair rangePair12 = RangePair(line12);
  rangePair12.getRangeQ().setNumChr(1);
  rangePair12.getRangeQ().setNameSeq("");
  rangePair12.getRangeS().setNumChr(2);
  rangePair12.getRangeS().setNameSeq("");
  rangePair12.setLength(604);

  rpList1.push_back(rangePair12);

  rangePairSet1.setPathDirectly(rpList1);

  return rangePairSet1;
}

RangePairSet Test_RangePairSetUtils::createExpRangePairSet2ForTest_diffQ(void)
{
    return createRangePairSet2ForTest_diffQ();
}

RangePairSet Test_RangePairSetUtils::createRangePairSet1For_test_equality(void)
{
  std::list<RangePair> rpList1;
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
  rpList1.push_back(rangePair11);
  // match part   
  SDGString line12 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair12 = RangePair(line12);
  rangePair12.getRangeQ().setNumChr(1);
  rangePair12.getRangeQ().setNameSeq("");
  rangePair12.getRangeS().setNumChr(1);
  rangePair12.getRangeS().setNameSeq("");
  rangePair12.setLength(269);
  // set match part to match
  rpList1.push_back(rangePair12);

  rangePairSet1.setPathDirectly(rpList1);

  return rangePairSet1;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet2For_test_equality(void)
{
  std::list<RangePair> rpList2;
  SDGString line2 = " \t105085\t110398\t \t3202\t2033\t0\t683\t80.6808";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(683);
  // match part   
  SDGString line21 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeS().setNameSeq("");
  rangePair21.setLength(414);
  rpList2.push_back(rangePair21);
  // match part   
  SDGString line22 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair22 = RangePair(line22);
  rangePair22.getRangeQ().setNumChr(1);
  rangePair22.getRangeQ().setNameSeq("");
  rangePair22.getRangeS().setNumChr(1);
  rangePair22.getRangeS().setNameSeq("");
  rangePair22.setLength(269);
  // set match part to match
  rpList2.push_back(rangePair22);

  rangePairSet2.setPathDirectly(rpList2);

  return rangePairSet2;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet1For_test_equality_not_equal(void)
{
	return createRangePairSet1For_test_equality();
}

RangePairSet Test_RangePairSetUtils::createRangePairSet2For_test_equality_not_equal(void)
{
  std::list<RangePair> rpList2;
  SDGString line2 = " \t105084\t110398\t \t3202\t2033\t0\t683\t80.6808";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(683);
  // match part   
  SDGString line21 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeS().setNameSeq("");
  rangePair21.setLength(414);
  rpList2.push_back(rangePair21);
  // match part   
  SDGString line22 = " \t105084\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair22 = RangePair(line22);
  rangePair22.getRangeQ().setNumChr(1);
  rangePair22.getRangeQ().setNameSeq("");
  rangePair22.getRangeS().setNumChr(1);
  rangePair22.getRangeS().setNameSeq("");
  rangePair22.setLength(269);
  // set match part to match
  rpList2.push_back(rangePair22);

  rangePairSet2.setPathDirectly(rpList2);

  return rangePairSet2;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet1For_test_equality_on_match_part(void)
{
  std::list<RangePair> rpList1;
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
  rpList1.push_back(rangePair11);
  // match part   
  SDGString line12 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair12 = RangePair(line12);
  rangePair12.getRangeQ().setNumChr(1);
  rangePair12.getRangeQ().setNameSeq("");
  rangePair12.getRangeS().setNumChr(1);
  rangePair12.getRangeS().setNameSeq("");
  rangePair12.setLength(269);
  // set match part to match
  rpList1.push_back(rangePair12);

  rangePairSet1.setPathDirectly(rpList1);

  return rangePairSet1;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet2For_test_equality_on_match_part(void)
{
  std::list<RangePair> rpList2;
  SDGString line2 = " \t105085\t110398\t \t3202\t2033\t0\t683\t80.6808";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(683);
  // match part   
  SDGString line21 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeS().setNameSeq("");
  rangePair21.setLength(414);
  rpList2.push_back(rangePair21);
  // match part   
  SDGString line22 = " \t105084\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair22 = RangePair(line22);
  rangePair22.getRangeQ().setNumChr(1);
  rangePair22.getRangeQ().setNameSeq("");
  rangePair22.getRangeS().setNumChr(1);
  rangePair22.getRangeS().setNameSeq("");
  rangePair22.setLength(269);
  // set match part to match
  rpList2.push_back(rangePair22);

  rangePairSet2.setPathDirectly(rpList2);

  return rangePairSet2;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet1For_test_equality_rps1_rps2_different_length(void)
{
	SDGString str1 = "1\t1\t900\t-1\t0\t0\t0\t900\t0\t0\t900\n";
	SDGString str11 = "\t1\t1\t99\t2\t181\t119\t9.4e-19\t99\t77.768\t99\n";
	SDGString str12 = "\t1\t500\t100\t1\t3351\t3238\t0\t5\t76.1062\t401\n";
	SDGString str13 = "\t1\t501\t599\t2\t499\t437\t9.4e-19\t99\t77.768\t222\n";
	SDGString str14 = "\t1\t850\t600\t1\t3351\t3238\t0\t5\t76.1062\t251\n";
	SDGString str15 = "\t1\t851\t900\t2\t689\t659\t9.4e-19\t50\t77.768\t50\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);	
	lPath.push_back(str12);	
	lPath.push_back(str13);	
	lPath.push_back(str14);	
	lPath.push_back(str15);	

	return createRangePairSet(lPath);

}

RangePairSet Test_RangePairSetUtils::createRangePairSet2For_test_equality_rps1_rps2_different_length(void)
{
	SDGString str1 = "1\t1\t900\t-1\t0\t0\t0\t900\t0\t0\t900\n";
	SDGString str11 = "\t1\t1\t99\t2\t181\t119\t9.4e-19\t99\t77.768\t98\n";
	SDGString str12 = "\t1\t500\t100\t1\t3351\t3238\t0\t5\t76.1062\t400\n";
	SDGString str13 = "\t1\t501\t599\t2\t499\t437\t9.4e-19\t99\t77.768\t221\n";
	SDGString str14 = "\t1\t850\t600\t1\t3351\t3238\t0\t5\t76.1062\t250\n";
	SDGString str15 = "\t1\t851\t900\t2\t689\t659\t9.4e-19\t50\t77.768\t49\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);	
	lPath.push_back(str12);	
	lPath.push_back(str13);	
	lPath.push_back(str14);	
	lPath.push_back(str15);	

	return createRangePairSet(lPath);

}


RangePairSet Test_RangePairSetUtils::createInputRangePairSetFor_test_computeScoreWithLength_one_match_part(void)
{
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
  return rangePairSet1;
}

RangePairSet Test_RangePairSetUtils::createExpRangePairSetFor_test_computeScoreWithLength_one_match_part(void)
{
  // match  
  SDGString line1 = " \t105091\t105206\t \t3351\t3238\t0\t8828\t76.1062";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(1);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(116);
  // match part
  SDGString line11 = " \t105091\t105206\t \t3351\t3238\t0\t8828\t76.1062";
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
  return rangePairSet1; 
}

RangePairSet Test_RangePairSetUtils::createInputRangePairSetFor_test_computeScoreWithLength_two_match_part(void)
{
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
  return rangePairSet1;
}

RangePairSet Test_RangePairSetUtils::createExpRangePairSetFor_test_computeScoreWithLength_two_match_part(void)
{
  // match 
  SDGString line1 = " \t105085\t110398\t \t3202\t2033\t0\t55055\t80.6808";
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
  return rangePairSet1;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_overalpQ_length(void)
{
  std::list<RangePair> rpList2;
  SDGString line2 = " \t105085\t110398\t \t3202\t2033\t0\t683\t80.6808";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(683);
  // match part   
  SDGString line21 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeS().setNameSeq("");
  rangePair21.setLength(414);
  rpList2.push_back(rangePair21);
  // match part   
  SDGString line22 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair22 = RangePair(line22);
  rangePair22.getRangeQ().setNumChr(1);
  rangePair22.getRangeQ().setNameSeq("");
  rangePair22.getRangeS().setNumChr(1);
  rangePair22.getRangeS().setNameSeq("");
  rangePair22.setLength(269);
  // set match part to match
  rpList2.push_back(rangePair22);

  rangePairSet2.setPathDirectly(rpList2);

  return rangePairSet2;

}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_overalpQ_length(void)
{
     // match 
  SDGString line1 = " \t109951\t110569\t \t119\t689\t9.4e-19\t619\t77.768";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(2);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(619);
  // match part
  SDGString line11 = " \t109951\t110569\t \t119\t689\t9.4e-19\t619\t77.768";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(2);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(619);

  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);

  rangePairSet1.setPathDirectly(rpList1);

  return rangePairSet1;

}

RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_mergeQ_score_rps1_greater_score_than_rps2(void)
{
  std::list<RangePair> rpList2;
  SDGString line2 = " \t105085\t110398\t \t3202\t2033\t0\t683\t80.6808";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(683);
  // match part   
  SDGString line21 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeS().setNameSeq("");
  rangePair21.setLength(414);
  rpList2.push_back(rangePair21);
  // match part   
  SDGString line22 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair22 = RangePair(line22);
  rangePair22.getRangeQ().setNumChr(1);
  rangePair22.getRangeQ().setNameSeq("");
  rangePair22.getRangeS().setNumChr(1);
  rangePair22.getRangeS().setNameSeq("");
  rangePair22.setLength(269);
  // set match part to match
  rpList2.push_back(rangePair22);

  rangePairSet2.setPathDirectly(rpList2);

  return rangePairSet2;


}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_mergeQ_score_rps1_greater_score_than_rps2(void)
{
     // match 
  SDGString line1 = " \t109951\t110569\t \t119\t689\t9.4e-19\t619\t77.768";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(2);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(619);
  // match part
  SDGString line11 = " \t109951\t110569\t \t119\t689\t9.4e-19\t619\t77.768";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(2);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(619);

  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);

  rangePairSet1.setPathDirectly(rpList1);

  return rangePairSet1;
}

RangePairSet Test_RangePairSetUtils::createExpRangePairSetForTest_mergeQ_score_rps1_greater_score_than_rps2(void)
{
  std::list<RangePair> rpList;
  SDGString line1 = " \t105085\t110569\t \t?\t?\t?\t?\t?";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(-1);
  rangePairSet1.getRangeS().setNameSeq("-1");
  rangePairSet1.setLength(888);
  rangePairSet1.setScore(888);

  // match part   
  SDGString line11 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(1);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(269);
  rpList.push_back(rangePair11);
  
  // match part
  SDGString line12 = " \t109951\t109984\t \t119\t149\t9.4e-19\t34\t77.768";
  RangePair rangePair12 = RangePair(line12);
  rangePair12.getRangeQ().setNumChr(1);
  rangePair12.getRangeQ().setNameSeq("");
  rangePair12.getRangeS().setNumChr(2);
  rangePair12.getRangeS().setNameSeq("");
  rangePair12.setLength(34);
  rpList.push_back(rangePair12);
 
  // match part   
  SDGString line13 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair13 = RangePair(line13);
  rangePair13.getRangeQ().setNumChr(1);
  rangePair13.getRangeQ().setNameSeq("");
  rangePair13.getRangeS().setNumChr(1);
  rangePair13.getRangeS().setNameSeq("");
  rangePair13.setLength(414);
  rpList.push_back(rangePair13);
 
  // match part
  SDGString line14 = " \t110399\t110569\t \t532\t689\t9.4e-19\t171\t77.768";
  RangePair rangePair14 = RangePair(line14);
  rangePair14.getRangeQ().setNumChr(1);
  rangePair14.getRangeQ().setNameSeq("");
  rangePair14.getRangeS().setNumChr(2);
  rangePair14.getRangeS().setNameSeq("");
  rangePair14.setLength(615);
  rpList.push_back(rangePair14);

  rangePairSet1.setPathDirectly(rpList);
  return rangePairSet1;

}

RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_mergeQ_score_rps1_lower_score_than_rps2(void)
{
  std::list<RangePair> rpList2;
  SDGString line2 = " \t105085\t110398\t \t3202\t2033\t0\t683\t80.6808";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(683);
  // match part   
  SDGString line21 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeS().setNameSeq("");
  rangePair21.setLength(414);
  rpList2.push_back(rangePair21);
  // match part   
  SDGString line22 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair22 = RangePair(line22);
  rangePair22.getRangeQ().setNumChr(1);
  rangePair22.getRangeQ().setNameSeq("");
  rangePair22.getRangeS().setNumChr(1);
  rangePair22.getRangeS().setNameSeq("");
  rangePair22.setLength(269);
  // set match part to match
  rpList2.push_back(rangePair22);

  rangePairSet2.setPathDirectly(rpList2);

  return rangePairSet2;

}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_mergeQ_score_rps1_lower_score_than_rps2(void)
{
     // match 
  SDGString line1 = " \t109951\t110650\t \t119\t700\t9.4e-19\t700\t77.768";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(2);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(700);
  // match part
  SDGString line11 = " \t109951\t110650\t \t119\t700\t9.4e-19\t700\t77.768";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(2);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(700);

  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);

  rangePairSet1.setPathDirectly(rpList1);

  return rangePairSet1;

}

RangePairSet Test_RangePairSetUtils::createExpRangePairSetForTest_mergeQ_score_rps1_lower_score_than_rps2(void)
{
  std::list<RangePair> rpList;
  SDGString line1 = " \t105085\t110650\t \t?\t?\t?\t?\t?";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(-1);
  rangePairSet1.getRangeS().setNameSeq("-1");
  rangePairSet1.setLength(77259);
  rangePairSet1.setScore(77259);

  // match part   
  SDGString line11 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(1);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(269);
  rpList.push_back(rangePair11);

  SDGString line12 = " \t109951\t110650\t \t119\t700\t9.4e-19\t700\t77.768";
  RangePair rangePair12 = RangePair(line12);
  rangePair12.getRangeQ().setNumChr(1);
  rangePair12.getRangeQ().setNameSeq("");
  rangePair12.getRangeS().setNumChr(2);
  rangePair12.getRangeS().setNameSeq("");
  rangePair12.setLength(700);
  rpList.push_back(rangePair12);
  
  rangePairSet1.setPathDirectly(rpList);
  return rangePairSet1;

}

RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_overlapQ_rps2_include_in_rps1_but_no_overlap(void)
{
  std::list<RangePair> rpList2;
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
  rpList2.push_back(rangePair21);
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

  return rangePairSet2;
}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_overlapQ_rps2_include_in_rps1_but_no_overlap(void)
{
  std::list<RangePair> rpList2;
  SDGString line2 = " \t105400\t105900\t \t3202\t2033\t0\t695\t82.8467";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(501);
  // match part   
  SDGString line21 = " \t105400\t105900\t \t3202\t2033\t0\t695\t82.8467";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeS().setNameSeq("");
  rangePair21.setLength(501);
  rpList2.push_back(rangePair21);

  rangePairSet2.setPathDirectly(rpList2);
  return rangePairSet2;
}


RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_mergeQ_rps2_overlap_on_rps1_frgt(void)
{
	SDGString str1 = "1\t1\t450\t1\t3351\t3238\t0\t441\t76.1062\t441\n";
	SDGString str11 = "\t1\t1\t100\t1\t3351\t3238\t0\t100\t76.1062\t100\n";
	SDGString str12 = "\t1\t110\t450\t1\t3351\t3238\t0\t441\t76.1062\t441\n";


	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);		
	lPath.push_back(str12);		

	return  createRangePairSet(lPath);	

}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_mergeQ_rps2_overlap_on_rps1_frgt(void)
{
	SDGString str2 = "1\t90\t500\t2\t119\t689\t9.4e-19\t411\t77.768\t411\n";
	SDGString str21 = "\t1\t90\t500\t2\t119\t689\t9.4e-19\t411\t77.768\t411\n";
	

	std::list<SDGString> lPath;
	lPath.push_back(str2);		
	lPath.push_back(str21);		

	return  createRangePairSet(lPath);	

}

RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_diffQ_rps2_overlap_on_rps1_frgt_case1(void)
{
	SDGString str1 = "1\t1\t450\t1\t3351\t3238\t0\t441\t76.1062\t441\n";
	SDGString str11 = "\t1\t1\t100\t1\t3351\t3238\t0\t100\t76.1062\t100\n";
	SDGString str12 = "\t1\t200\t450\t1\t3351\t3238\t0\t441\t76.1062\t441\n";


	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);		
	lPath.push_back(str12);		

	return  createRangePairSet(lPath);	

}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_diffQ_rps2_overlap_on_rps1_frgt_case1(void)
{
	SDGString str2 = "1\t50\t500\t2\t119\t689\t9.4e-19\t411\t77.768\t411\n";
	SDGString str21 = "\t1\t50\t500\t2\t119\t689\t9.4e-19\t411\t77.768\t411\n";
	

	std::list<SDGString> lPath;
	lPath.push_back(str2);		
	lPath.push_back(str21);		

	return  createRangePairSet(lPath);	

}

RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_diffQ_rps2_overlap_on_rps1_frgt_case2(void)
{
	SDGString str1 = "1\t1\t450\t1\t3351\t3238\t0\t441\t76.1062\t441\n";
	SDGString str11 = "\t1\t1\t100\t1\t3351\t3238\t0\t100\t76.1062\t100\n";
	SDGString str12 = "\t1\t111\t450\t1\t3351\t3238\t0\t441\t76.1062\t441\n";


	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);		
	lPath.push_back(str12);		

	return  createRangePairSet(lPath);	

}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_diffQ_rps2_overlap_on_rps1_frgt_case2(void)
{
	SDGString str2 = "1\t50\t500\t2\t119\t689\t9.4e-19\t411\t77.768\t411\n";
	SDGString str21 = "\t1\t50\t500\t2\t119\t689\t9.4e-19\t411\t77.768\t411\n";
	

	std::list<SDGString> lPath;
	lPath.push_back(str2);		
	lPath.push_back(str21);		

	return  createRangePairSet(lPath);	

}


RangePairSet Test_RangePairSetUtils::createRangePairSet1ForTest_overlap_first_and_second_fragment(void)
{
	SDGString str1 = "1\t100\t850\t1\t3351\t3238\t0\t10\t76.1062\t652\n";
	SDGString str11 = "\t1\t100\t500\t1\t3351\t3238\t0\t5\t76.1062\t401\n";
	SDGString str12 = "\t1\t600\t850\t1\t3351\t3238\t0\t5\t76.1062\t251\n";
	
	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);		
	lPath.push_back(str12);		

	return  createRangePairSet(lPath);	

}

RangePairSet Test_RangePairSetUtils::createRangePairSet2ForTest_overlap_first_and_second_fragment(void)
{

	SDGString str2 = "1\t1\t900\t2\t119\t689\t9.4e-19\t900\t77.768\t900\n";
	SDGString str21 = "\t1\t1\t900\t2\t119\t689\t9.4e-19\t900\t77.768\t900\n";
	
	std::list<SDGString> lPath;
	
	lPath.push_back(str2);		
	lPath.push_back(str21);		

	return createRangePairSet(lPath);

}

RangePairSet Test_RangePairSetUtils::generateInputs_for_test_orientSubjects_on_minus_strand()
{
	SDGString str1 ="1\t1\t500\t-1\t0\t0\t0\t401\t0\t0\t401\n";
	SDGString str11="\t1\t1\t100\t1\t3351\t3238\t0\t100\t76.1062\t100\n";
	SDGString str12="\t1\t200\t450\t1\t3351\t3238\t0\t251\t76.1062\t251\n";
	SDGString str13="\t1\t451\t500\t2\t576\t689\t9.4e-19\t50\t77.768\t50\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);		
	lPath.push_back(str12);		
	lPath.push_back(str13);		

	return createRangePairSet(lPath);
}


RangePairSet Test_RangePairSetUtils::generateOutputs_for_test_orientSubjects_on_minus_strand()
{
	SDGString str1 ="1\t1\t500\t-1\t0\t0\t0\t401\t0\t0\t401\n";
	SDGString str11="\t1\t1\t100\t1\t3351\t3238\t0\t100\t76.1062\t100\n";
	SDGString str12="\t1\t200\t450\t1\t3351\t3238\t0\t251\t76.1062\t251\n";
	SDGString str13="\t1\t451\t500\t2\t689\t576\t9.4e-19\t50\t77.768\t50\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);		
	lPath.push_back(str12);		
	lPath.push_back(str13);		

	return createRangePairSet(lPath);
}

RangePairSet Test_RangePairSetUtils::generateInputs_for_test_orientSubjects_on_plus_strand()
{
 	SDGString str1 ="1\t1\t500\t-1\t0\t0\t0\t401\t0\t0\t401\n";
	SDGString str11="\t1\t1\t100\t1\t3238\t3351\t0\t100\t76.1062\t100\n";
	SDGString str12="\t1\t200\t450\t1\t3238\t3351\t0\t251\t76.1062\t251\n";
	SDGString str13="\t1\t451\t500\t2\t689\t576\t9.4e-19\t50\t77.768\t50\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);		
	lPath.push_back(str12);		
	lPath.push_back(str13);		


	return createRangePairSet(lPath);
}

RangePairSet Test_RangePairSetUtils::generateOutputs_for_test_orientSubjects_on_plus_strand()
{
 	SDGString str1 ="1\t1\t500\t-1\t0\t0\t0\t401\t0\t0\t401\n";
	SDGString str11="\t1\t1\t100\t1\t3238\t3351\t0\t100\t76.1062\t100\n";
	SDGString str12="\t1\t200\t450\t1\t3238\t3351\t0\t251\t76.1062\t251\n";
	SDGString str13="\t1\t451\t500\t2\t576\t689\t9.4e-19\t50\t77.768\t50\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);		
	lPath.push_back(str12);		
	lPath.push_back(str13);		

	return createRangePairSet(lPath);
}

RangePairSet Test_RangePairSetUtils::generateInputs_for_test_orientSubjects_path_size_is_1()
{
 	SDGString str1 ="1\t1\t500\t-1\t0\t0\t0\t401\t0\t0\t401\n";
	SDGString str11="\t1\t1\t100\t1\t3238\t3351\t0\t100\t76.1062\t100\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);	
	return createRangePairSet(lPath);
}

RangePairSet Test_RangePairSetUtils::generateOutputs_for_test_orientSubjects_path_size_is_1()
{
	SDGString str1 ="1\t1\t500\t-1\t0\t0\t0\t401\t0\t0\t401\n";
	SDGString str11="\t1\t1\t100\t1\t3238\t3351\t0\t100\t76.1062\t100\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);	

	return createRangePairSet(lPath);
}

RangePairSet Test_RangePairSetUtils::generateInputs_for_test_orientSubjects_orient_more_than_1_fragment()
{
	SDGString str1 = "1\t1\t900\t-1\t0\t0\t0\t900\t0\t0\t900\n";
	SDGString str11 = "\t1\t1\t99\t2\t119\t181\t9.4e-19\t99\t77.768\t99\n";
	SDGString str12 = "\t1\t500\t100\t1\t3351\t3238\t0\t401\t76.1062\t401\n";
	SDGString str13 = "\t1\t501\t599\t2\t437\t499\t9.4e-19\t99\t77.768\t99\n";
	SDGString str14 = "\t1\t850\t600\t1\t3351\t3238\t0\t251\t76.1062\t251\n";
	SDGString str15 = "\t1\t851\t900\t2\t659\t689\t9.4e-19\t50\t77.768\t50\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);	
	lPath.push_back(str12);	
	lPath.push_back(str13);	
	lPath.push_back(str14);	
	lPath.push_back(str15);	
	return createRangePairSet(lPath);
}

RangePairSet Test_RangePairSetUtils::generateOutputs_for_test_orientSubjects_orient_more_than_1_fragment()
{
	SDGString str1 = "1\t1\t900\t-1\t0\t0\t0\t900\t0\t0\t900\n";
	SDGString str11 = "\t1\t1\t99\t2\t181\t119\t9.4e-19\t99\t77.768\t99\n";
	SDGString str12 = "\t1\t500\t100\t1\t3351\t3238\t0\t401\t76.1062\t401\n";
	SDGString str13 = "\t1\t501\t599\t2\t499\t437\t9.4e-19\t99\t77.768\t99\n";
	SDGString str14 = "\t1\t850\t600\t1\t3351\t3238\t0\t251\t76.1062\t251\n";
	SDGString str15 = "\t1\t851\t900\t2\t689\t659\t9.4e-19\t50\t77.768\t50\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);	
	lPath.push_back(str12);	
	lPath.push_back(str13);	
	lPath.push_back(str14);	
	lPath.push_back(str15);	
	return createRangePairSet(lPath);
}

RangePairSet Test_RangePairSetUtils::generateInputs_for_test_orientSubjects_no_fragments_to_orient()
{
	SDGString str1 = "1\t1\t900\t-1\t0\t0\t0\t900\t0\t0\t900\n";
	SDGString str11 = "\t1\t1\t99\t2\t181\t119\t9.4e-19\t99\t77.768\t99\n";
	SDGString str12 = "\t1\t500\t100\t1\t3351\t3238\t0\t5\t76.1062\t401\n";
	SDGString str13 = "\t1\t501\t599\t2\t499\t437\t9.4e-19\t99\t77.768\t222\n";
	SDGString str14 = "\t1\t850\t600\t1\t3351\t3238\t0\t5\t76.1062\t251\n";
	SDGString str15 = "\t1\t851\t900\t2\t689\t659\t9.4e-19\t50\t77.768\t50\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);	
	lPath.push_back(str12);	
	lPath.push_back(str13);	
	lPath.push_back(str14);	
	lPath.push_back(str15);	

	return createRangePairSet(lPath);
}

RangePairSet Test_RangePairSetUtils::generateOutputs_for_test_orientSubjects_no_fragments_to_orient()
{
	SDGString str1 = "1\t1\t900\t-1\t0\t0\t0\t900\t0\t0\t900\n";
	SDGString str11 = "\t1\t1\t99\t2\t181\t119\t9.4e-19\t99\t77.768\t99\n";
	SDGString str12 = "\t1\t500\t100\t1\t3351\t3238\t0\t5\t76.1062\t401\n";
	SDGString str13 = "\t1\t501\t599\t2\t499\t437\t9.4e-19\t99\t77.768\t222\n";
	SDGString str14 = "\t1\t850\t600\t1\t3351\t3238\t0\t5\t76.1062\t251\n";
	SDGString str15 = "\t1\t851\t900\t2\t689\t659\t9.4e-19\t50\t77.768\t50\n";

	std::list<SDGString> lPath;
	lPath.push_back(str1);		
	lPath.push_back(str11);	
	lPath.push_back(str12);	
	lPath.push_back(str13);	
	lPath.push_back(str14);	
	lPath.push_back(str15);	

	return createRangePairSet(lPath);
}

RangePairSet Test_RangePairSetUtils::createRangePairSet(std::list<SDGString> inputList)
{

	RangePairSet currentRps;
	for(std::list<SDGString>::iterator input_it=inputList.begin(); input_it!=inputList.end();input_it++){
	
		size_t posFirstTab = input_it->find_first_of("\t");
		size_t posLastTab = input_it->find_last_of("\t");
		size_t posEndLine = input_it->find_first_of("\\n");
		size_t lengthSeq = posEndLine - posLastTab;
		size_t lengthInit = posLastTab - posFirstTab;
		
		// isARangePairSetToCreate		
		if (posFirstTab != 0){
			// createARangePairSet 
			currentRps = RangePairSet(input_it->substr(0, lengthInit));
			SDGString sNameSeq = currentRps.getRangeS().getNameSeq();
			SDGString qNameSeq = currentRps.getRangeQ().getNameSeq();
			currentRps.getRangeS().setNumChr(atol(sNameSeq));
			currentRps.getRangeQ().setNumChr(atol(qNameSeq));
			currentRps.setLength(atol(input_it->substr(posLastTab, lengthSeq)));
			// empty path
			std::list<RangePair> path;
			currentRps.setPathDirectly(path);

		// isARangePairToCreate	
		}else{	
			// createARangePair
			RangePair rp = RangePair(input_it->substr(posFirstTab, lengthInit));
			SDGString sNameSeq = rp.getRangeS().getNameSeq();
			SDGString qNameSeq = rp.getRangeQ().getNameSeq();
			rp.getRangeS().setNumChr(atol(sNameSeq));
			rp.getRangeQ().setNumChr(atol(qNameSeq));
			rp.setLength(atol(input_it->substr(posLastTab, lengthSeq)));
			
			// updatePath
			std::list<RangePair> path = currentRps.getPath();
			std::list<RangePair> newPath;
			for(std::list<RangePair>::iterator lrp_it=path.begin(); lrp_it!=path.end();lrp_it++){
				RangePair newRp = *lrp_it;	
            			newPath.push_back(newRp); 
    			}
			newPath.push_back(rp);
			currentRps.setPathDirectly(newPath);

		}

	}
	return currentRps;
}
