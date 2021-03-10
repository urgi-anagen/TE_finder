#include "Test_BLRMatchPath.h"
#include "SDGString.h"
#include "FileUtils.h"
#include "../../matcher/BLRMatcherParameter.h"
#include "BLRMatchPath.h"
#include "Range.h"
#include "RangePair.h"
#include "RangePairSet.h"
#include "Test_BLRUtils.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_BLRMatchPath );
//---------------------------------------------------------------------------
void Test_BLRMatchPath::setUp()
{
}
//---------------------------------------------------------------------------
void Test_BLRMatchPath::tearDown()
{
}
//---------------------------------------------------------------------------
void Test_BLRMatchPath::test_load()
{
   BLRJoinParameter para = Test_BLRUtils::createParameter();
   BLRMatchPath matchPath;
   SDGString path_file = "match.path";

	std::ostringstream inputData;
    inputData<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
    inputData<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
    inputData<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
    inputData<<"2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100\n";
    inputData<<"3\tchunk2\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80\n";
    inputData<<"3\tchunk2\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90\n";

	std::ofstream obs(path_file);
	obs<<inputData.str();
	obs.close();

	matchPath.load(para,path_file);

	std::ostringstream sout_obs;
	matchPath.write(sout_obs);
	matchPath.writeAttribute(sout_obs);

	std::ostringstream sout_exp;
    sout_exp<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
    sout_exp<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
    sout_exp<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
    sout_exp<<"2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100\n";
    sout_exp<<"3\tchunk2\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80\n";
    sout_exp<<"3\tchunk2\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90\n";
    sout_exp<<"[1\tchunk1\t100\t2000\trefTE_1\t1\t3000\t0\t1725\t100]\n";
    sout_exp<<"[2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100]\n";
    sout_exp<<"[3\tchunk2\t10000\t11500\trefTE_2\t1\t1500\t0\t1380\t82.8602]\n";



	CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
	FileUtils::removeFile(path_file);

	CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
}
//---------------------------------------------------------------------------
void Test_BLRMatchPath::test_read(void)
{
	BLRJoinParameter para = Test_BLRUtils::createParameter();
	BLRMatchPath matchPath;

	std::ostringstream inputData;
    inputData<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
    inputData<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
    inputData<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
    inputData<<"2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100\n";
    inputData<<"3\tchunk2\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80\n";
    inputData<<"3\tchunk2\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90\n";

	std::istringstream inputDataStream(inputData.str());
	matchPath.read(para,inputDataStream, 0);

	std::ostringstream sout_obs;
	matchPath.write(sout_obs);
	matchPath.writeAttribute(sout_obs);

	std::ostringstream sout_exp;
    sout_exp<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
    sout_exp<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
    sout_exp<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
    sout_exp<<"2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100\n";
    sout_exp<<"3\tchunk2\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80\n";
    sout_exp<<"3\tchunk2\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90\n";
    sout_exp<<"[1\tchunk1\t100\t2000\trefTE_1\t1\t3000\t0\t1725\t100]\n";
    sout_exp<<"[2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100]\n";
    sout_exp<<"[3\tchunk2\t10000\t11500\trefTE_2\t1\t1500\t0\t1380\t82.8602]\n";


	CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
}
//---------------------------------------------------------------------------
void Test_BLRMatchPath::test_setFromRpsList(void)
{
    BLRJoinParameter para = Test_BLRUtils::createParameter();
    BLRMatchPath matchPath,matchPath2;

    std::ostringstream inputData;
    inputData<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
    inputData<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
    inputData<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
    inputData<<"2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100\n";
    inputData<<"3\tchunk2\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80\n";
    inputData<<"3\tchunk2\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90\n";

    std::istringstream inputDataStream(inputData.str());
    matchPath.read(para,inputDataStream, 0);

    std::list<RangePairSet> rps_list=matchPath.getRpsListFromMapPath();
    matchPath2.setFromRpsList(para,rps_list,0);

    std::ostringstream sout_obs;
    matchPath2.write(sout_obs);
    matchPath2.writeAttribute(sout_obs);

    std::ostringstream sout_exp;
    sout_exp<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
    sout_exp<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
    sout_exp<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
    sout_exp<<"2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100\n";
    sout_exp<<"3\tchunk2\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80\n";
    sout_exp<<"3\tchunk2\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90\n";
    sout_exp<<"[1\tchunk1\t100\t2000\trefTE_1\t1\t3000\t0\t1725\t100]\n";
    sout_exp<<"[2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100]\n";
    sout_exp<<"[3\tchunk2\t10000\t11500\trefTE_2\t1\t1500\t0\t1380\t82.8602]\n";

    CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
}
//---------------------------------------------------------------------------
void Test_BLRMatchPath::test_setFromRpsList2(void)
{
    BLRJoinParameter para = Test_BLRUtils::createParameter();
    BLRMatchPath matchPath1,matchPath2,matchPath3;

    std::ostringstream inputData1;
    inputData1<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
    inputData1<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
    inputData1<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
    inputData1<<"2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100\n";
    inputData1<<"3\tchunk2\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80\n";
    inputData1<<"3\tchunk2\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90\n";

    std::istringstream inputDataStream1(inputData1.str());
    matchPath1.read(para,inputDataStream1, 0);

    std::ostringstream inputData2;
    inputData2<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
    inputData2<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
    inputData2<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
    inputData2<<"2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100\n";
    inputData2<<"3\tchunk3\t10000\t11000\trefTE_3\t1\t1000\t0\t1000\t80\n";
    inputData2<<"3\tchunk3\t11100\t11500\trefTE_3\t1100\t1500\t0\t400\t90\n";

    std::istringstream inputDataStream2(inputData2.str());
    matchPath2.read(para,inputDataStream2, 0);

    std::list<RangePairSet> rps_list1=matchPath1.getRpsListFromMapPath();
    std::list<RangePairSet> rps_list2=matchPath2.getRpsListFromMapPath();
    matchPath3.setFromRpsList(para,rps_list1,0);
    matchPath3.setFromRpsList(para,rps_list2,0);

    std::ostringstream sout_obs;
    matchPath3.write(sout_obs);
    matchPath3.writeAttribute(sout_obs);

    std::ostringstream sout_exp;
    sout_exp<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
    sout_exp<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
    sout_exp<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
    sout_exp<<"2\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
    sout_exp<<"2\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
    sout_exp<<"2\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
    sout_exp<<"3\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100\n";
    sout_exp<<"4\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100\n";
    sout_exp<<"5\tchunk2\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80\n";
    sout_exp<<"5\tchunk2\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90\n";
    sout_exp<<"6\tchunk3\t10000\t11000\trefTE_3\t1\t1000\t0\t1000\t80\n";
    sout_exp<<"6\tchunk3\t11100\t11500\trefTE_3\t1100\t1500\t0\t400\t90\n";
    sout_exp<<"[1\tchunk1\t100\t2000\trefTE_1\t1\t3000\t0\t1725\t100]\n";
    sout_exp<<"[2\tchunk1\t100\t2000\trefTE_1\t1\t3000\t0\t1725\t100]\n";
    sout_exp<<"[3\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100]\n";
    sout_exp<<"[4\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100]\n";
    sout_exp<<"[5\tchunk2\t10000\t11500\trefTE_2\t1\t1500\t0\t1380\t82.8602]\n";
    sout_exp<<"[6\tchunk3\t10000\t11500\trefTE_3\t1\t1500\t0\t1380\t82.8602]\n";

    CPPUNIT_ASSERT_EQUAL(sout_exp.str(), sout_obs.str());
}
//---------------------------------------------------------------------------
void Test_BLRMatchPath::test_writeBED(void){

	BLRJoinParameter para = Test_BLRUtils::createParameter();
	BLRMatchPath matchPath;

	std::ostringstream inputData;
	inputData << "10\tCHR1v01212004\t100\t250\tTNAT1A\t36\t523\t2e-128\t460\t81.56\n";
	inputData << "10\tCHR1v01212004\t800\t1000\tTNAT1A\t36\t523\t2e-128\t460\t81.56\n";
	inputData << "20\tCHR2v01212004\t1050\t2000\tTNAT1A\t580\t480\t7e-78\t87\t79.14\n";

	std::istringstream inputDataStream(inputData.str());
	matchPath.read(para,inputDataStream, 0);



	std::ostringstream obs;
	SDGString color = "255,127,0";

	matchPath.writeBED(obs);

	std::ostringstream exp;
	exp << "CHR1v01212004\t100\t1000\tTNAT1A,TNAT1A\t868\t+\t100\t1000\t2\t151,201,\t0,700,\n";
	exp << "CHR2v01212004\t1050\t2000\tTNAT1A\t87\t-\t1050\t2000\n";

	CPPUNIT_ASSERT_EQUAL(exp.str(), obs.str());
}




















