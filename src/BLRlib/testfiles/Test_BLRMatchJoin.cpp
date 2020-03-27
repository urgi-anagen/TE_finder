#include "Test_BLRMatchJoin.h"
#include "SDGString.h"
#include "FileUtils.h"
#include "BLRJoinParameter.h"
#include "BLRMatchJoin.h"
#include "Range.h"
#include "RangePair.h"
#include "RangePairSet.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_BLRMatchJoin );
//---------------------------------------------------------------------------
void Test_BLRMatchJoin::setUp()
{
}
//---------------------------------------------------------------------------
void Test_BLRMatchJoin::tearDown()
{
}
//---------------------------------------------------------------------------
void Test_BLRMatchJoin::test_join()
{
	std::ostringstream inputData;
	inputData<<"chunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100.00\n";
	inputData<<"chunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100.00\n";
	inputData<<"chunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100.00\n";
	inputData<<"chunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100.00\n";
	inputData<<"chunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80.00\n";
	inputData<<"chunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90.00\n";
    inputData<<"chunk321\t18184\t18354\trefTE_365\t655\t472\t0\t110\t78.8\n";
    inputData<<"chunk321\t18704\t18882\trefTE_365\t180\t2\t0\t137\t78\n";
    inputData<<"chunk321\t18188\t18882\trefTE_369\t709\t3\t0\t1105\t98.4\n";


	BLRJoinParameter para = Test_BLRUtils::createParameter();
	BLRMatchJoin matchJoin(para);

	std::istringstream inputDataStream(inputData.str());
	BLRMatchAlign match_align;
	match_align.read(para,inputDataStream,0);

	BLRMatchPath match_path;

	matchJoin.join(match_align,match_path,0);

	std::ostringstream obs;
	match_path.write(obs);
	match_path.writeAttribute(obs);

	std::ostringstream exp;
	exp<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
	exp<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
	exp<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
	exp<<"2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100\n";
	exp<<"3\tchunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80\n";
	exp<<"3\tchunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90\n";
    exp<<"4\tchunk321\t18184\t18354\trefTE_365\t655\t472\t0\t110\t78.8\n";
    exp<<"4\tchunk321\t18704\t18882\trefTE_365\t180\t2\t0\t137\t78\n";
    exp<<"5\tchunk321\t18188\t18882\trefTE_369\t709\t3\t0\t1105\t98.4\n";
	exp<<"[1\tchunk1\t100\t2000\trefTE_1\t1\t3000\t0\t1725\t100]\n";
	exp<<"[2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100]\n";
	exp<<"[3\tchunk1\t10000\t11500\trefTE_2\t1\t1500\t0\t1380\t82.8602]\n";
    exp<<"[4\tchunk321\t18184\t18882\trefTE_365\t655\t2\t0\t215\t78.4055]\n";
    exp<<"[5\tchunk321\t18188\t18882\trefTE_369\t709\t3\t0\t1105\t98.4]\n";


	CPPUNIT_ASSERT_EQUAL(exp.str(), obs.str());
}
//---------------------------------------------------------------------------
void Test_BLRMatchJoin::test_nojoin()
{
	std::ostringstream inputData;
	inputData<<"chunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100.00\n";
	inputData<<"chunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100.00\n";
	inputData<<"chunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100.00\n";
	inputData<<"chunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100.00\n";
	inputData<<"chunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80.00\n";
	inputData<<"chunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90.00\n";

	BLRJoinParameter para = Test_BLRUtils::createParameter();
	BLRMatchJoin matchJoin(para);

	std::istringstream inputDataStream(inputData.str());
	BLRMatchAlign match_align;
	match_align.read(para,inputDataStream,0);

	BLRMatchPath match_path;

	matchJoin.noJoin(match_align,match_path,0);

	std::ostringstream obs;
	match_path.write(obs);
	match_path.writeAttribute(obs);

	std::ostringstream exp;
	exp<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
	exp<<"2\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
	exp<<"3\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
	exp<<"4\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100\n";
	exp<<"5\tchunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80\n";
	exp<<"6\tchunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90\n";
	exp<<"[1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100]\n";
	exp<<"[2\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100]\n";
	exp<<"[3\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100]\n";
	exp<<"[4\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t1000\t100]\n";
	exp<<"[5\tchunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80]\n";
	exp<<"[6\tchunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90]\n";


	CPPUNIT_ASSERT_EQUAL(exp.str(), obs.str());
}
//---------------------------------------------------------------------------
void Test_BLRMatchJoin::test_clean_conflicts()
{
	std::ostringstream inputData;
	inputData<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
	inputData<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
	inputData<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
	inputData<<"2\tchunk1\t400\t800\trefTE_2\t1\t400\t0\t400\t80\n";
    inputData<<"3\tchunk1\t100\t500\trefTE_3\t1\t400\t0\t400\t90\n";
    inputData<<"3\tchunk1\t600\t1000\trefTE_3\t500\t900\t0\t400\t90\n";
    inputData<<"3\tchunk1\t1001\t2000\trefTE_3\t2000\t3000\t0\t1000\t90\n";

	BLRJoinParameter para = Test_BLRUtils::createParameter();
	BLRMatchJoin matchJoin(para);

	std::istringstream inputDataStream(inputData.str());
	BLRMatchPath match_path;
	match_path.read(para,inputDataStream,0);

	matchJoin.clean_conflicts(match_path,0);

	std::ostringstream obs;
	match_path.write(obs);
	match_path.writeAttribute(obs);

	std::ostringstream exp;
	exp<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100\n";
	exp<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100\n";
	exp<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
	exp<<"2\tchunk1\t501\t599\trefTE_2\t102\t200\t0\t99\t80\n";
	exp<<"[1\tchunk1\t100\t2000\trefTE_1\t1\t3000\t0\t1725\t100]\n";
	exp<<"[2\tchunk1\t501\t599\trefTE_2\t102\t200\t0\t99\t80]\n";

	CPPUNIT_ASSERT_EQUAL(exp.str(), obs.str());
}
//---------------------------------------------------------------------------
void Test_BLRMatchJoin::test_merge(void){
	std::ostringstream inputData;
	inputData<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t400\t100.00\n";
	inputData<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t400\t100.00\n";
	inputData<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100.00\n";
	inputData<<"2\tchunk1\t1500\t2500\trefTE_2\t1500\t2500\t0\t999\t100.00\n";
	inputData<<"3\tchunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80.00\n";
	inputData<<"3\tchunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90.00\n";

	BLRJoinParameter para = Test_BLRUtils::createParameter();
	BLRMatchJoin matchJoin(para);

	std::istringstream inputDataStream(inputData.str());
	BLRMatchPath match_path;
	match_path.read(para,inputDataStream,0);

	matchJoin.merge(match_path,0);

	std::ostringstream obs;
	match_path.write(obs);
	match_path.writeAttribute(obs);

	std::ostringstream exp;
	exp<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t401\t100\n";
	exp<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t401\t100\n";
	exp<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
	exp<<"1\tchunk1\t2001\t2500\trefTE_2\t2001\t2500\t0\t500\t100\n";
	exp<<"2\tchunk1\t10000\t11000\trefTE_2\t1\t1000\t0\t1000\t80\n";
	exp<<"2\tchunk1\t11100\t11500\trefTE_2\t1100\t1500\t0\t400\t90\n";
	exp<<"[1\tchunk1\t100\t2500\t-1\t0\t0\t0\t2302\t100]\n";
	exp<<"[2\tchunk1\t10000\t11500\trefTE_2\t1\t1500\t0\t1380\t82.8602]\n";

	CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}
//---------------------------------------------------------------------------
void Test_BLRMatchJoin::test_split(void){
	std::ostringstream inputData;
	inputData<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t401\t100\n";
	inputData<<"1\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t401\t100\n";
	inputData<<"1\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
	inputData<<"2\tchunk1\t501\t599\trefTE_2\t102\t200\t0\t99\t80\n";

	BLRJoinParameter para = Test_BLRUtils::createParameter();
	BLRMatchJoin matchJoin(para);

	std::istringstream inputDataStream(inputData.str());
	BLRMatchPath match_path;
	match_path.read(para,inputDataStream,0);

	matchJoin.split(match_path);

	std::ostringstream obs;
	match_path.write(obs);
	match_path.writeAttribute(obs);

	std::ostringstream exp;
	exp<<"1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t401\t100\n";
    exp<<"2\tchunk1\t501\t599\trefTE_2\t102\t200\t0\t79\t80\n";
	exp<<"3\tchunk1\t600\t1000\trefTE_1\t500\t900\t0\t401\t100\n";
	exp<<"3\tchunk1\t1001\t2000\trefTE_1\t2000\t3000\t0\t1000\t100\n";
	exp<<"[1\tchunk1\t100\t500\trefTE_1\t1\t400\t0\t401\t100]\n";
    exp<<"[2\tchunk1\t501\t599\trefTE_2\t102\t200\t0\t79\t80]\n";
	exp<<"[3\tchunk1\t600\t2000\trefTE_1\t500\t3000\t0\t1401\t100]\n";



	CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}























