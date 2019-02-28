/*
 * \file Test_FileUtils.cpp
 */

#include <vector>
#include <ctime>
#include <cstdio>
using namespace std;

#include "Test_FileUtils.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_FileUtils );

void Test_FileUtils::setUp()
{
}

void Test_FileUtils::tearDown()
{
	line.str("");
}

void Test_FileUtils::test_doesFileExist( void )
{
	time_t rawTime;
	time( &rawTime );
	struct tm * timeInfo;
	timeInfo = localtime( &rawTime );
	char buffer [80];
	strftime( buffer, 80, "%Y%m%d%H%M%S", timeInfo );

	string file = "skgb7362_";
	file.push_back( *buffer );
	bool exp = false;
	bool obs = FileUtils::doesFileExist( file );
	CPPUNIT_ASSERT_EQUAL( exp, obs );
	remove( file.c_str() );
}

void Test_FileUtils::test_areTwoFilesIdentical( void )
{
	string inFile1 = "dummyFile1";
	ofstream inFileStream1;
	inFileStream1.open( inFile1.c_str() );
	line<<"HAV)[_*^HVDCZEY\n";
	line<<"BbjkucegzjnBJKU\n";
	line<<"hbv1625738783E4\n";
	inFileStream1 << line.str();
	inFileStream1.close();

	line.str("");

	string inFile2 = "dummyFile2";
	ofstream inFileStream2;
	inFileStream2.open( inFile2.c_str() );
	line<<"HAV)[_*^HVDCZEY\n";
	line<<"BbjkucegzjnBJKU\n";
	line<<"hbv1625738783E4\n";
	inFileStream2 << line.str();
	inFileStream2.close();

	bool exp = true;
	bool obs = FileUtils::areTwoFilesIdentical( inFile1, inFile2 );
	CPPUNIT_ASSERT_EQUAL( exp, obs );

	remove( inFile1.c_str() );
	remove( inFile2.c_str() );
}

void Test_FileUtils::test_removeLinesFromFile( void )
{
	string inFile = "dummyInFile";
	ofstream inFileStream;
	inFileStream.open( inFile.c_str() );
	line<<"HAV)[_*^HVDCZEY\n";
	line<<"BHDXV19736HZE\n";
	line<<"BbjkucegzjnBJKU\n";
	line<<"hbv1625738783E4\n";
	line<<"hbv1625738783E4KKKKKKKKKKKK\n";
	inFileStream << line.str();
	inFileStream.close();

	line.str("");

	string expFile = "dummyExpFile";
	ofstream expFileStream;
	expFileStream.open( expFile.c_str() );
	line<<"HAV)[_*^HVDCZEY\n";
	line<<"BbjkucegzjnBJKU\n";
	line<<"hbv1625738783E4KKKKKKKKKKKK\n";
	expFileStream << line.str();
	expFileStream.close();

	string obsFile = "dummyObsFile";
	vector< unsigned long> vLinesToRmv;
	vLinesToRmv.push_back( 1 );
	vLinesToRmv.push_back( 3 );
	FileUtils::removeLinesFromFile( inFile, obsFile, vLinesToRmv );

	bool exp = true;
	bool obs = FileUtils::areTwoFilesIdentical( expFile, obsFile );
	CPPUNIT_ASSERT_EQUAL( exp, obs );

	remove( inFile.c_str() );
	remove( expFile.c_str() );
	remove( obsFile.c_str() );
}


void Test_FileUtils::test_removeFile( void )
{
	string inFile = "dummyInFile";
	ofstream inFileStream;
	inFileStream.open( inFile.c_str() );
	line<<"HAV)[_*^HVDCZEY\n";
	inFileStream << line.str();
	inFileStream.close();

	FileUtils::removeFile( inFile );

	bool exp = false;
	bool obs = FileUtils::doesFileExist( inFile );
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}
