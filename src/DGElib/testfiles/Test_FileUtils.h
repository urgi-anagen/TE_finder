/**
 * \file Test_FileUtils.h
 * \brief Unitary tests for the class FileUtils
 */

#ifndef TEST_FILEUTILS_H
#define TEST_FILEUTILS_H

#include <sstream>
using namespace std;

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "FileUtils.h"

/**
 * \class Test_FileUtils
 * \brief Unitary tests for the class FileUtils
 */
class Test_FileUtils: public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( Test_FileUtils );
	CPPUNIT_TEST( test_doesFileExist );
	CPPUNIT_TEST( test_areTwoFilesIdentical );
	CPPUNIT_TEST( test_removeLinesFromFile );
	CPPUNIT_TEST_SUITE_END();

public:
	stringstream line;
	void setUp();
	void tearDown();

protected:
	void test_doesFileExist( void );
	void test_areTwoFilesIdentical( void );
	void test_removeLinesFromFile( void );
	void test_removeFile( void );
};

#endif
