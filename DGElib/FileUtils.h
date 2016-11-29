/**
 * \file FileUtils.h
 * \brief Header file for the class FileUtils
 */

#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <string>
#include <vector>
#include <iostream>
#include "SDGString.h"

/**
 * \class FileUtils
 * \brief Static methods useful when handling files
 */
class FileUtils
{
public:
	static bool doesFileExist( SDGString file );
	static bool areTwoFilesIdentical( SDGString file1, SDGString file2 );
	static void removeLinesFromFile( SDGString inFile, SDGString outFile, std::vector<unsigned long> vLinesToRmv );
	static void removeFile( SDGString file );
	static void openFile( SDGString file, std::ofstream &fileStream );
};

#endif
