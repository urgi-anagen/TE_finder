/*
 * \file FileUtils.cpp
 */

#include <string>
#include <fstream>
#include <sys/stat.h>
#include <algorithm>
#include <iostream>

#include "FileUtils.h"

bool FileUtils::doesFileExist( SDGString file )
{
	struct stat stFileInfo;
	int intStat;
	intStat = stat( file.c_str(), &stFileInfo );
	if( intStat == 0 )
		return( true );
	else
		return( false );
}

bool FileUtils::areTwoFilesIdentical( SDGString file1, SDGString file2 )
{
	bool result = true;
	std::ifstream stream1;
	stream1.open( file1.c_str() );
	std::ifstream stream2;
	stream2.open( file2.c_str() );
	char ch1, ch2;
	while( stream1.get( ch1 ) )
		if( ! stream2.get( ch2 ) || ( ch1 != ch2 ) )
			result = false;
	if( stream2.get( ch2 ) )
		result = false;
	stream1.close();
	stream2.close();
	return( result );
}

void FileUtils::removeLinesFromFile( SDGString inFile, SDGString outFile, std::vector<unsigned long> vLinesToRmv )
{
	std::ifstream inFileStream;
	inFileStream.open( inFile.c_str() );
	std::ofstream outFileStream;
	outFileStream.open( outFile.c_str() );
	SDGString line;
	unsigned long countLines = 0;
	std::vector<unsigned long>::iterator it;
	while( ! inFileStream.eof() )
	{
		getline( inFileStream, line );
		if( line == "" )
			break;
		countLines += 1;
		it = find( vLinesToRmv.begin(), vLinesToRmv.end(), countLines-1 );
		if( it == vLinesToRmv.end() )
			outFileStream << line << "\n";
	}
	inFileStream.close();
	outFileStream.close();
}

void FileUtils::removeFile( SDGString file )
{
	if( doesFileExist( file ) )
		remove( file.c_str() );
	else
		std::cout<<"WARNING: can't remove file '"<<file<<"'"<<std::endl;
}

void FileUtils::openFile( SDGString file, std::ofstream &fileStream )
{
	fileStream.open( file.c_str(), std::ofstream::out | std::ofstream::trunc );
	if ( ! fileStream.is_open() )
	{
		std::cerr<<"ERROR: unable to open file '"<<file<<"'"<<std::endl;
		std::exit( EXIT_FAILURE );
	}
}
