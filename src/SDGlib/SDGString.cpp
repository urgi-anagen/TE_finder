/*
 * \file SDGString.cpp
 */

#include "SDGString.h"

std::vector<SDGString> SDGString::tokenize( const SDGString& str, const SDGString& delimiters )
{
	std::vector<SDGString> tokens;
	// skip delimiters at beginning
	SDGString::size_type lastPos = str.find_first_not_of( delimiters, 0 );
	// find first "non-delimiter"
	SDGString::size_type pos = str.find_first_of( delimiters, lastPos );
	while( SDGString::npos != pos || SDGString::npos != lastPos )
		{
			// found a token, add it to the vector.
			tokens.push_back( str.substr( lastPos, pos - lastPos ) );
			// skip delimiters. Note the "not_of".
			lastPos = str.find_first_not_of( delimiters, pos );
			// find next "non-delimiter"
			pos = str.find_first_of( delimiters, lastPos );
		}
	return( tokens );
}
