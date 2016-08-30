/**
 * \file BLRBioSeqDB.cpp
 */

#include <fstream>
#include "BLRBioSeqDB.h"

void BLRBioSeqDB::init( SDGString dbfilename )
{
	std::ifstream file( dbfilename.start() );
	char buff[2048];
	while( file )
	{
		file.getline( buff, 2048 );
		if( *(buff) == '>' )
		{
			SDGString chr_name(buff), str_num;
			str_num = chr_name.beforematch( "[ \t]+" );
			str_num = str_num.substr(1);

			SDGString range_str = chr_name.aftermatch( "\\{Cut\\}[ \t]+" );
			SDGString start_str = range_str.beforematch( "\\.\\." );
			unsigned long start = atol( start_str.start() );

			SDGString end_str = range_str.aftermatch( "\\.\\." );
			unsigned long end = atol( end_str.start() );

			chr_name = chr_name.aftermatch( "[0-9]+[ \t]+" );
			chr_name = chr_name.beforematch( "[ \t]+\\{Cut\\}" );

			RangeSeq range( str_num, chr_name, start, end );

			num2infoseq[ atol(str_num.start()) ] = range;
		}
	}
	file.close();
};
