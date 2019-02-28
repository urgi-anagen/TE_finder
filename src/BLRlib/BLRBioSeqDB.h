/**
 * \file BLRBioSeqDB.h
 * \brief Handle input fasta files given to and cut by BLRBlaster.
 */

#ifndef BLRBIOSEQDB_H
#define BLRBIOSEQDB_H

#include <iostream>
#include <stdlib.h>
#include <SDGString.h>
#include <map>
#include <RangeSeq.h>

class BLRBioSeqDB
{
	std::map<unsigned long,RangeSeq> num2infoseq;

	public:

		//! A constructor.
		BLRBioSeqDB( void ) {};

		//! A constructor.
		BLRBioSeqDB( SDGString dbfilename ) { init( dbfilename ); };

		/** \fn void init( SDGString dbfilename );
		 * \brief Load the fasta file that was cut
		 * \param dbfilename name of the fasta file
		 */
		void init( SDGString dbfilename );

		/** \fn RangeSeq& getRange( unsigned long num )
		 * \brief Return the Range object associated to the given identifier
		 * \param num identifier of a given sequence
		 * \return Range instance
		 */
		RangeSeq& getRange( unsigned long num )
		{
			if( num2infoseq.find(num) == num2infoseq.end() )
			{
				std::cout<<num<<" not found!"<<std::endl;
				exit( EXIT_FAILURE );
			}
			return num2infoseq[ num ];
		};

		bool operator==(const BLRBioSeqDB& db) const
		{
			return num2infoseq==db.num2infoseq;
		};
		bool operator!=(const BLRBioSeqDB& db) const
		{
			return num2infoseq==db.num2infoseq;
		};
};

#endif
