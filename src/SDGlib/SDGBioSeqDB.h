/**
 * \file SDGBioSeqDB.h
 * \brief Header file for the class SDGBioSeqDB
 */

#ifndef SDGBIOSEQDB
#define SDGBIOSEQDB

#include <vector>
#include <map>
#include <SDGReference.h>
#include <SDGString.h>
#include <SDGMemBioSeq.h>

/**
 * \class SDGBioSeqDB
 * \brief Handle a fasta file as a vector of SDGBioSeq instances
 */
class SDGBioSeqDB: public std::vector<SDGBioSeq>
{
  std::map<SDGString,unsigned> name2pos;
  SDGString filename;

 public:

  SDGBioSeqDB()= default;
  explicit SDGBioSeqDB(SDGString fichier, int verbose=0){load(fichier,verbose);};
  void load( SDGString fichier, int verbose=0 );

  SDGBioSeq& find(SDGString name)
    { return operator[](name2pos[name]); };

//  void load_idx(void){};
//  void save_idx(void){};
  static SDGString name() { return "SDGBioSeqDB"; };
  unsigned long getSize() { return name2pos.size(); };

  bool operator==(const SDGBioSeqDB& bd)
	{
	  return filename==bd.filename;
	};
  bool operator!=(const SDGBioSeqDB& bd)
 	{
 	  return filename!=bd.filename;
 	};
};

#endif
