/**
 *
 * BLRMatchList.h
 *
 */

#ifndef BLRMATCHLIST_H
#define BLRMATCHLIST_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <iterator>
#include <iomanip>
#include <BlastMatch.h>
#include "../blaster/BLRBlasterParameter.h"

class BLRMatchList 
{
 private:
  BlastMatch match;
 
 public: 
  std::list<BlastMatch> match_list;
  //! default constructor
  BLRMatchList(): match_list(){};
  //! destructor
  virtual ~BLRMatchList();
  //! - set empty match_list 
  void clear(); 
  //! - Function to remove the matches that do not need
  void remove(BLRBlasterParameter para);
  //! - Function to fill matchlist with matches found 
  void fill_ncbilist(SDGString& filename);
  void fill_wulist(SDGString& filename);
  //!- Function for developement to show the alignlist
  void readncbiblastfield(std::istream& in);
  void readwublastfield(std::istream& in);
  void show_list();
  //!- Function to save the matchlist
  void save_list(SDGString, int flag);
  //!- function to get size of the membre of matchlist
  void savebin(std::ostream& out);
  void writetxt(std::ostream& out);
  void readtxt(std::istream& in);
  void readtxt(char* data);

  unsigned getSize()
    {
      return match_list.size();
    };

  void* clone() const
    {
      return (void*) new BLRMatchList(*this);
    };
  
}; // BLRMatchList

#endif




