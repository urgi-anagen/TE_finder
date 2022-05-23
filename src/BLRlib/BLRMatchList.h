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
#include <BLRJoinParameter.h>

class BLRMatchList 
{
 private:
  BlastMatch match;
 
 public: 
  std::list<BlastMatch> match_list;
  //! default constructor
  BLRMatchList(void): match_list(){};
  //! destructor
  virtual ~BLRMatchList();
  //! - set empty match_list 
  void clear(); 
  //! - Function to remove_self_hits the matches that do not need
  void read_blast_results(const std::string& match_filename, double max_eval, unsigned min_length, double min_identity, bool is_wublast);
  void remove_self_hits(unsigned cut_length, unsigned cut_over);
  void show_list();
  //!- Function to save the matchlist
  void save_list(SDGString, unsigned flag);
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




