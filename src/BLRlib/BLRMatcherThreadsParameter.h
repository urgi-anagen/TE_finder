/**
 * \file BLRMatcherParameter.h
 * \brief Header file for the class BLRMatcherParameter
 */

#ifndef BLRMATCHERPARAMETER_H
#define BLRMATCHERPARAMETER_H

#include <cstdlib>
#include <iostream>
#include <time.h>
#include <fstream>
#include <SDGString.h>
#include <list>
#include "BLRJoinParameter.h"


/**
 * \class BLRMatcherThreadsParameter
 * \brief Parameters for the MATCHER program
 */
class BLRMatcherThreadsParameter: public BLRJoinParameter
{
 private:
  bool clean_before, clean_after, merge;
  int verbose;
  unsigned nb_sets;

 public:
  //! -Constructor
  BLRMatcherThreadsParameter(void){reset();};
  
  //! -Destructor

  virtual ~BLRMatcherThreadsParameter(void){};

  void reset(void)
    {
      BLRJoinParameter::reset();
      clean_after=false;
      clean_before=true;
      merge=false;
      nb_sets=1;
      verbose=0;
    };
  void parseOptArg (int numarg, char *tabarg[]);
  void help (void);
  void view(std::ostream& out) const;
  void write(const SDGString& filename) const;

  void setMerge(bool set_merge){ merge=set_merge;};

  bool getCleanBefore(void) {return  clean_before;};
  void setCleanBefore(bool clean) {clean_before = clean;};

  bool getCleanAfter(void) {return  clean_after;};
  void setCleanAfter(bool clean) {clean_after = clean;};

  unsigned getNbSets(void) {return nb_sets;};
  bool getMerge(void) {return merge;};

  int getVerbose(void) const {return verbose;};

};

#endif
