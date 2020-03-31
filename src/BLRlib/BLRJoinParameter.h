/**
 * \file BLRJoinParameter.h
 * \brief Header file for the class BLRJoinParameter
 */

#ifndef BLRJOINPARAMETER_H
#define BLRJOINPARAMETER_H

#include <cstdlib>
#include <iostream>
#include <time.h>
#include <fstream>
#include <SDGString.h>
#include <list>

/**
 * \class BLRBJoinParameter
 * \brief Parameters for the "join" procedure used in MATCHER and GROUPER
 */
class BLRJoinParameter 
{
 protected:
  bool parameter_file;

  SDGString match_filename;        //!< file name of matches to treat
  SDGString path_filename;
  SDGString bank_name;        //!< file name of data bank
  SDGString bankQuery_name;   //!< file name of data query
  double gap_pen;            //!< gap penalty to join 2 HSP
  double dist_pen;            //!< distance penalty to join 2 HSP
  unsigned overlap;           
  double eval_filter;               //!< Threshold of e-value
  double id_filter;               //!< Threshold of identity
  unsigned len_filter;               //!< Threshold of length
  SDGString parameter_filename;  // Name of parameter file
  SDGString prefix_filename;  // Name of prefix
  bool join_frag;
  int id_param;
  bool time_stamp;
  double idtol;
  bool loadPath;
  int nbthread; // number of thread to use

 public:
  //! -Constructor
  BLRJoinParameter(void){reset();}; 
  
  void reset(void)
    {
      path_filename ="<not set>";
      bank_name="<not set>";
      bankQuery_name="<not set>";
      join_frag=false;
      gap_pen=0.05;
      dist_pen=0.2;
      overlap=20;
      eval_filter=1e-10;
      id_filter=0.0;
      len_filter=20;
      parameter_filename="<not set>";
      prefix_filename="<not set>";
      time_t t;
      t=time(&t);
      id_param=((int)t)-1043322138;
      time_stamp=true;
      idtol=2.0;
      loadPath=false;
      nbthread=1;
    };

  //! -Destructor
  virtual ~BLRJoinParameter(void){};
  
  
 public:

  /*! Function to parse the commande line 
    \param numarg: argument number must to parse
    \param tabarg: aray of arguments
    \return booleen true if parse is all right*/
  void parseOptArg (int numarg, char *tabarg[]);

  /*! Function to get parameter with interactive menu*/
  void help (void);

  void write(const SDGString& filename) const;

  void view(std::ostream& out) const;

  void setAlignFile(SDGString filename) {match_filename=filename;};

  /*! Function to get file name which use as subject by blast*/

  SDGString getBank(void) const              {return bank_name;};  

  /*! Function to get file name which use as query by blast */

  SDGString getQuery(void)  const            {return bankQuery_name;};

  SDGString getParameterFileName(void) const     {return parameter_filename;};
  SDGString getMatchFileName(void) const     {return match_filename;};
  SDGString getPrefixFileName(void) const     {return prefix_filename;};

  int getIdParam(void) const {return id_param;};

  /*! Function to get filter value used to build group
    \return filter \sa blaster.cpp*/ 
  double getEvalFilter(void) const           {return eval_filter;};
  double getIdFilter(void) const           {return id_filter;};
  unsigned getLenFilter(void) const           {return len_filter;};
  double getIdTolerance(void) const {return idtol;};

  /*! Function to set filter value used to build group
     \sa blaster.cpp*/ 
  void setEvalFilter(double f)                 {eval_filter=f;};
  void setIdFilter(double f)                 {id_filter=f;};
  void setLenFilter(unsigned f)                 {len_filter=f;};
  
  bool getJoin_frag(void) const             {return join_frag;};
  void setJoin_frag(bool join)             {join_frag=join;};
  double getGap_pen(void) const             {return gap_pen;};
  double getDist_pen(void) const             {return dist_pen;};
  unsigned   getOverlap(void) const             {return overlap;};
  void setMatch_filename(SDGString  filename) {match_filename=filename;};
  void setPrefix_filename(SDGString filename) {prefix_filename=filename;};
  void setPath_filename(SDGString filename){path_filename = filename;};
  SDGString getPath_filename (void) const {return path_filename;};
  bool getLoad_path (void) const {return loadPath;};
  void setLoad_path (bool load){loadPath=load;};
  int getNbThread(void) const {return nbthread;};
  void setNbThread(int n){nbthread=n;};
};
const static BLRJoinParameter defaultBLRJoinParameters;
#endif
