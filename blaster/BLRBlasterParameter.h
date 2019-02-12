/**
 * \file BLRBlasterParameter.h
 * \brief Header file for the class BLRBlasterParameter
 */

#ifndef BLRBLASTERPARAMETER_H
#define BLRBLASTERPARAMETER_H

#include <string.h>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include <fstream>
#include <list>
#include "SDGString.h"
#include "Cutter.h"

const int dim=2;             //!< dimension of array-type

/**
 * \class BLRBBlasterParameter
 * \brief Parameters for the BLASTER program
 */
class BLRBlasterParameter
{
 private:
  bool parameter_file;
  int verbose;

  SDGString bank_name;        //!< file name of data bank
  SDGString bankQuery_name;   //!< file name of data query
  SDGString bank_cut;         //!< file name of cut-out data bank
  SDGString bankQuery_cut;    //!< file name of cut-out data query
  SDGString type_blast;       //!< blast program can be run
  SDGString option_blast;     //!< option of blast command line
  unsigned long numdo;        //!< number of blast will do
  unsigned short nb_process;  //!< number of concurrent blast processes
  Cutter cutter;           //!< parameter of bank cut-out
  double eval_filter;               //!< Threshold of e-value
  double id_filter;               //!< Threshold of identity
  unsigned len_filter;               //!< Threshold of length
  SDGString parameter_filename;  // Name of blaster parameter file
  SDGString blaster_filename;    // Name of blaster files
  unsigned sensitivity_lvl;
  bool restart;
  bool first_run;
  int id_param;
  bool time_stamp;
  bool all_by_all;
  bool is_ncbiBlast;
  bool is_wuBlast;
  bool is_ncbiBlastPlus;
  bool prepare;
  bool cleanTmpFiles;

 public:
  //! -Constructor
  BLRBlasterParameter(void){reset();};

  BLRBlasterParameter(const SDGString& filename)
    {
      reset();
      load(filename);
    };

  void reset(void)
    {
      bank_name="<not set>";
      bankQuery_name="<not set>";
      bank_cut="<not set>";
      bankQuery_cut="<not set>";
      type_blast="blastn";
      option_blast="";
      numdo=0;
      nb_process=1;
      sensitivity_lvl=0;
      eval_filter=1e-10;
      id_filter=0.0;
      len_filter=20;
      restart=false;
      first_run=true;
      parameter_filename="<not set>";
      blaster_filename="<not set>";
      parameter_file=false;
      cutter.reset();
      time_t t;
      t=time(&t);
      id_param=((int)t)-1043322138;
      time_stamp=true;
      all_by_all=false;
      is_wuBlast=false;
      is_ncbiBlastPlus=false;
      is_ncbiBlast=false;
      prepare=false;
      cleanTmpFiles=false;
      verbose=0;
    };

  //! -Destructor
  virtual ~BLRBlasterParameter(void){};

 private:
  /*! Function to parse the command line
    \param numarg: argument number must to parse
    \param tabarg: array of arguments
    \return boolean true if parse is all right*/
  void parseOptArg (int numarg, char *tabarg[]);

  /*! Function to get parameter with interactive menu*/
  void help (void);
  void loadPrevious(BLRBlasterParameter& last);

  /*! Function to count the blast what will be do and check file
    \return mumber of blast shall do it*/
  int blastcount (void);


  /*! overload operator >>  */
  friend std::istream& operator>>(std::istream&, BLRBlasterParameter&);
  /*! overload operator <<  */
  friend std::ostream& operator<<(std::ostream&, const BLRBlasterParameter&);

  /*! overload operator == */
  friend bool operator==(const BLRBlasterParameter&,const BLRBlasterParameter&);
  /*! overload operator != */
  friend bool operator!=(const BLRBlasterParameter&,const BLRBlasterParameter&);

  //! Function to set length in Cutter to cut-out
  void setLength(unsigned long l)              {cutter.setLength(l);};

  //! Function to set over in Cutter to cut-out
  void setOver(unsigned short o)                {cutter.setOver(o);};

  //! Function to set number of maximum unknown subsequence length in Cutter to cut-out
  void setWord(unsigned short w)                {cutter.setWord(w);};

  //! set extension file of cut-out in Cutter
  void setExtention(SDGString e)     {cutter.setExtention(e);};

 public:

  void write(const SDGString& filename) const;

  void view(std::ostream& out) const;

  void load(const SDGString& filename)
    {
      if(filename.afterlast(".")!="param")
	{
	  std::string msg="BLRBlasterParameter::load():";
	  msg=filename+" with no suffix '.param' !";
	  throw SDGException(this,msg);
	}
      std::ifstream fin(filename);
      if(!fin)
	throw SDGException(this,SDGString("BLRBlasterParameter::load():")
			   +filename+" cannot be open!");
      fin>>*this;
      parameter_filename=filename;
    };

  void save(void)
    {
      if(time_stamp)
	{
	  std::ofstream fout("last_time_stamp.log");
	  fout<<id_param;
	  fout.close();
	  parameter_filename=blaster_filename
	    +"-"+SDGString(id_param)+".param";
	}
      else
	parameter_filename=blaster_filename+".param";

      std::ofstream fout(parameter_filename);
      if(!fout)
	throw SDGException(this,SDGString("BLRBlasterParameter::save():")
			   +parameter_filename+" cannot be open!");
      view(fout);
    };

  //! Function to check the parameter
  /*!\param numarg: argument number must to parse
    \param tabarg: array of arguments
    \return boolean true if parse is all right */
  bool start (int numarg, char* tabarg[]);

  /*! Function to get file name which use as subject by blast
    \return bank_cut */
  SDGString getBankCut(void) const              {return bank_cut;};

  SDGString getBank(void) const              {return bank_name;};

  /*! Function to get file name which use as query by blast
    \return bankQuery_cut*/

  SDGString getQueryCut(void)  const            {return bankQuery_cut;};

  SDGString getQuery(void)  const            {return bankQuery_name;};
  /*! Function to get file name which use by blaster
    \return blaster_filename*/

  bool getAll_by_all(void) const {return all_by_all;};

  SDGString getBlasterFileName(void) const     {return blaster_filename;};

  SDGString getParameterFileName(void) const     {return parameter_filename;};

  int getIdParam(void) const {return id_param;};

  /*! Function to get type of blast will be use
    \return type_blast*/
  SDGString getType(void) const             {return type_blast;};

  /*! Function to get blast options will be use
    \return option_blast*/
  SDGString getOption(void) const            {return option_blast;};

  /*! Function to get number blast will be do
    \return numdo*/
  unsigned long getNumdo(void) const                  {return numdo;};

  /*! Function to get number of concurrent blast process
    \return nb_process*/
  unsigned short getProcess(void) const       {return nb_process;};

  /*! Function to get length of cut-out
    \return length \sa Cutter*/
  unsigned long  getLength(void) const        {return cutter.getLength();};

  /*! Function to get over of cut-out
    \return over \sa Cutter*/
  unsigned short getOver(void) const         {return cutter.getOver();};

  /*! Function to get word of cut-out
    \return word \sa Cutter*/
  unsigned short getWord(void) const        {return cutter.getWord();};

  /*! Function to get extention file of cut-out
    \return extention \sa Cutter*/
  SDGString getCutExtention(void) const      {return cutter.getExtention();};

  unsigned getSensitivityLvl(void) const           {return sensitivity_lvl;};
  /*! Function to get filter value used to build group
    \return filter \sa blaster.cpp*/
  double getEvalFilter(void) const           {return eval_filter;};
  double getIdFilter(void) const           {return id_filter;};
  unsigned getLenFilter(void) const           {return len_filter;};

  /*! Function to set filter value used to build group
     \sa blaster.cpp*/
  void setEvalFilter(double f)                 {eval_filter=f;};
  void setIdFilter(double f)                 {id_filter=f;};
  void setLenFilter(unsigned f)                 {len_filter=f;};

  bool get_is_wuBlast(void) {return is_wuBlast;};
  bool get_is_ncbiBlast(void) {return is_ncbiBlast;};
  bool get_is_ncbiBlastPlus(void) {return is_ncbiBlastPlus;};
  bool getPrepare(void) {return prepare;};
  bool getCleanTmpFiles(void) {return cleanTmpFiles;};
  int getVerbose(void) {return verbose;};

};

#endif
