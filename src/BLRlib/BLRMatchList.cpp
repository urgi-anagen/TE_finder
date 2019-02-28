/**
 *
 * BLRMatchList.cpp
 *
 */

#include <BLRMatchList.h>

BLRMatchList::~BLRMatchList()
{
//! -- destructor --
  clear();
}

void BLRMatchList::clear()
{
//! -- set empty match_list  --
  match_list.clear();
}

void BLRMatchList::remove(BLRBlasterParameter para)
{
  //  ------------------------------------------------- 
  //! - Function to remove the matches that do not need
  /*  - Arguments: length is length of cut-out
      -            over is the overlap of cut-out
      -            same_bank_name is boolen, true if bank and querybank are same
      ---------------------------------------------------*/
  bool same_bank_name;
  if (para.getBankCut()==para.getQueryCut())
    same_bank_name=true;
  else
    same_bank_name=false;

  if(!same_bank_name) return;
  std::list<BlastMatch>::iterator iter_list=match_list.begin();  /*iterator to build list of matches*/
  
  while (iter_list!=match_list.end())
    {
	  bool remove=false;
	  if ( iter_list->getE_value() > para.getEvalFilter() 
	       || iter_list->getIdentity() < para.getIdFilter()
	       || iter_list->getLength() < para.getLenFilter()
	       )
	  {
//		  std::cout<<"remove 1"<<std::endl;
		  remove=true;
	  }
	  // Removing alignement of the same query/subject 
	  else if((iter_list->getQuery_num()==iter_list->getSubject_num())
	     &&(iter_list->getQuery_strt()==iter_list->getSubject_strt())
	     &&(iter_list->getQuery_end()==iter_list->getSubject_end()))
	  {
//		  std::cout<<"remove 2"<<std::endl;
		  remove=true;
	  }
	  else if(
			   ( ( iter_list->getQuery_num() == iter_list->getSubject_num() - 1 )
			     &&
			     ( iter_list->getQuery_strt() == (para.getLength()-para.getOver()+1)
				   && iter_list->getQuery_end() == para.getLength() )
			     &&
			     ( iter_list->getSubject_strt() == 1
				   && iter_list->getSubject_end() == para.getOver() )
			     &&
			     ( para.getLength()>0 && para.getOver()>0 ) )
		       ||
		       ( ( iter_list->getQuery_num()-1 == iter_list->getSubject_num() )
			     &&
			     ( iter_list->getSubject_strt() == (para.getLength()-para.getOver()+1)
			       && iter_list->getSubject_end() == para.getLength() )
			     &&
			     ( iter_list->getQuery_strt()==1
			       && iter_list->getQuery_end()==para.getOver() )
			     &&
			     ( para.getLength()>0 && para.getOver()>0 ) )
		     )
	  {
//		  std::cout<<"remove 3"<<std::endl;
		  remove=true;
	  }
	  if(remove)
	    iter_list=match_list.erase(iter_list);
	  else
	    iter_list++;
    }
}

void BLRMatchList::fill_ncbilist(SDGString& filename)
{
  //  --------------------------------------------------- 
  //! - Function to fill matchlist with alignments found 
  /*   ---------------------------------------------------*/
  std::ifstream match_file(filename);           /*file where input matchList*/
  if (match_file.bad())
    {
      std::cerr<<"ERROR: "<<filename<<" could not be open!"<<std::endl;
      exit(EXIT_FAILURE);
    }
  readncbiblastfield(match_file);
  match_file.close();
}

void BLRMatchList::readncbiblastfield(std::istream& in)
{
  bool empty=true;
  while(in)
    {
      try
	{
	  match.readncbiblastfield(in);
	}    
      catch(BlastMatch::Empty_Parser_output e)
	{
	  if(empty)
	    std::cout<<"empty output"<<std::endl;
	  break;
	}
      empty=false;
      match_list.push_back(match);
    }
}

void BLRMatchList::fill_wulist(SDGString& filename)
{
  //  --------------------------------------------------- 
  //! - Function to fill matchlist with alignments found 
  /*   ---------------------------------------------------*/
  std::ifstream match_file(filename);           /*file where input matchList*/
  if (match_file.bad())
    {
      std::cerr<<"ERROR: "<<filename<<" could not be open!"<<std::endl;
      exit(EXIT_FAILURE);
    }
  readwublastfield(match_file);
  match_file.close();
}

void BLRMatchList::readwublastfield(std::istream& in)
{
  bool empty=true;
  while(in)
    {
      try
	{
	  match.readwublastfield(in);
	}    
      catch(BlastMatch::Empty_Parser_output e)
	{
	  if(empty)
	    std::cout<<"empty output"<<std::endl;
	  break;
	}
      empty=false;
      match_list.push_back(match);
    }
}

void BLRMatchList::show_list()
{
  // --------------------------------------------------
  //!- Function for developement to show the matchlist
  /* --------------------------------------------------*/
  std::list<BlastMatch>::iterator id=match_list.begin();
  std::cout<<"-----------------"<<std::endl; 
  while(id!=match_list.end())
  {
    id->view();
    id++;
  }
}

void BLRMatchList::save_list(SDGString name, int flag)
{
  // -----------------------------------------------------
  /*!- Function to save the matchlist
    /param  name is file to save
    /param flag is boolean true if is a new run of blast*/
  //-----------------------------------------------------
  std::ofstream file;
  if (flag)
    // to write in file end
    file.open(name.start(), std::ios::app);
  else
    // to write in file begin
    file.open(name.start(), std::ios::trunc);
  if (file.bad())
    {
      std::cout<<"ERROR: "<<name.start()<<" could not open"<<std::endl;
      return;
    }

  savebin(file);

  file.close();
}

//! fonction to save the hash_align in binary file
void BLRMatchList::savebin(std::ostream& out)
{
  std::list<BlastMatch>::iterator iter_list=match_list.begin();  //iterator to visit match_list
   while (iter_list!=match_list.end())
    {
      iter_list->writelst(out);
      iter_list++;
    }
}

//! fonction to write the hash_align in ascii format
void BLRMatchList::writetxt(std::ostream& out)
{
  std::list<BlastMatch>::iterator iter_list=match_list.begin();  //iterator to visit match_list
  while (iter_list!=match_list.end())
    {
      out<<*iter_list;
      iter_list++;
    }
}

//! fonction to read the hash_align in ascii format
void BLRMatchList::readtxt(std::istream& in)
{
  BlastMatch match;            // current alignment 
  while (in)
    {
      in>>match;
      if(match.getQuery_num()==0) break;
      match_list.push_back(match);
    }
}
void BLRMatchList::readtxt(char * data)
{
  std::string sdata=data;
  std::istringstream sin(sdata);
  BlastMatch match;            // current alignment 
  while (sin)
    {
      sin>>match;
      if(match.getQuery_num()==0) break;
      match_list.push_back(match);
    }
}
