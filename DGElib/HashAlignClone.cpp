#include "HashAlignClone.h"

HashAlignClone::HashAlignClone(){
}

HashAlignClone::~HashAlignClone(){
}


void HashAlignClone::search(unsigned word_size)
{
  //bool b1 = sequence1.getDE() == "BlastclustCluster2Mb22_chunk1617 (dbseq-nr 1) [174895,179677]";
  //bool b2 = sequence2.getDE() == "BlastclustCluster2Mb35_chunk2038 (dbseq-nr 2) [84935,89713]";

  //bool verbose = b1 && b2;
  //if (verbose)
	//std::cout<<"F O U N D"<<std::endl<<std::flush;


  SDGString seq1=sequence1.toString();
  SDGString seq2=sequence2.toString();
  unsigned len1=seq1.length();
  unsigned len2=seq2.length();
  std::map<SDGString,std::list<unsigned> > map_word;
  //std::map<int,DiagClone> diag_map;
   
  std::map<int,std::list<Range> > diag_map;
  //if (verbose)
	//std::cout<<"L O O P 1"<<std::endl<<std::flush;
 
  for(unsigned i=0;i<len1-word_size+1;i++)
    {
      SDGString word=seq1.substr(i,word_size);
      unsigned pos=word.rfind('N');
      if(pos>=word_size)
	  map_word[word].push_back(i);
    }

  //if (verbose)
	//std::cout<<"L O O P 2"<<std::endl<<std::flush;
 
  for(unsigned i=0;i<len2-word_size+1;i++)
    {
      SDGString word=seq2.substr(i,word_size);
      unsigned pos=word.rfind('N');
      if(pos>=word_size)
	{
	  std::list<unsigned>& l=map_word[word];
	  for(std::list<unsigned>::iterator it=l.begin();
	      it!=l.end();it++)
	    {
	      //diag_map[i-(*it)].add(Range((*it),(*it)+word_size-1));
	      Range rToAdd = Range((*it),(*it)+word_size-1);
	      addToDiagMap(rToAdd, i-(*it), diag_map);
	    }
	}
    }
  
  sort_frag.clear();
   
  //if (verbose)
	//std::cout<<"L O O P 3"<<std::endl<<std::flush;
  
  //for(std::map<int,DiagClone>::iterator i=diag_map.begin();
      //i!=diag_map.end();i++)
  for(std::map<int,std::list<Range> >::iterator i=diag_map.begin();
      i!=diag_map.end();i++)

    {
      for(DiagClone::iterator j=i->second.begin();
	  j!=i->second.end();j++)
	{
	  RangeAlign r1(0,j->getStart(),j->getEnd());
	  RangeAlign r2(0,i->first+j->getStart(),
			i->first+j->getEnd());

	  RangePair rp(r1,r2);
	  rp.setScore(r1.getLength());
	  sort_frag.insert(std::pair<unsigned long,RangePair>(rp.getScore(),rp));
	}
    }
  //if (verbose)
	//std::cout<<"E X I T"<<std::endl<<std::flush;
 
}


void HashAlignClone::addToDiagMap(const Range& r, int key, std::map<int, std::list<Range> >& diag_map)
{
    if(diag_map[key].back().getEnd()==r.getEnd()-1)
	{
	  diag_map[key].back().setEnd(r.getEnd());
	}
      else
	{
	  if(diag_map[key].back().getStart()==r.getStart()+1)
	    {
	      diag_map[key].back().setStart(r.getStart());
	    }
	  else
	    {      	
	      diag_map[key].push_back(r);
	    }
	}
    
}



