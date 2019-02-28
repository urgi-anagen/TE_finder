/**
 *
 * BLRGroup.cpp
 *
 **/

#include <map>

#include "BLRGroup.h"

/*--------------------------------------------
  - Function add  Member in group
  - Arguments: gi-> reference on current group
  -            mem-> reference on Member object
  --------------------------------------------*/
void BLRGroup::addMember(GROUPLIST::iterator& gi,Member& mem)
{
  mem.id=count_memb;
  vec_memb.push_back(mem);
  gi->insert(gi->begin(),mem.id);

  memIdx->insert_member_idx(mem,gi);
  count_memb++;

  /* debug */
   //std::cout<<"Members #:"<<count_memb<<" "<<vec_memb.size()<<std::endl;
   //std::cout<<"Members id";
   //for(unsigned m=0;m<count_memb;m++)
   //  std::cout<<":"<<m<<"-"<<vec_memb[m].id<<
   //    std::flush<<"("<<vec_idx[m]->id_member<<")"<<std::flush;
   //std::cout<<std::endl;

}

/*-----------------------------------------------------
  - Function add group in group list
  - Arguments: memQ-> Member created by query alignment
  -            memS-> Member created by subject alignment
  -------------------------------------------------------*/
void BLRGroup::addGroup(Member& memQ,Member& memS)
{
  std::list<unsigned> g;
  GROUPLIST::iterator group_it=group_list.insert(group_list.begin(),g);
  addMember(group_it,memQ);
  addMember(group_it,memS);
}

/*--------------------------------------------------------
  - Function merge two group
  - Arguments: iQ-> current group find by query alignment
  -            iS-> current group find by subject alignment
  --------------------------------------------------------*/
void BLRGroup::mergeGroup(GROUPLIST::iterator& iQ,GROUPLIST::iterator& iS)
{
  for(std::list<unsigned>::iterator i=iS->begin(); i!=iS->end(); i++)
    memIdx->set_groupit_member_idx(*i,iQ);
  iQ->splice(iQ->begin(),*iS,iS->begin(),iS->end());
  group_list.erase(iS);
}

/*---------------------------------------------------------------------
  - Function merge two Member
  - Arguments: a-> reference on current alignment
  -            b-> reference on current  Member
  - c-> reference on current complement alignment a.S<-a.Q || a.Q<-a.S.
  ----------------------------------------------------------------------*/
void BLRGroup::mergeMember(GROUPLIST::iterator& ga,unsigned ma,
			   GROUPLIST::iterator& gb,unsigned mb,
			   GROUPLIST::iterator& gc )
{
  //std::cout<<"mergeMembers:"<<std::endl
  //	   <<vec_memb[ma]<<std::endl<<vec_memb[mb]<<std::endl;

  // Case when c is assigned but not a
  // roam_group must assign a group at current Member,
  // function Compare must manage this:
  if (ga==group_list.end())
    {
      std::cout<<"###VERIFIER COMPARE QUAND Si==group_list.end()###";
      return;
    }
  // merge if a and b is in different group
  if(ga!=gb)
    {
      for(std::list<unsigned>::iterator i=gb->begin(); i!=gb->end(); i++)
	memIdx->set_groupit_member_idx(*i,ga);
      ga->splice(ga->begin(),*gb,gb->begin(),gb->end());
      if(gc==gb) gc=ga; // To change gc before erase gb
      group_list.erase(gb);
      gb=ga;
    }
  // Normally case where current Member (a) or/and this complement (c) found group before
  // merge if a and b is in the same group
  if (ga==gb && ma!=mb)
    {
      vec_memb[ma].merge(vec_memb[mb]);
      vec_memb[ma].idlist.splice(vec_memb[ma].idlist.end(),vec_memb[mb].idlist);
      if(vec_memb[ma].getMin()>=memIdx->get_start_member_idx(ma)
	 && vec_memb[ma].getMax()<=memIdx->get_end_member_idx(ma))
	return;
      memIdx->adjust_member_idx(ma,vec_memb[ma].getMin(),vec_memb[ma].getMax());
      memIdx->erase_member_idx(mb);
      ga->remove(mb);
      mb=ma;
    }
  //std::cout<<"\tmerged:"<<vec_memb[ma]<<std::endl;
}

/*--------------------------------------------------------
  - Function merge a RangeAlignSet to an existing Member
  - Arguments: align-> RangeAlignSet instance to merge
  -            member_id-> member id
  --------------------------------------------------------*/
void BLRGroup::mergeAlign( RangeAlignSet& align,long long id, unsigned member_id)
{
  //std::cout<<"mergeAlign"<<std::endl;
  //std::cout<<vec_memb[member_id]<<std::endl;
  //std::cout<<align<<std::endl;
  vec_memb[member_id].idlist.push_back(id);
  vec_memb[member_id].merge(align);

  if(align.getMin()>=memIdx->get_start_member_idx(member_id)
     && align.getMax()<=memIdx->get_end_member_idx(member_id))
    return;
  memIdx->adjust_member_idx(member_id,vec_memb[member_id].getMin(),vec_memb[member_id].getMax());
  //std::cout<<"\tmerged"<<*(vec_idx[member_id])<<std::endl;
}

/*----------------------------------------------------
  - Function inverse strand of all Member of this group
  - Argument: g-> group to inverse
  -----------------------------------------------------*/
void BLRGroup::Inverse(GROUPLIST::iterator& g)
{
  for(std::list<unsigned>::iterator member_iter=g->begin();
      member_iter!=g->end();member_iter++)
    vec_memb[*member_iter].reverse();
}

/*----------------------------------------------------
  - Function compare Member whith Query alignment
  - Arguments: ai->current alignment
  -            mi->current Member
  -            Qi,Si->group assigne to Q and S
  -            Gi->current group
  - Return bool true if Qi and mi are same
  ------------------------------------------------------*/
bool BLRGroup::compareQ(RangeAlignSet& alignQ,RangeAlignSet& alignS,
			long long id,
			unsigned mi,
			GROUPLIST::iterator& Qi,GROUPLIST::iterator& Si,
			GROUPLIST::iterator& Gi)
{
  id=id*(-1);

  // same num ?
  if (alignQ.getNumChr()!=vec_memb[mi].getNumChr()) return false;

  // fit member to Query if there are overlap, and modifie strand if necesary.
  if( vec_memb[mi].overlap(alignQ) )/*test for sequence include into Member*/
    {
      double lenQ=vec_memb[mi].overlap_length(alignQ);
      unsigned lenM=std::max(vec_memb[mi].getLengthSet(),alignQ.getLengthSet());

      double coverQ=lenQ/lenM;
      //std::cout<<"coverQ:"<<coverQ<<std::endl;
      if(coverQ<cover_limit) return false;
      if (vec_memb[mi].isPlusStrand()!=alignQ.isPlusStrand())
	{
	  if (Qi==group_list.end() && Si==group_list.end())
	    {
	      // alignment did not find a group before
	      // could change strand of subject and  query.
	      alignQ.reverse();
	      alignS.reverse();
	    }
	  else if (Si!=Gi && Qi!=Gi)
	    {
	      // Query or subject find another group before !!
	      // to push it into this group, strand of other group's Member
	      // must be inversed

	      Inverse(Gi);
	    }
	}
      mergeAlign(alignQ,id,mi);
      return true;
   }
  return false;
}

bool BLRGroup::compareS(RangeAlignSet& alignQ,RangeAlignSet& alignS,
			long long id,
	 		unsigned mi,
			GROUPLIST::iterator& Qi,GROUPLIST::iterator& Si,
			GROUPLIST::iterator& Gi)
{
  // same num ?
  if (alignS.getNumChr()!=vec_memb[mi].getNumChr()) return false;

  // fit Member to Query if there are overlap, and modifie strand if necesary.
  if( vec_memb[mi].overlap(alignS) )/*test for sequence include into Member*/
    {
      double lenS=vec_memb[mi].overlap_length(alignS);
      unsigned lenM=std::max(vec_memb[mi].getLengthSet(),alignS.getLengthSet());

      double coverS=lenS/lenM;
      //std::cout<<"coverS:"<<coverS<<std::endl;
      if(coverS<cover_limit) return false;
      if (vec_memb[mi].isPlusStrand()!=alignS.isPlusStrand())
	{
	  if (Qi==group_list.end() && Si==group_list.end())
	    {
	      // alignment did not find a group before
	      // could change strand of subject and  query.
	      alignS.reverse();
	      alignQ.reverse();
	    }
	  else if (Si!=Gi && Qi!=Gi)
	    {
	      // Query or subject find another group before !!
	      // to push it into this group, strand of other group's Member
	      // must be inversed

	      Inverse(Gi);
	    }
	}
      mergeAlign(alignS,id,mi);
      return true;
   }
  return false;
}

/*----------------------------
  - Function all over groups -
  - Parameter: alignment     -
  ----------------------------*/
void BLRGroup::roam_group( RangePairSet& align, int verbose )
{
  long long id=align.getId();

  RangeAlignSet alignQ(align.getRangeAlignSetQ());
  RangeAlignSet alignS(align.getRangeAlignSetS());
  
  GROUPLIST::iterator group_iter=group_list.end();
  GROUPLIST::iterator group_iterQ=group_list.end();
  GROUPLIST::iterator group_iterS=group_list.end();

  unsigned member_idQ=0,member_idS=0 ;
  if(verbose>0)
  	  std::cout<<alignS<<std::endl<<std::flush;

  std::vector<unsigned> idQ=memIdx->search_member(alignQ.getNumChr(),
					  alignQ.getMin(),
					  alignQ.getMax());
  for(std::vector<unsigned>::iterator i=idQ.begin();
      i!=idQ.end();i++)
    {
      unsigned member_id=*i;
      group_iter=memIdx->get_groupit_member_idx(member_id);
      if (compareQ(alignQ,alignS,id,member_id,group_iterQ,group_iterS,
		   group_iter))
	{
	  // alignment query and member overlap
	  if(group_iterQ==group_list.end())
	    {
	      // first group found
	      group_iterQ=group_iter;
	      member_idQ=member_id;
	    }
	  else
	    {
	      // query found a member before
	      mergeMember(group_iterQ,member_idQ,group_iter,
				  member_id,group_iterS);
	    }
	}
    }

  if(verbose>0)
	  std::cout<<alignS<<std::endl<<std::flush;
  std::vector<unsigned> idS=memIdx->search_member(alignS.getNumChr(),
					  alignS.getMin(),
					  alignS.getMax());


  for(std::vector<unsigned>::iterator i=idS.begin();
      i!=idS.end();i++)
    {
      unsigned member_id=*i;
      group_iter=memIdx->get_groupit_member_idx(member_id);
      if (compareS(alignQ,alignS,id,member_id,group_iterQ,group_iterS,
		   group_iter))
	{
	  // alignment query and member overlap
	  if(group_iterS==group_list.end())
	    {
	      // first group found
	      group_iterS=group_iter;
	      member_idS=member_id;
	    }
	  else
	    {
	      // query found a member before
	      mergeMember(group_iterS,member_idS,group_iter,
				  member_id,group_iterQ);
	    }
	}
    }
  build_group(alignQ,alignS,id,group_iterQ,group_iterS);
}

/*-----------------------------------------
  - Function to build group
  - Arguments: ai-> current alignment
  -            iQ-> group found by Query
  -            iS-> group found by Subject
  ------------------------------------------*/
void BLRGroup::build_group(RangeAlignSet& alignQ,RangeAlignSet& alignS,
			   long long id,
			   GROUPLIST::iterator& iQ, GROUPLIST::iterator& iS)
{
  Member mQ(alignQ,id*(-1));
  Member mS(alignS,id);

  if(iQ!=group_list.end() && iS==group_list.end())
    {
      addMember(iQ,mS);
    }
  else if(iQ==group_list.end() && iS!=group_list.end())
    {
      addMember(iS,mQ);
    }
  else if(iQ!=group_list.end() && iS!=group_list.end() && iQ!=iS)
    {
      mergeGroup(iQ,iS);
    }
  else if(iQ==group_list.end() && iS==group_list.end())
    {
      addGroup(mQ,mS);
    }
}

std::list<RangePairSet> BLRGroup::getRpsListAfterLoad( int verbose )
{
      unsigned long nbAlignments=0;      // number of alignment found and treated

      std::list<RangePairSet> rp_list;
     
      if (!grouper_parameter->getLoad_path()){

      	if(verbose>0)
    	  	std::cout<<"Load the matches..."<<std::endl<<std::flush;

      	match_map.load();

      	if(verbose>0)
    	  	std::cout<<"Matches were loaded."<<std::endl;

      	if(grouper_parameter->getJoin_frag())
      	{
    	  	if(verbose>0)
    		  	std::cout<<"Connect the fragments..."<<std::endl<<std::flush;
   	 	  	match_map.mapPath(true, false, false, false, 2);
    	  	if(verbose>0)
    		  	std::cout<<"Fragments were connected."<<std::endl;
      	}
      	else
    	  match_map.mapPath(false, false, false, false, 0);

      	if(verbose>0)
    	  	std::cout<<"Roam groups (coverage "<<grouper_parameter->getCoverage()<<")..."<<std::endl<<std::flush;

      	cover_limit=grouper_parameter->getCoverage();

      	for( BLRMatchMap::MapPath::iterator m=match_map.path_begin();
      	m!=match_map.path_end(); m++ )
      	{
    	  	while(!m->second.empty())
    	  	{
    		  	RangePairSet rp=m->second.back();
    		  	m->second.pop_back();
    		  	rp_list.push_back(rp);
    	  	}
      	}
      }else{
     	if(verbose>0)
    	  	std::cout<<"Load the matches from path file..."<<std::endl<<std::flush;
	match_map.loadPath();
	rp_list = match_map.getRpsList();
      }

      memIdx=new BLRMemIdx(match_map.getNbQseq());
      match_map.clear();
      
      rp_list.sort(RangePair::greaterLengthQ);
      return rp_list;
}

/*---------------------------------
  - Function to build the groups
  - Argument: verbose
  ----------------------------------*/
void BLRGroup::group( int verbose )
{
      unsigned long nbAlignments=0;      // number of alignment found and treated

      std::list<RangePairSet> rp_list;
     
      if (!grouper_parameter->getLoad_path()){

      	if(verbose>0)
    	  	std::cout<<"Load the matches..."<<std::endl<<std::flush;

      	match_map.load();

      	if(verbose>0)
    	  	std::cout<<"Matches were loaded."<<std::endl;

      	if(grouper_parameter->getJoin_frag())
      	{
    	  	if(verbose>0)
    		  	std::cout<<"Connect the fragments..."<<std::endl<<std::flush;
   	        
		match_map.mapPath(true, false, false, false, verbose-1);
			if (verbose>0)
				match_map.writePath(grouper_parameter->getPrefixFileName() +  ".gpath");
		if(verbose>0)
    		  	std::cout<<"Fragments were connected."<<std::endl;
      	}
      	else
    	  match_map.mapPath(false, false, false, false, 0);
      	
	unsigned long id = 1;
 	for( BLRMatchMap::MapPath::iterator m=match_map.path_begin();
      	m!=match_map.path_end(); m++ )
      	{
    	  	while(!m->second.empty())
    	  	{
    		  	RangePairSet rp=m->second.front();
    		  	m->second.pop_front();
			rp.setId(id);
    		  	rp_list.push_back(rp);
			id++;
    	  	}
      	}
      }else{
     	if(verbose>0)
    	  	std::cout<<"Load the matches from path file..."<<std::endl<<std::flush;
	//match_map.loadPath();
	match_map.loadPath(grouper_parameter->getPath_filename(),1);
	rp_list = match_map.getRpsList();
      }

      memIdx=new BLRMemIdx(match_map.getNbQseq());
      match_map.clear();
  
      if(verbose>0)
    	  	std::cout<<"Roam groups (coverage "<<grouper_parameter->getCoverage()<<")..."<<std::endl<<std::flush;

      cover_limit=grouper_parameter->getCoverage();
	
      
	      
      rp_list.sort(RangePair::greaterLengthQ);
   
      if (verbose > 0){
      	if (!grouper_parameter->getLoad_path())
      		writePath(grouper_parameter->getPrefixFileName() + ".rpsListNoLoadPath",rp_list,1);
      	else
      		writePath(grouper_parameter->getPrefixFileName() + ".rpsListLoadPath",rp_list,1);
      } 	

      std::list<RangePairSet>::iterator r=rp_list.begin();
      while(r!=rp_list.end())
      {
    	  nbAlignments++;
    	  r->setId(nbAlignments);
    	  roam_group( *r, verbose-1 );
    	  rp_list.pop_front();
    	  r=rp_list.begin();
      }
     
      if(verbose>0)
      {
    	  std::cout<<"number of alignments: "<<nbAlignments<<std::endl;
    	  //std::cout<<"number of members: "<<count_memb<<std::endl;    //overestimated (Tim, 19.05.2009)
    	  std::cout<<"number of members: "<<getNbMembers()<<std::endl;
	  std::cout<<"number of groups: "<<group_list.size()<<std::endl;
    	  std::cout<<"Roaming done."<<std::endl;
      }

      if(nbAlignments==0)
      {
    	  std::cout<<"Error when building groups: no alignment !"<<std::endl;
    	  exit(EXIT_FAILURE);
      }

}

void BLRGroup::writePath(const SDGString& filename,std::list<RangePairSet> l, int verbose)
{
  if(verbose>0)
    std::cout<<"writing 'path' file..."<<std::flush;

  std::ofstream fout(filename);
  std::ofstream foutRpsAttr(filename + ".attr");  

  unsigned path_id=0;
 
  for(std::list<RangePairSet>::iterator iter_list
	    =l.begin();iter_list!=l.end();
	  iter_list++)
	{
		unsigned long id = iter_list->getId();
		iter_list->write(fout,id,iter_list->first.getNameSeq(),iter_list->second.getNameSeq());
		iter_list->writeRpsAttr(foutRpsAttr,id,iter_list->first.getNameSeq(),iter_list->second.getNameSeq());
	}
    
  if(verbose>0)
    std::cout<<" done"<<std::endl;
}

/*---------------------------------
  - Function to print out group_list
  - Argument: out-> file of out
  ----------------------------------*/
void BLRGroup::show_group(std::ostream& out)
{
  int count_group=1;             //! count number of groups found
  int count_mb=1;                //! count number of members found
  SDGString chr_name;

  GROUPLIST::iterator group_iter=group_list.begin();

  bool samedb=false;
  if(grouper_parameter->getBank()==grouper_parameter->getQuery())
    samedb=true;

  while(group_iter!=group_list.end())
    {
      std::list<unsigned>::iterator member_iter=group_iter->begin();
      out<<"GROUP NUMBER "<<count_group<<" CLUSTER:"<<gr2clust[count_group]
	 <<" Contains "<<group_iter->size()<<" members"<<std::endl;
      out<<"\tmember \tSeq_num\tstart \tend  \tfrag \tlength \ttotal \talign_list"<<std::endl;
	  while(member_iter!=group_iter->end())
	    {
	      Member memb=vec_memb[*member_iter];
	      if (memb.idlist.front()>0) //subject
		if(samedb)
		  chr_name=match_map.numQ2name(memb.getNumChr());
		else
		  chr_name=match_map.numS2name(memb.getNumChr());
	      else  // query
		chr_name=match_map.numQ2name(memb.getNumChr());

	      out<<count_group<<"\t"<<count_mb<<"\t";
	      out<<chr_name<<"\t";
	      out<<memb.getStart()<<"\t";
	      out<<memb.getEnd()<<"\t";
	      out<<memb.getRangeSet().size()<<"\t";
	      out<<memb.getLengthSet()<<"\t";
	      out<<memb.idlist.size()<<"\t";
	      for(std::list<long long>::iterator i=memb.idlist.begin();
		  i!=memb.idlist.end();i++)
		out<<*i<<" ";
	      out<<std::endl;
	      count_mb++;
	      member_iter++;
	    }
	  count_group++;
	  group_iter++;
    }

  unsigned count=0;
  for(std::vector< std::vector<unsigned> >::iterator l=clust.begin();l!=clust.end();l++)
    {
      out<<"cluster #"<<++count<<" size="<<l->size()<<" :";
      for(std::vector<unsigned>::iterator i=l->begin();i!=l->end();i++)
	{
	  out<<" "<<*i;
	}
      out<<std::endl;
    }
  out<<"Members found="<<--count_mb<<std::endl;
  out<<"Groups found="<<--count_group<<std::endl;
  out<<"Clusters found="<<clust.size()<<std::endl;
}

/*----------------------------------------
  - Function to save groups in fasta format
  -----------------------------------------*/
void BLRGroup::save( int verbose )
{
  unsigned long count_gp;            // Number of group
  unsigned long count_mb;            // Number of group member
  unsigned long count_all;                 // Total number of members
  GROUPLIST::iterator group_iter;   // counter of loop to visit group_list
  std::list<unsigned>::iterator member_iter; // counter of loop to visit member of group

  RangeMap membmap;

  group_iter=group_list.begin();

  count_gp=1;
  count_all=0;

  bool samedb=false;
  if(grouper_parameter->getBank()==grouper_parameter->getQuery())
    samedb=true;

  // all over group
  while(group_iter!=group_list.end())
    {
	  if(verbose>0)
		  std::cout<<"Working on group number "<<count_gp<<":"<<std::endl;
      member_iter=group_iter->begin();
      count_mb=0;

      // all over member
      while(member_iter!=group_iter->end())
	{
	  count_mb++;
	  count_all++;
	  SDGString subseqname;
	  std::string chr_name;
	  Member memb=vec_memb[*member_iter];
	  if(verbose>0)
	  {
		  std::cout<<"\tsubseq from "<<memb.getNumChr()<<std::flush;
		  std::cout<<" "<<memb.getStart()<<".."<<memb.getEnd()
		  <<" length member:"<<memb.getLengthSet()<<std::endl;
	  }
	  if (memb.idlist.front()>0) //subject
	    {
	      if(samedb)
		chr_name=match_map.numQ2name(memb.getNumChr());
	      else
		chr_name=match_map.numS2name(memb.getNumChr());
	      subseqname="MbS"+SDGString(count_all)
		+"Gr"+SDGString(count_gp)
		+"Cl"+SDGString(gr2clust[count_gp]);
	    }
	  else  // query
	    {
	      chr_name=match_map.numQ2name(memb.getNumChr());
	      subseqname="MbQ"+SDGString(count_all)
		+"Gr"+SDGString(count_gp)
		+"Cl"+SDGString(gr2clust[count_gp]);
	    }
	  if(verbose>0)
	  {
		  std::cout<<"\tseqname:"<<chr_name<<std::endl;
		  std::cout<<"\tsubseqname:"<<subseqname<<std::endl;
	  }
	  membmap.add(RangeSeq(memb,subseqname,chr_name));
	  member_iter++;
	}
      count_gp++;
      group_iter++;
    }

  std::ostringstream filename;
  filename<<grouper_parameter->getPrefixFileName()<<".group.c"
	  <<grouper_parameter->getCoverage();

  if(verbose>0)
	  std::cout<<"Writing 'map' file..."<<std::flush;
  SDGString mapfile=filename.str()+".map";
  membmap.save(mapfile);
  if(verbose>0)
	  std::cout<<" done"<<std::endl;

  if(verbose>0)
	  std::cout<<"Writing 'set' file..."<<std::flush;
  SDGString setfile=filename.str()+".set";
  membmap.saveSet(setfile);
  if(verbose>0)
	  std::cout<<" done"<<std::endl;

  if(verbose>0)
	  std::cout<<"Writing 'fasta' file..."<<std::flush;
  SDGString seqfile=filename.str()+".fa";
  membmap.writeSeq(seqfile,grouper_parameter->getBank());
  if(grouper_parameter->getBank()!=grouper_parameter->getQuery())
    membmap.writeSeq(seqfile,grouper_parameter->getQuery());
  if(verbose>0)
	  std::cout<<" done"<<std::endl;

  if(verbose>0)
	  std::cout<<"Writing 'param' file..."<<std::flush;
  SDGString parafile=filename.str()+".param";
  std::ofstream fout(parafile);
  grouper_parameter->view(fout);
  if(verbose>0)
	  std::cout<<" done"<<std::endl;
}

/*----------------------------------------
  - Function to cluster the groups
  -----------------------------------------*/
void BLRGroup::cluster( int verbose )
{
  Graph<unsigned> graph;

  std::vector<unsigned> group_size;
  unsigned nb_group=0;
  for(GROUPLIST::iterator i=group_list.begin(); i!=group_list.end();
      i++,nb_group++)
  {
      group_size.push_back(i->size());
      nb_group++;
    };

  if(verbose>0)
	  std::cout<<"Building graph..."<<std::flush;
  // all over group
  unsigned long count_gp1=0;
  unsigned long count_gp2=0;
  unsigned long count_mb1=0;
  GROUPLIST::iterator group_iter1=group_list.begin();
  while(group_iter1!=group_list.end())
    {
      std::map<unsigned long ,unsigned > count;
      std::list<unsigned>::iterator member_iter1=group_iter1->begin();
      count_gp1++;
      // all over member
      while(member_iter1!=group_iter1->end())
	{
	  Member memb1=vec_memb[*member_iter1];
	  count_mb1++;
	  std::vector<unsigned> member_id
	    =memIdx->search_member(memb1.getNumChr(),memb1.getMin(), memb1.getMax());
	  for(std::vector<unsigned>::iterator i=member_id.begin();
	      i!=member_id.end(); i++)
	    {
	      Member memb2=vec_memb[*i];
	      int cover=0;
	      if(// member1 included in member2
		 memb1.isIncluded(memb2))
		{
		  cover=(unsigned)
		    memb1.overlap_length(memb2);
		}
	      else if(//member2 included in member1
		      memb1.isContained(memb2))
		{
		  cover=(unsigned)
		    memb1.overlap_length(memb2);
		}
	      else if(// member1 and member2 overlap and member1 first
		      memb1.overlap(memb2)
		      && memb1>memb2)
		{
		  cover=(unsigned)
		    memb1.overlap_length(memb2);
		}
	      else if(// member1 and member2 overlap and member 2 first
		      memb1.overlap(memb2)
		      && (memb1<memb2
			  || memb1==memb2))
		{
		  cover=(unsigned)
		    memb1.overlap_length(memb2);
		}
	      if( cover>=grouper_parameter->getGraphfilter())
		{
		  count_gp2=0;
		  for(GROUPLIST::iterator gi=group_list.begin();
		      gi!=group_list.end(); gi++)
		    {
		      count_gp2++;
		      if(memIdx->get_groupit_member_idx(memb2.id)==gi)
			break;
		    }

		  count[count_gp2]++;
		}
	    }
       	  member_iter1++;
	}
      graph.add_node(count_gp1);
      for(std::map<unsigned long,unsigned>::iterator i=count.begin();i!=count.end();i++)
	{
	  graph.add_edge(count_gp1,i->first);
	}
      group_iter1++;
    }
  if(verbose>0)
	  std::cout<<" done"<<std::endl;

  if(verbose>0)
	  std::cout<<"Search connected components..."<<std::flush;
  graph.connexComp(clust);
  if(verbose>0)
	  std::cout<<" done"<<std::endl;

  std::ostringstream clus_filename;
  clus_filename<<grouper_parameter->getPrefixFileName()<<".group.c"
	  <<grouper_parameter->getCoverage()<<".cluster.dot";

  graph.toDot(clus_filename.str());
  unsigned count=0;
  for(std::vector< std::vector<unsigned> >::iterator l=clust.begin();l!=clust.end();l++)
    {
      count++;
      for(std::vector<unsigned>::iterator i=l->begin();i!=l->end();i++)
	  gr2clust[*i]=count;
    }
}

//---------------------------------------------------------------------------
void BLRGroup::include_filter( int verbose )
{
	if(verbose>0)
		std::cout<<"Including filter..."<<std::flush;
	GROUPLIST::iterator group_iter1=group_list.begin();
	std::list<GROUPLIST::iterator> gr2erase;

	while(group_iter1!=group_list.end())
	{
		std::list<unsigned>::iterator member_iter1=group_iter1->begin();
		// all over member of a group
		while(member_iter1!=group_iter1->end())
		{
			Member memb1=vec_memb[*member_iter1];
			std::vector<unsigned> member_id
			=memIdx->search_member(memb1.getNumChr(),memb1.getMin(), memb1.getMax());
			for(std::vector<unsigned>::iterator i=member_id.begin();
			i!=member_id.end(); i++)
			{
				if(*i!=*member_iter1)
				{
					Member memb2=vec_memb[*i];
					if(//member2 included in member1
							memb1.isIncluded(memb2))
					{
						if(std::find(gr2erase.begin(),gr2erase.end(),
								memIdx->get_groupit_member_idx(*i))
						==gr2erase.end())
							gr2erase.push_back(memIdx->get_groupit_member_idx(*i));
					}
				}
			}
			member_iter1++;
		}
		group_iter1++;
	}

	for(std::list<GROUPLIST::iterator>::iterator i=gr2erase.begin();
	i!=gr2erase.end();i++)
	{
		for(std::list<unsigned>::iterator m=(*i)->begin();
		m!=(*i)->end();m++)
		{
			//std::cout<<"remove member:"<<*m<<std::endl;
			memIdx->erase_member_idx(*m);
		}
		group_list.erase(*i);
	}
	if(verbose>0)
		std::cout<<" done"<<std::endl;
}

/*-------------------------------------------------
  - Function to select groups according to their size
  - Argument: limiting size
  --------------------------------------------------*/
void BLRGroup::group_size_filter( int verbose )
{
	if(verbose>0)
		std::cout<<"Filtering groups..."<<std::flush;
  GROUPLIST::iterator group_iter=group_list.begin();
  while(group_iter!=group_list.end())
    {
      if (group_iter->size()<grouper_parameter->getSizefilter())
	{
	  for(std::list<unsigned>::iterator m=(group_iter)->begin();
	      m!=(group_iter)->end();m++)
	    memIdx->erase_member_idx(*m);
	  group_iter=group_list.erase(group_iter);
	}
      else
	group_iter++;
    }
  if(verbose>0)
	  std::cout<<" done"<<std::endl;
}

unsigned BLRGroup::getNbMembers(void)
{
	unsigned nbMembers = 0;
	for( GROUPLIST::iterator group_iter=group_list.begin();
	group_iter!=group_list.end(); group_iter++ )
		nbMembers += (*group_iter).size();
	return nbMembers;
}

