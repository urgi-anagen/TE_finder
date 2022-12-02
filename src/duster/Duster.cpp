#include "Duster.h"


//-------------------------------------------------------------------------
void Duster::search(const BioSeq& sequence, unsigned numseq, unsigned start, unsigned end, bool repeat, std::vector< std::pair<unsigned,unsigned> >& fmerged)
{
	clock_t clock_begin, clock_end;
	clock_begin = clock();
	std::cout<<"hashing query sequence..."<<std::flush;
	Diag_map diag_map(nbseqS+1);
    matchKmersHole(sequence, numseq, start, end, repeat, diag_map);
	std::cout<<"ok"<<std::endl;
	std::cout<<diag_map.size()<<" hits found";
	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;

	clock_begin = clock();
	std::cout<<"search fragments..."<<std::flush;
	std::vector< std::pair<unsigned,unsigned> > frag;
	diagSearchDist(diag_map,(kmer_size+1)*wdist,kmer_size,frag);
	diag_map.clear();
	std::cout<<"ok"<<std::endl;
	std::cout<<frag.size()<<" fragments found";
	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;

	clock_begin = clock();
	std::cout<<"merge fragments..."<<std::flush;
	fragMerge(frag,fdist,fmerged);
	std::cout<<"ok"<<std::endl;
	std::cout<<fmerged.size()<<" ranges found";
	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
}
//-------------------------------------------------------------------------
// merge found fragments
void Duster::fragMerge(std::vector< std::pair<unsigned,unsigned> >& frag,
                           unsigned connect_dist,
                           std::vector< std::pair<unsigned,unsigned> >& fmerged)
{

    sort(frag.begin(),frag.end());
    unsigned start=0;
    unsigned end=0;
    unsigned size=frag.size();

    if(size>=2)
      {
		//search for consecutive kmer matches on direct strand
		unsigned prev_start=frag[0].first;
		unsigned prev_end=frag[0].second;
		unsigned curr_start;
		unsigned curr_end;
		for( unsigned i=1; i<size; ++i)
		{
			  curr_start=frag[i].first;
			  curr_end=frag[i].second;
			  if(prev_end+connect_dist>=curr_start) // consecutive or overlapping matches
				{
				  if(start!=0) //extend the current match range
					{
					  end=std::max(end,curr_end);
					  prev_start=0; // means that previous range must be forgotten
					}
				  else //create match range
					{
					  start=prev_start;
					  prev_start=0; // means that previous range must be forgotten
					  end=std::max(prev_end,curr_end);
					}
				}
			  else
			  {
				  if(start!=0 && end-start>=min_size) //match range created
				  {
					fmerged.emplace_back(start,end);
					start=0;
				  }
				  else
				  {
					  if(prev_end-prev_start>=min_size && prev_start!=0)
						  fmerged.emplace_back(prev_start,prev_end);
					  start=0;
				  }
			  }
			prev_start=curr_start;
            if (curr_end > end) { prev_end = curr_end; }
            else if (start != 0) prev_end = end;
            else prev_end = curr_end;
		}
		if(start!=0) //end of sequence reached, must save the match range if one is in extension
		{
			if(end-start>=min_size)
				fmerged.emplace_back(start,end);
		}
		else
		{
			if(prev_end-prev_start>=min_size)
				fmerged.emplace_back(prev_start,prev_end);
		}
      }
}
//-------------------------------------------------------------------------
void Duster::writeBED(const SDGString& qname, const std::vector< std::pair<unsigned,unsigned> >& frag, std::ostream& out)
{
  unsigned size=frag.size();
    for (unsigned i = 0; i < size; ++i) {
        if (frag[i].first < frag[i].second) {
            out << qname << "\t" // chromosome
                << frag[i].first  //chrom start
                << "\t" << frag[i].second  //chrom end
                << "\t" << "duster" // name
                << "\t" << frag[i].second - frag[i].first //score
                << "\t" << "+" //strand
                << std::endl;
        } else {
            out << qname << "\t" // chromosome
                << frag[i].second  //chrom start
                << "\t" << frag[i].first  //chrom end
                << "\t" << "duster" // name
                << "\t" << frag[i].first - frag[i].second //score
                << "\t" << "-" //strand
                << std::endl;
        }
    }
}
//-------------------------------------------------------------------------
unsigned Duster::compute_coverage(const std::vector< std::pair<unsigned,unsigned> >& frag)
{
  unsigned size=frag.size();
  unsigned coverage=0;
  for(unsigned i=0; i<size; ++i)
    {
        if (frag[i].first < frag[i].second)
            coverage += frag[i].second - frag[i].first;
        else
            coverage += frag[i].first - frag[i].second;
    }
  return coverage;
}
//-------------------------------------------------------------------------
void Duster::get_sequences(const std::vector< std::pair<unsigned,unsigned> >& frag, BioSeq& seq, FastaOstream& out)
{

    unsigned size = frag.size();
    for (unsigned i = 0; i < size; ++i) {
        std::ostringstream name;
        name << seq.header << ":" << frag[i].first << ".." << frag[i].second;
        if (frag[i].first < frag[i].second) {
            BioSeq subseq = seq.subseq(frag[i].first, frag[i].second - frag[i].first + 1);
            subseq.header =  name.str();
            out << subseq;
        } else {
            BioSeq subseq = seq.subseq(frag[i].second, frag[i].first - frag[i].second + 1);
            BioSeq subseq_comp = subseq.complement();
            subseq_comp.header=name.str();
            out << subseq_comp;
        }
    }
}





