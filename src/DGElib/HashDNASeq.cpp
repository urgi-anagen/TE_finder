#include <cstring>
#include <set>
#include <memory>
#include <limits>
#include <SDGMemBioSeq.h>
#include "HashDNASeq.h"

//-------------------------------------------------------------------------
void HashDNASeq::load(const SDGString& filenameS, unsigned kmerSize, unsigned mask_hole_period, unsigned mask_hole_length, unsigned kmer_window,
                      unsigned bkmerSize, unsigned mkmerSize, double count_cutoff, double diversity_cutoff,
                      unsigned min_count, bool & valid_idx_file, bool first_iter, bool filter_ssr)
{

  nbseqS=0;

  clock_t begin, end;
  begin = clock();


    if (hash_algorithm == 1 || hash_algorithm == 3)
        std::cout<<"Kmer mask :"<<hseq.getMask()<<" -> Effective Kmer size (without kmer mask size):"<<hseq.getEffectiveKmerSize()<<std::endl;
    else
        std::cout<<"Effective Kmer size (without kmer mask size):"<<kmerSize<<std::endl;

    if(valid_idx_file)
	  valid_idx_file=read_idx(filenameS, count_cutoff, diversity_cutoff, min_count, mask_hole_period, mask_hole_length);
  if(!valid_idx_file)
    {
	  //Count kmers and filters
      std::vector<unsigned> kmer_count;
        if (hash_algorithm == 1 || hash_algorithm == 3)
            kmer_count.resize((unsigned) pow(4, hseq.getEffectiveKmerSize()), 0);
        else
            kmer_count.resize((unsigned) pow(4,kmerSize),0);
	  unsigned nb_kmer=0;
	  std::list< Info_kmer > list_infokmer;
	  Info_kmer kmer_threshold;
	  kmer_analysis(filenameS, kmerSize, mask_hole_period, mask_hole_length, kmer_window, bkmerSize, mkmerSize, count_cutoff,
                    diversity_cutoff, kmer_count, nb_kmer, list_infokmer, kmer_threshold);
	  if(filter_ssr) kmer_ssr_filter(kmerSize, kmer_count);
      kmer_filter(list_infokmer, kmer_threshold, min_count, kmer_count, first_iter);
      end = clock();
      std::cout<<" --> Time spent: "<<(double)(end-begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;


      //Prepare hash table of kmer positions
      begin = clock();
      std::cout<<"Prepare hash table"<<std::endl;
      kmer_pos.clear();
      kmer_pos.resize(nb_kmer);
      unsigned k=0;
      auto last_it=kmer_pos.begin();
      for(unsigned int & i : kmer_count)
		{
		  hash2wpos[k]=last_it;
		  last_it=last_it+i;
		  k++;
		}
      hash2wpos[k]=last_it;

      //Load hash table with kmer positions
      std::cout<<"Load kmer positions in hash table"<<std::endl;
      hash_ptr=hash2wpos;
      FastaIstream inS2(filenameS);
      unsigned count_seq=0;
      while(inS2)
		{
		  BioSeq sS;
		  if(inS2)
          {
              inS2>>sS;
              count_seq++;
              if (count_seq==std::numeric_limits<unsigned>::max())
                  throw Unsigned_Out_Of_Range("Number of subject sequence is out of range: ",count_seq);
          }
            if(hash_algorithm==1)
                hashSubjectSeqPosWHole(sS, kmer_count);
            else if(hash_algorithm==0)
                hashSubjectSeqPos(sS, count_seq, kmerSize, kmer_count);
            else if(hash_algorithm==2)
                hashSubjectSeqPosMinimizer(sS, count_seq, kmerSize, kmer_count);
            else if(hash_algorithm==3)
                hashSubjectSeqPosWHoleMinimizer(sS, count_seq, kmerSize, kmer_count);
		}
      inS2.close();
      hash_ptr.clear();
//      save_idx(filenameS,count_cutoff,diversity_cutoff,min_count,mask_hole_period,mask_hole_length,kmer_count);
//      valid_idx_file=true;
    }

  std::cout<<"end preparation of subjects";
  end = clock();
  std::cout<<" --> Time spent: "<<(double)(end-begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
}
//-------------------------------------------------------------------------
void HashDNASeq::kmer_analysis(const SDGString& filenameS,
                               unsigned kmerSize, unsigned mask_hole_period, unsigned mask_hole_length,
                               unsigned kmer_window,
                               unsigned bkmerSize, unsigned mkmerSize,
                               double count_cutoff, double diversity_cutoff,
                               std::vector<unsigned>& kmer_count, unsigned& nb_kmer,
                               std::list< Info_kmer >& list_infokmer, Info_kmer& kmer_threshold)
{

	  std::vector<unsigned> background_count((unsigned)pow(4, bkmerSize), 0);
	  std::vector<unsigned> model_count((unsigned)pow(4, mkmerSize), 0);
	  std::vector<unsigned> nuc_count(4,0);

  unsigned nb_bkmer=0;
  unsigned nb_mkmer=0;
  unsigned nb_nuc=0;
  nb_kmer=0;

  kmer_counts(filenameS, kmerSize, bkmerSize, mkmerSize,
              kmer_count, nb_kmer,
              background_count, nb_bkmer,
              model_count, nb_mkmer,
              nuc_count, nb_nuc);

  if(nb_kmer==0)
      throw "Error: No more kmer to use!";


  kmer_count[0]=0; //remove_self_hits kmers AAAAAA... NNNNNN.... XXXXX....
  kmer_prob(kmerSize, bkmerSize, mkmerSize, mask_hole_period, mask_hole_length, kmer_window,
            kmer_count, nb_kmer,
            background_count, nb_bkmer,
            model_count, nb_mkmer,
            nuc_count, nb_nuc,
            list_infokmer);

    kmer_occurrence_percentiles(list_infokmer, count_cutoff, kmer_threshold);

  double cutoff_entropy=1.0;
  kmer_entropy_percentiles(list_infokmer, cutoff_entropy, kmer_threshold);

  kmer_diversity_percentiles(list_infokmer, diversity_cutoff, kmer_threshold);
  kmer_goodkmer_percentiles(list_infokmer);

}
//-------------------------------------------------------------------------
void HashDNASeq::kmer_counts(const SDGString& filenameS, unsigned kmerSize, unsigned bkmerSize, unsigned mkmerSize,
                             std::vector<unsigned>& wcount, unsigned& nb_kmer,
                             std::vector<unsigned>& bcount, unsigned& nb_bkmer,
                             std::vector<unsigned>& mcount, unsigned& nb_mkmer,
                             std::vector<unsigned>& ncount, unsigned& nb_nuc)
{
  std::cout << "Counting kmer of size " << kmerSize << " ... " << std::flush;

  try {
      FastaIstream inS(filenameS);
      unsigned count_seq=0;
      BioSeq sS;
      while (inS) {
          if (inS) inS >> sS;
          std::cout << "\n" << ++count_seq << "->" << sS.header << std::flush;
          if (hash_algorithm == 1 || hash_algorithm == 3)
              nb_kmer += hashSeqCountWHole(sS, kmerSize, wcount);
          else
              nb_kmer += hashSeqCount(sS, kmerSize, wcount);
          nb_bkmer += hashSeqCount(sS, bkmerSize, bcount);
          nb_mkmer += hashSeqCount(sS, mkmerSize, mcount);
          nb_nuc += hashSeqCount(sS,1, ncount);
      }
      inS.close();
      std::cout<<"\nNumber of kmers="<<nb_kmer<<std::endl;
  }
  catch (const std::ifstream::failure& e) {
      std::cout<<"Error opening file "<<filenameS<<std::endl;
  }
}
//-------------------------------------------------------------------------
//
void HashDNASeq::kmer_occurrence_percentiles(const std::list< Info_kmer >& list_infokmer, double cutoff,
                                             Info_kmer& kmer_threshold)
{
	std::cout<<"\nKmer occurrence percentiles:"<<std::endl;

    std::vector<unsigned> distr;
	for(const auto & i : list_infokmer)
    	distr.emplace_back(i.count);
    sort(distr.begin(),distr.end());

    std::vector<unsigned> clean_distr;
    copy(lower_bound(distr.begin(),distr.end(),(unsigned)1), distr.end(),back_inserter(clean_distr));
    std::cout<<" min="<<clean_distr.front()<<std::endl;
    std::cout<<" 1st quartile="
	       << clean_distr[(unsigned)floor(0.25*(double)clean_distr.size())]<<std::endl;
    std::cout<<" median="
	       << clean_distr[(unsigned)floor(0.5*(double)clean_distr.size())]<<std::endl;
    std::cout<<" 3st quartile="
	       << clean_distr[(unsigned)floor(0.75*(double)clean_distr.size())]<<std::endl;
    std::cout<<" 95 percentile="
	       << clean_distr[(unsigned)floor(0.95*(double)clean_distr.size())]<<std::endl;
    std::cout<<" 99 percentile="
	       << clean_distr[(unsigned)floor(0.99*(double)clean_distr.size())]<<std::endl;
    std::cout<<" max="<< clean_distr.back()<<std::endl;
    long threshold_index=(long)(cutoff*(double)distr.size())-1;
    if(threshold_index<0) threshold_index=0;
    kmer_threshold.count=distr[threshold_index];

	unsigned count=0;
	for(unsigned int & it : distr)
	{
		if(it>kmer_threshold.count)
			count++;
	}
    std::cout<<"=>cut-off="<<cutoff<<". Kmers occurring more than "<<kmer_threshold.count
	       <<" will be removed! Will remove_self_hits "<<count<<" kmers."<<std::endl;
}
//-------------------------------------------------------------------------
//
void HashDNASeq::kmer_entropy_percentiles(const std::list< Info_kmer >& list_infokmer, double cutoff_entropy,
		Info_kmer& kmer_threshold)
{

	std::cout<<"\nEntropy percentiles:"<<std::endl;

	std::vector<double> distr;
	for(const auto & i : list_infokmer)
    	distr.emplace_back(i.entropy);
    sort(distr.begin(),distr.end());

    std::cout<<" min="<<distr.front()<<std::endl;
    std::cout<<" 1st quartile="
	       << distr[(unsigned)floor(0.25*(double)distr.size())]<<std::endl;
    std::cout<<" median="
	       << distr[(unsigned)floor(0.5*(double)distr.size())]<<std::endl;
    std::cout<<" 3st quartile="
	       << distr[(unsigned)floor(0.75*(double)distr.size())]<<std::endl;
    std::cout<<" 95 percentile="
	       << distr[(unsigned)floor(0.95*(double)distr.size())]<<std::endl;
    std::cout<<" 99 percentile="
	       << distr[(unsigned)floor(0.99*(double)distr.size())]<<std::endl;
    std::cout<<" max="<< distr.back()<<std::endl;
    long threshold_index=(long)(cutoff_entropy*(double)distr.size())-1;
    if(threshold_index<0) threshold_index=0;
    kmer_threshold.entropy=distr[threshold_index];

	unsigned count=0;
	for(double & it : distr)
	{
		if(it>kmer_threshold.entropy)
			count++;
	}
    std::cout<<"=>cut-off="<<cutoff_entropy<<". Kmers with entropy greater than "<<kmer_threshold.entropy
	       <<" will be removed! Will remove_self_hits "<<count<<" kmers."<<std::endl;
}
//-------------------------------------------------------------------------
//
void HashDNASeq::kmer_diversity_percentiles(const std::list< Info_kmer >& list_infokmer, double diversity_cutoff,
		Info_kmer& kmer_threshold)
{

	std::cout<<"\nDiversity (# of different background kmers / # of background kmer in kmer) percentiles:"<<std::endl;

	std::vector<double> distr;
	for(const auto & i : list_infokmer)
    	distr.emplace_back(i.diversity);
    sort(distr.begin(),distr.end());

    std::cout<<" min="<<distr.front()<<std::endl;
    std::cout<<" 1st quartile="
	       << distr[(unsigned)floor(0.25*(double)distr.size())]<<std::endl;
    std::cout<<" median="
	       << distr[(unsigned)floor(0.5*(double)distr.size())]<<std::endl;
    std::cout<<" 3st quartile="
	       << distr[(unsigned)floor(0.75*(double)distr.size())]<<std::endl;
    std::cout<<" 95 percentile="
	       << distr[(unsigned)floor(0.95*(double)distr.size())]<<std::endl;
    std::cout<<" 99 percentile="
	       << distr[(unsigned)floor(0.99*(double)distr.size())]<<std::endl;
    std::cout<<" max="<< distr.back()<<std::endl;


	kmer_threshold.diversity=diversity_cutoff;
	unsigned count=0;
	for(double & it : distr)
	{
		if(it>=kmer_threshold.diversity)
			break;
		count++;
	}
	std::cout<<"=>diversity_cutoff="<<kmer_threshold.diversity<<". Kmers with diversity less than "<<kmer_threshold.diversity
		   <<" will be removed! Will remove_self_hits "<<count<<" kmers."<<std::endl;

}
//-------------------------------------------------------------------------
//
void HashDNASeq::kmer_goodkmer_percentiles(const std::list< Info_kmer >& list_infokmer)
{

	std::cout<<"\nGood kmer (kmer model log-likelihood) percentiles:"<<std::endl;

	std::vector<double> distr;
	for(const auto & i : list_infokmer)
    	distr.emplace_back(i.good_kmer);
    sort(distr.begin(),distr.end());

    std::cout<<" min="<<distr.front()<<std::endl;
    std::cout<<" 1st quartile="
	       << distr[(unsigned)floor(0.25*(double)distr.size())]<<std::endl;
    std::cout<<" median="
	       << distr[(unsigned)floor(0.5*(double)distr.size())]<<std::endl;
    std::cout<<" 3st quartile="
	       << distr[(unsigned)floor(0.75*(double)distr.size())]<<std::endl;
    std::cout<<" 95 percentile="
	       << distr[(unsigned)floor(0.95*(double)distr.size())]<<std::endl;
    std::cout<<" 99 percentile="
	       << distr[(unsigned)floor(0.99*(double)distr.size())]<<std::endl;
    std::cout<<" max="<< distr.back()<<std::endl;


	unsigned count=0;
	for(double & it : distr)
	{
		if(it>0.0)
			break;
		count++;
	}
	std::cout<<"Kmers with poor model probability will be removed ! Will remove_self_hits "<<count<<" kmers."<<std::endl;

}
//-------------------------------------------------------------------------
//
void HashDNASeq::kmer_prob(unsigned wsize, unsigned bwsize, unsigned mwsize, unsigned mask_hole_period, unsigned mask_hole_length, unsigned kmer_window,
                           const std::vector<unsigned>& wcount, unsigned nb_kmer,
                           const std::vector<unsigned>& bcount, unsigned nb_bkmer,
                           const std::vector<unsigned>& mcount, unsigned nb_mkmer,
                           const std::vector<unsigned>& ncount, unsigned nb_nuc,
                           std::list< Info_kmer >& list_infokmer) const
{
	std::cout<<"\nCompute kmer stats ... "<<std::endl;
    HashDNASeq h(wsize, mask_hole_period, mask_hole_length, kmer_window, hash_algorithm, bwsize, 0);
	std::cout << "kmer size=" << wsize
              << "\tmask_hole_period period=" << mask_hole_period
            <<"\tmask_hole_period hole size="<<mask_hole_length<<std::flush;
    if (hash_algorithm == 1) { std::cout << "\teffective kmer size=" << h.hseq.getEffectiveKmerSize() << std::endl; }
    else { std::cout << "\teffective kmer size=" << h.mseq.getEffectiveKmerSize() << std::endl; }

	std::cout<<"\nSequence composition:\n";

	//Compute nucleotide probabilities
	auto max_key_bkmer=(unsigned)pow(4,bwsize);
	auto max_key_mkmer=(unsigned)pow(4,mwsize);

    std::vector<double> bprob(max_key_bkmer,0),prob_ind(max_key_bkmer,0),bprob_cond(max_key_bkmer,0),nprob(4,0),
    		mprob(max_key_mkmer,0),mprob_cond(max_key_mkmer,0);

    for(unsigned k=0; k<4;k++)
     {
     	nprob[k]=(double)ncount[k]/nb_nuc;
     	std::cout<<"freq("<<h.nhseq.reverse_hash(k)<<")="<<nprob[k]<<std::endl;
     }
    std::cout<<"\nAT%="<<(nprob[h.nhseq.hash("A")]+nprob[h.nhseq.hash("T")])*100<<std::endl;
    std::cout<<"GC%="<<(nprob[h.nhseq.hash("G")]+nprob[h.nhseq.hash("C")])*100<<std::endl;

    //Compute kmer background probabilities
    for(unsigned k=0; k<max_key_bkmer;k++)
    {
    	bprob[k]=(double)(bcount[k]+0.0001)/nb_bkmer; //use pseudo count!!
    }

    //Compute kmer model probabilities
    for(unsigned k=0; k<max_key_mkmer;k++)
    {
    	mprob[k]=(double)(mcount[k]+0.0001)/nb_mkmer; //use pseudo count!!
    }

    //Compute kmer background independance and conditional probabilities
    std::list< std::pair <double,unsigned> > list_prob;
    for(unsigned k=0; k<max_key_bkmer;k++)
    {
    	unsigned order=bwsize-1;
    	std::string kmer=h.bhseq.reverse_hash(k);
    	double sum_prob=0;

    	std::string kmer2=kmer.substr(0,order)+"A";
    	const char* s=kmer2.c_str();
    	sum_prob+=bprob[h.bhseq(s)];

    	kmer2=kmer.substr(0,order)+"C";
    	s=kmer2.c_str();
    	sum_prob+=bprob[h.bhseq(s)];

    	kmer2=kmer.substr(0,order)+"G";
    	s=kmer2.c_str();
    	sum_prob+=bprob[h.bhseq(s)];

    	kmer2=kmer.substr(0,order)+"T";
    	s=kmer2.c_str();
    	sum_prob+=bprob[h.bhseq(s)];

    	bprob_cond[k]=bprob[k]/sum_prob;
    	list_prob.emplace_back(bprob[k],k);

    	prob_ind[k]=1;
    	for(unsigned n=0; n<bwsize; n++)
    	{
    		kmer2=kmer.substr(n,1);
    		const char* n1=kmer2.c_str();
    		prob_ind[k]*=nprob[h.nhseq(n1)];
    	}
    }

    list_prob.sort();
    unsigned prob2display=0,maxProb2display=20;
    std::cout<<"\nBackground frequency and probabilities\n";
    for(auto i=list_prob.rbegin()  ; i!=list_prob.rend() && prob2display < maxProb2display; i++)
         {
    		unsigned k=i->second;
    		std::cout<<"freq("<<h.bhseq.reverse_hash(k)<<")="<<bprob[k]
                  <<"\tprob_ind("<<h.bhseq.reverse_hash(k)<<")="<<prob_ind[k]
                  <<"\tprob_cond("<<h.bhseq.reverse_hash(k)<<")="<<bprob_cond[k]
       	          <<std::endl;
    		prob2display++;
         }
    if(prob2display == maxProb2display)
    	std::cout<<"..... best "<<maxProb2display;
    std::cout<<std::endl;

    std::cout<<"Model frequency and probabilities ..."<<std::flush;
    // Compute model conditional probability
    for(unsigned k=0; k<max_key_mkmer;k++)
     {
     	unsigned order=mwsize-1;
     	std::string kmer=h.mhseq.reverse_hash(k);
     	double sum_prob=0;

     	std::string kmer2=kmer.substr(0,order)+"A";
     	const char* s=kmer2.c_str();
     	sum_prob+=mprob[h.mhseq(s)];

     	kmer2=kmer.substr(0,order)+"C";
     	s=kmer2.c_str();
     	sum_prob+=mprob[h.mhseq(s)];

     	kmer2=kmer.substr(0,order)+"G";
     	s=kmer2.c_str();
     	sum_prob+=mprob[h.mhseq(s)];

     	kmer2=kmer.substr(0,order)+"T";
     	s=kmer2.c_str();
     	sum_prob+=mprob[h.mhseq(s)];

     	mprob_cond[k]=mprob[k]/sum_prob;
     }
    std::cout<<"done!"<<std::endl;

    //Compute background and model probability of kmers
    std::cout<<"Observed Kmer frequency and probabilities ..."<<std::flush;
    long double khi2=0;
    for(unsigned k=0;k<max_key;k++)
    {
    	if(wcount[k]==0) continue;

    	std::string kmer=h.hseq.reverse_hash(k);
    	const char* s=kmer.c_str();

    	std::list<unsigned> list_skmer;
    	list_skmer.push_back(h.bhseq(s));
    	double pi=bprob[h.bhseq(s)];
    	double entropy=pi*log2(pi);
    	double prob_bkmer=bprob[h.bhseq(s++)];
    	for(unsigned i=1; i<=wsize-bwsize && *s!='\0'; i++)
    	{
    		list_skmer.push_back(h.bhseq(s));
    		pi=bprob[h.bhseq(s)];
    		entropy+=pi*log2(pi);
    		prob_bkmer*=bprob_cond[h.bhseq(s++)];
    	}
    	double ecount=prob_bkmer*nb_kmer;
    	list_skmer.sort();
    	list_skmer.unique();
       	double diversity=(double)list_skmer.size()/(wsize-bwsize+1);

       	s=kmer.c_str();

    	double prob_mkmer=mprob[h.mhseq(s++)];
    	for(unsigned i=1; i<=wsize-mwsize && *s!='\0'; i++)
    	{
    		prob_mkmer*=mprob_cond[h.mhseq(s++)];
    	}
    	double belong2model=log10(prob_mkmer)-log10(prob_bkmer);

    	list_infokmer.emplace_back(k,wcount[k],ecount,entropy*-1,diversity,belong2model);
    	khi2+=pow((wcount[k]-ecount),2)/(ecount);
    	if(wcount[k]==0) continue;
    	//count_kmer++;
    }
    std::cout<<"done!"<<std::endl;
    std::cout<<"\nKhi2 between kmers count and expected:"<<khi2<<"\n";

    list_infokmer.sort(Compare_Info_kmer());
    std::cout<<"\nDisplay most abundant kmers (count/expected in parenthesis):\n";
    long double khi2_most_abundant=0;
    unsigned count2display=0,maxCount2display=20;
    for(auto i=list_infokmer.rbegin(); i!=list_infokmer.rend() && count2display < maxCount2display; i++)
      {
       	unsigned c=i->count;
        double bc=i->expected_count;
        khi2_most_abundant+=pow((c-bc),2)/(bc);
        std::cout<<h.hseq.reverse_hash(i->hash_key)<<" "<<list_infokmer.size()-count2display
        		<<"="<<c<<" ("<<c/bc<<") entropy="<<i->entropy<<" diversity "<<bwsize<<"-mer="<<i->diversity
        		<<" good kmer ratio="<<i->good_kmer<<std::endl;
        count2display++;
      }
    std::cout<<"\nKhi2 between most abundant kmers count and expected:"<<khi2_most_abundant<<"\n";

    std::cout<<"\nDisplay less abundant kmers (count/expected in parenthesis):\n";
    long double khi2_less_abundant=0;
    count2display=0,maxCount2display=20;
    for(auto i=list_infokmer.begin(); i!=list_infokmer.end() && count2display < maxCount2display; i++)
      {
    	unsigned c=i->count;
    	double bc=i->expected_count;
    	khi2_less_abundant+=pow((c-bc),2)/(bc);
    	std::cout<<h.hseq.reverse_hash(i->hash_key)<<" "<<count2display+1
    			<<"="<<c<<" ("<<c/bc<<") entropy="<<i->entropy<<" diversity "<<bwsize<<"-mer="<<i->diversity
    			<<" good kmer ratio="<<i->good_kmer<<std::endl;
    	count2display++;
      }
    std::cout<<"\nKhi2 between less abundant kmers count and expected:"<<khi2_less_abundant<<"\n";
}
//-------------------------------------------------------------------------
//
void HashDNASeq::kmer_filter(const std::list< Info_kmer >& list_infokmer, const Info_kmer& kmer_threshold, unsigned min_count,
			std::vector<unsigned>& wcount, bool first_iter) const
{
    std::cout<<"\nFilter kmers ... "<<std::flush;
    unsigned nb_kmer=0;
    unsigned nb_removed=1; // 1 for AAAAAA...
    for(const auto & i : list_infokmer)
    {
    	unsigned k=i.hash_key;
    	if(wcount[k]==0) continue;
    	nb_kmer++;
    	if(i.count<min_count*(i.expected_count))
    	{
    		wcount[k]=0;
    		nb_removed++;
    	}
    	if(i.count>kmer_threshold.count)
    	{
    		wcount[k]=0;
    		nb_removed++;
    	}
    	if(i.entropy>kmer_threshold.entropy)
    	{
    		wcount[k]=0;
    		nb_removed++;
    	}
    	if(i.diversity<kmer_threshold.diversity)
    	{
    		wcount[k]=0;
    		nb_removed++;
    	}
    	if(!first_iter && i.good_kmer<0.0)
    	{
    		wcount[k]=0;
    		nb_removed++;
    	}
    }
    std::cout<<"removed kmers:"<<nb_removed<<" over "<<nb_kmer<<" existing ("<<(double)nb_removed/nb_kmer<<") and "<<max_key+1<<" possible."<<std::endl;
}
//-------------------------------------------------------------------------
//
void HashDNASeq::kmer_ssr_filter(unsigned wsize, std::vector<unsigned>& wcount) {
    std::cout << "\nFilter SSR kmers ... " << std::endl;
    unsigned nb_removed = 0;

    std::list<std::string> list_ssr;
    HashDNASeq h2(2,0,2, 0);
    for(unsigned i=0; i<16; i++)
    {list_ssr.push_back(h2.hseq.reverse_hash(i));}
    HashDNASeq h3(3,0,3, 0);
    for(unsigned i=0; i<64; i++)
    {list_ssr.push_back(h3.hseq.reverse_hash(i));}

    for(auto & it : list_ssr)
    {
        unsigned l = 0;
        std::ostringstream oseq;
        while (l < 18) {
            oseq << it;
            l += it.size();
        }
        SDGBioSeq seq_ssr = newSDGMemBioSeq(oseq.str());
        unsigned len = seq_ssr.length();
        // std::cout << "SSR:" << seq_ssr.toString() <<" len="<<len<<std::endl;
        unsigned last_pos = len - wsize;
        char *s = new char[len + 1];
        char *ptr=s;
        std::strcpy(s, seq_ssr.toString().c_str());
        for (unsigned i = 0; i < last_pos; i++) {
            unsigned hval=hseq(s);
            if(wcount[hval]!=0){
                std::cout <<"kmer:"<< hseq.reverse_hash(hval) <<" count:"<< wcount[hval]<<" removed !" << std::endl;
                wcount[hval]=0;
                nb_removed++;
            }
            s++;
        }
        delete [] ptr;
    }

    std::cout << "removed kmers:" << nb_removed << std::endl;
}
//-------------------------------------------------------------------------
// Read a hash index
bool HashDNASeq::read_idx(const SDGString& filename, double count_cutoff , double diversity_cutoff, unsigned min_count, unsigned kmask, unsigned mask_hole_length)
{
  std::ifstream in(filename+".kidx");
  if(!in) return false;
  std::vector<unsigned> kmer_count;
  unsigned val;
  in>>val;
  if(val!=kmer_size) return false;
  in>>val;
  if(val!=kmask) return false;
  in>>val;
  if(val!=mask_hole_length) return false;
  in>>val;
  if(val!=min_count) return false;
  double fval;
  in>>fval;
  if(fval!=count_cutoff) return false;
  in>>fval;
  if(fval!=diversity_cutoff) return false;

  in>>val;
  if(val!=nbseqS) return false;
  nbseqS=val;

  std::cout<<"Read existing idx file:"<<filename+".kidx"<<std::endl;
  in>>val;
  for(unsigned i=0;i<val;i++)
    {
      unsigned p, n;
      in>>p>>n;
      KmerSpos wpos(p,n);
      kmer_pos.push_back(wpos);
    }
  kmer_count.clear();
  in>>val;
  for(unsigned i=0;i<val;i++)
    {
      unsigned n;
      in>>n;
      kmer_count.push_back(n);
    }
  in.close();

  std::cout<<"Prepare hash table pointers"<<std::endl;
  unsigned k=0;
  auto last_it=kmer_pos.begin();
  for(unsigned int & i : kmer_count)
    {
      hash2wpos[k]=last_it;
      last_it=last_it+i;
      k++;
    }
  hash2wpos[k]=last_it;

  return true;
}
//-------------------------------------------------------------------------
// Save a hash index
void HashDNASeq::save_idx(const SDGString& filename, double count_cutoff, double diversity_cutoff, unsigned min_count, unsigned kmask, unsigned mask_hole_length, const std::vector<unsigned>& wcount)
{
  std::ofstream fout(filename+".kidx");
  std::cout<<"Save idx"<<std::endl;
  fout<<kmer_size<<std::endl;
  fout<<kmask<<std::endl;
  fout<<mask_hole_length<<std::endl;
  fout<<min_count<<std::endl;
  fout<<count_cutoff<<std::endl;
  fout<<diversity_cutoff<<std::endl;
  fout<<nbseqS<<std::endl;
  fout<<kmer_pos.size()<<std::endl;
  for(auto & kmer_po : kmer_pos)
    {
      fout<<kmer_po.pos<<" "<<kmer_po.numSeq<<"\t";
    }
  fout<<wcount.size()<<std::endl;
  for(unsigned int i : wcount)
    {
      fout<<i<<"\t";
    }
}
//-------------------------------------------------------------------------
// Count kmers
unsigned HashDNASeq::hashSeqCountWHole(const BioSeq& seq, unsigned wsize, std::vector<unsigned>& wcount)
{
  unsigned len=seq.length();
  unsigned nb_kmer=0;
  if(len<=wsize) return 0;
  unsigned last_pos=len-wsize;
  const char *s=seq.c_str();
  for(unsigned i=0;i<=last_pos;i++)
    {
	  nb_kmer++;
      wcount[hseq(s)]++;
      s++;
    }
  return nb_kmer;
}
//-------------------------------------------------------------------------
// Count kmers
unsigned HashDNASeq::hashSeqCount(const BioSeq& seq, unsigned wsize, std::vector<unsigned>& wcount)
{
    unsigned len=seq.length();
    unsigned nb_kmer=0;
    if(len<=wsize) return 0;
    const char *s=seq.c_str();
    unsigned h=0;
    unsigned val_nuc=0;
    unsigned i=0;
    unsigned mask=3;
    mask<<=(wsize-1)*2;
    while (*s != '\0') {
        switch (*s) {
            case 'C': {
                val_nuc = 3;
                break;
            }
            case 'G': {
                val_nuc = 2;
                break;
            }
            case 'T': {
                val_nuc = 1;
                break;
            }
            case 'c': {
                val_nuc = 3;
                break;
            }
            case 'g': {
                val_nuc = 2;
                break;
            }
            case 't': {
                val_nuc = 1;
                break;
            }
            default: {
                val_nuc = 0;
                break;
            }
        }

        //Add new nucleotide value
        h <<= 2;
        h |= val_nuc;
        i++;
        if(i>=wsize) {
            wcount[h]++;
            //Suppress old boundary value
            unsigned old_val = h & mask;
            h  ^= old_val;
            nb_kmer++;
        }
        s++;
    }
    return nb_kmer;
}

//-------------------------------------------------------------------------
void HashDNASeq::hashSubjectSeqPosWHole(const BioSeq& seq, const std::vector<unsigned>& wcount)
{
    nbseqS++;
    unsigned len = seq.length();
    if (len <= kmer_size) return;
    unsigned last_pos = len - kmer_size;
    unsigned key;
    const char *s=seq.c_str();
    for (unsigned i = 0; i <= last_pos; i++) {
        key = hseq(s);

        if (wcount[key] != 0) {
            *(hash_ptr[key]) = KmerSpos(i, nbseqS);
            hash_ptr[key]++;
        }
        s++;
    }
}
//-------------------------------------------------------------------------
void HashDNASeq::hashSubjectSeqPos(const BioSeq &seq, unsigned num_seq, unsigned wsize, const std::vector<unsigned> &wcount)
{
    unsigned nb_kmer=0;
    const char *s=seq.c_str();
    unsigned h=0;
    unsigned val_nuc=0;
    unsigned i=0;
    unsigned mask=3;
    mask<<=(wsize-1)*2;
    while (*s != '\0') {
        switch (*s) {
            case 'C': {
                val_nuc = 3;
                break;
            }
            case 'G': {
                val_nuc = 2;
                break;
            }
            case 'T': {
                val_nuc = 1;
                break;
            }
            case 'c': {
                val_nuc = 3;
                break;
            }
            case 'g': {
                val_nuc = 2;
                break;
            }
            case 't': {
                val_nuc = 1;
                break;
            }
            default: {
                val_nuc = 0;
                break;
            }
        }

        //Add new nucleotide value
        h <<= 2;
        h |= val_nuc;
        i++;
        if(i>=wsize) {
            if (wcount[h] != 0) {
                *(hash_ptr[h]) = KmerSpos(i - wsize, num_seq);
                hash_ptr[h]++;
            }
            //Suppress old boundary value
            unsigned old_val = h & mask;
            h  ^= old_val;
            nb_kmer++;
        }
        s++;
    }
}
//-------------------------------------------------------------------------
void HashDNASeq::hashSubjectSeqPosMinimizer(const BioSeq &seq, unsigned num_seq, unsigned wsize, const std::vector<unsigned> &wcount)
{
    std::list<std::pair<unsigned,unsigned>> kmer_pos_list;
    std::set<std::pair<unsigned,unsigned>> minimized_kmer_pos_list;
    hashSeqPos(seq,wsize,kmer_pos_list);
    minimize(window_size,kmer_pos_list,minimized_kmer_pos_list);

    for(auto &iter_pos : minimized_kmer_pos_list) {
        unsigned pos = iter_pos.second;
        unsigned key_d = iter_pos.first;
        if (wcount[key_d] != 0) {
            *(hash_ptr[key_d]) = KmerSpos(pos, num_seq);
            hash_ptr[key_d]++;
        }

    }
}
//-------------------------------------------------------------------------
void HashDNASeq::hashSubjectSeqPosWHoleMinimizer(const BioSeq &seq, unsigned num_seq, unsigned wsize, const std::vector<unsigned> &wcount)
{
    std::list<std::pair<unsigned,unsigned>> kmer_pos_list;
    std::set<std::pair<unsigned,unsigned>> minimized_kmer_pos_list;
    hashSeqPosWHole(seq,wsize,kmer_pos_list);
    minimize(window_size,kmer_pos_list,minimized_kmer_pos_list);

    for(auto &iter_pos : minimized_kmer_pos_list) {
        unsigned pos = iter_pos.second;
        unsigned key_d = iter_pos.first;
        if (wcount[key_d] != 0) {
            *(hash_ptr[key_d]) = KmerSpos(pos, num_seq);
            hash_ptr[key_d]++;
        }

    }
}
//-------------------------------------------------------------------------
void HashDNASeq::hashSeqPos(const BioSeq &seq, unsigned wsize, std::list<std::pair<unsigned, unsigned> > &kmer_pos_list) {
    unsigned nb_kmer = 0;
    const char *s = seq.c_str();
    unsigned h = 0;
    unsigned val_nuc = 0;
    unsigned i = 0;
    unsigned mask = 3;
    mask <<= (wsize - 1) * 2;
    while (*s != '\0') {
        switch (*s) {
            case 'C': {
                val_nuc = 3;
                break;
            }
            case 'G': {
                val_nuc = 2;
                break;
            }
            case 'T': {
                val_nuc = 1;
                break;
            }
            case 'c': {
                val_nuc = 3;
                break;
            }
            case 'g': {
                val_nuc = 2;
                break;
            }
            case 't': {
                val_nuc = 1;
                break;
            }
            default: {
                val_nuc = 0;
                break;
            }
        }

        //Add new nucleotide value
        h <<= 2;
        h |= val_nuc;
        i++;
        if (i >= wsize)
            kmer_pos_list.emplace_back(std::pair<unsigned, unsigned>(h,i - wsize));
        //Suppress old boundary value
        unsigned old_val = h & mask;
        h ^= old_val;
        nb_kmer++;
        s++;
    }
}
//-------------------------------------------------------------------------
void HashDNASeq::hashSeqPosWHole(const BioSeq &seq, unsigned wsize, std::list<std::pair<unsigned, unsigned> > &kmer_pos_list) {

    unsigned len=seq.length();
    if(len<=wsize) return;
    unsigned last_pos=len-wsize;
    const char *s=seq.c_str();
    for(unsigned i=0;i<=last_pos;i++)
    {
        if (i >= wsize)
            kmer_pos_list.emplace_back(std::pair<unsigned, unsigned>(hseq(s),i - wsize));
        s++;
    }
}
//-------------------------------------------------------------------------
void HashDNASeq::search(const BioSeq &sequence, unsigned numseq, unsigned start, unsigned end, bool repeat,
                        std::vector<std::pair<unsigned, unsigned> > &frag, unsigned verbose)
{
	clock_t clock_begin, clock_end;
    if(verbose>0) {
        clock_begin = clock();
        std::cout << "hashing query sequence..." << std::flush;
    }
	Diag_map diag_map(nbseqS);

    if(hash_algorithm==1)
        matchKmersHole(sequence, numseq, start, end, repeat, diag_map);
    else if(hash_algorithm==0)
        matchKmers(sequence, numseq, start, end, repeat, diag_map);
    else if(hash_algorithm==2)
        matchKmersMinimizer(sequence, numseq, start, end, repeat, diag_map);
    else if(hash_algorithm==3)
        matchKmersWHoleMinimizer(sequence, numseq, start, end, repeat, diag_map);
    if(verbose>0) {
        std::cout << "ok" << std::endl;
        std::cout << diag_map.size() << " hits found";
        clock_end = clock();
        std::cout << " --> Time spent: " << (double) (clock_end - clock_begin) / CLOCKS_PER_SEC << " seconds"
                  << std::endl;

        clock_begin = clock();
        std::cout << "search fragments..." << std::flush;
    }
    diagSearchDist(diag_map,(kmer_size+1)*wdist,kmer_size,frag);
	diag_map.clear();
    if(verbose>0) {
        std::cout << "ok" << std::endl;
        std::cout << frag.size() << " fragments found";
        clock_end = clock();
        std::cout << " --> Time spent: " << (double) (clock_end - clock_begin) / CLOCKS_PER_SEC << " seconds"
                  << std::endl;
    }
}
//-------------------------------------------------------------------------
// Search for alignments with kmer matches
void HashDNASeq::matchKmersHole(const BioSeq& sequence,unsigned qSeq,
                                unsigned start, unsigned end, bool repeat,
                                Diag_map& diag_map)
{
    unsigned last_pos=end-kmer_size;
    if(end<=kmer_size) return;

    std::string str=sequence.substr(start,end-start+1);
    const char* seq=str.c_str();

    unsigned key_d,dirhit=0;
    unsigned i=start;
    while(i<=last_pos) {
        bool found=false;
        key_d = hseq(seq);
        auto begin_d = hash2wpos[key_d];
        auto end_d = hash2wpos[key_d + 1];
        for (auto j = begin_d; j != end_d; j++) {
            unsigned sSeq=j->numSeq;
            unsigned sPos=j->pos;
            if (sSeq == 0) continue;
            if (sSeq > 0) {
                long diag = long(i) - sPos;
                if (!repeat || (repeat && i < sPos && qSeq<=sSeq)){
                    dirhit++;
                    diag_map[sSeq].emplace_back(Diag(diag, sPos, sSeq));
                    found = true;
                }
            }
        }
        if(found){
            seq += step_q;
            i+=step_q;
        } else {
            seq += 1;
            i += 1;
        }
    }
    std::cout<<dirhit<<" hits found / ";
}
//-------------------------------------------------------------------------
// Search for alignments with kmer matches
void HashDNASeq::matchKmers(const BioSeq& sequence, unsigned qSeq,
                            unsigned start, unsigned end, bool repeat,
                            Diag_map& diag_map)
{
    unsigned dirhit=0;

    std::list<std::pair<unsigned,unsigned>> kmer_pos_list;
    BioSeq subseq=sequence.substr(start,end-start+1);
    hashSeqPos(subseq,kmer_size,kmer_pos_list);

    for(auto &iter_pos : kmer_pos_list){
        unsigned pos=start+iter_pos.second;
        unsigned key_d=iter_pos.first;

        auto begin_d = hash2wpos[key_d];
        auto end_d = hash2wpos[key_d + 1];
        for (auto j = begin_d; j != end_d; j++) {
            unsigned sSeq=j->numSeq;
            unsigned sPos=j->pos;
            if (sSeq == 0) continue;
            if (sSeq > 0) {
                long diag = long(pos) - sPos;
                if (!repeat || (repeat && pos <sPos && qSeq<=sSeq)){
                    dirhit++;
                    diag_map[sSeq].emplace_back(Diag(diag, sPos, sSeq));
                }
            }
        }
    }
    std::cout<<dirhit<<" hits found / ";
}
//-------------------------------------------------------------------------
// Search for alignments with kmer matches
void HashDNASeq::matchKmersMinimizer(const BioSeq& sequence, unsigned qSeq,
                            unsigned start, unsigned end, bool repeat,
                                     Diag_map& diag_map)
{
    unsigned dirhit=0;

    std::list<std::pair<unsigned,unsigned>> kmer_pos_list;
    std::set<std::pair<unsigned,unsigned>> minimized_kmer_pos_list;

    BioSeq subseq=sequence.substr(start,end-start+1);
    hashSeqPos(subseq,kmer_size,kmer_pos_list);
    minimize(window_size,kmer_pos_list,minimized_kmer_pos_list);

    for(auto &iter_pos : minimized_kmer_pos_list){
        unsigned pos=start+iter_pos.second;
        unsigned key_d=iter_pos.first;

        auto begin_d = hash2wpos[key_d];
        auto end_d = hash2wpos[key_d + 1];
        for (auto j = begin_d; j != end_d; j++) {
            unsigned sSeq=j->numSeq;
            unsigned sPos=j->pos;
            if (sSeq == 0) continue;
            if (sSeq > 0) {
                long diag = long(pos) - sPos;
                if (!repeat || (repeat && pos < sPos && qSeq<=sSeq)){
                    dirhit++;
                    diag_map[sSeq].emplace_back(Diag(diag, sPos, sSeq));
                }
            }
        }
    }
    std::cout<<dirhit<<" hits found / ";
}
//-------------------------------------------------------------------------
// Search for alignments with kmer matches
void HashDNASeq::matchKmersWHoleMinimizer(const BioSeq& sequence, unsigned qSeq,
                                     unsigned start, unsigned end, bool repeat,
                                     Diag_map& diag_map)
{
    unsigned dirhit=0;

    std::list<std::pair<unsigned,unsigned>> kmer_pos_list;
    std::set<std::pair<unsigned,unsigned>> minimized_kmer_pos_list;

    BioSeq subseq=sequence.substr(start,end-start+1);
    hashSeqPosWHole(subseq,kmer_size,kmer_pos_list);
    minimize(window_size,kmer_pos_list,minimized_kmer_pos_list);

    for(auto &iter_pos : minimized_kmer_pos_list){
        unsigned pos=start+iter_pos.second;
        unsigned key_d=iter_pos.first;

        auto begin_d = hash2wpos[key_d];
        auto end_d = hash2wpos[key_d + 1];
        for (auto j = begin_d; j != end_d; j++) {
            unsigned sSeq=j->numSeq;
            unsigned sPos=j->pos;
            if (sSeq == 0) continue;
            if (sSeq > 0) {
                long diag = long(pos) - sPos;
                if (!repeat || (repeat && pos < sPos && qSeq<=sSeq)){
                    dirhit++;
                    diag_map[sSeq].emplace_back(Diag(diag, sPos, sSeq));
                }
            }
        }
    }
    std::cout<<dirhit<<" hits found / ";
}
//-------------------------------------------------------------------------
void HashDNASeq::minimize(unsigned window_size, std::list<std::pair<unsigned, unsigned>> kmer_pos_list,
                          std::set<std::pair<unsigned, unsigned>> &minimized_kmer_pos_list) {

    std::set<std::pair<unsigned,unsigned>> minimizer;

    for(auto &iter_pos : kmer_pos_list){
        minimizer.emplace(iter_pos);
        if(minimizer.size()==window_size){
            minimized_kmer_pos_list.emplace(*minimizer.begin());
            std::set<std::pair<unsigned,unsigned>>::iterator k=minimizer.begin();
            while(k != minimizer.end()){
                if(k->second==iter_pos.second-(window_size-1)){
                    k=minimizer.erase(k);
                }
                else{k++;}
            }
        }
    }
}
//-------------------------------------------------------------------------
// Search for diagonal of kmer matches
void HashDNASeq::diagSearchDist(Diag_map & diag_map,
                                unsigned connect_dist, unsigned kmerSize,
                                std::vector< std::pair<unsigned,unsigned> >& frag)
{
    for (auto &iter_seq : diag_map) { // iter seq
        if (iter_seq.size() > 2) {
            iter_seq.sort();
            auto iter_diag = iter_seq.begin();
            Diag prev_d = *iter_diag;
            unsigned start = 0;
            unsigned end = 0;
            long diag = 0;
            while (++iter_diag != iter_seq.end()) {
                Diag curr_d = *iter_diag;
                if (prev_d.diag == curr_d.diag // Same diagonal
                    && prev_d.wpos.numSeq == curr_d.wpos.numSeq // Same sequence
                    &&
                    prev_d.wpos.pos + connect_dist >= curr_d.wpos.pos) // hits overlaps or close enough to be joined
                {
                    if (start != 0) // already extending a diagonal
                    {
                        end = curr_d.wpos.pos;
                    } else //create a diagonal
                    {
                        diag = prev_d.diag;
                        start = prev_d.wpos.pos;
                        end = curr_d.wpos.pos;
                    }
                } else // too far to be joined
                {
                    if (start != 0) // Create a diagonal
                    {
                        frag.emplace_back(diag + start + 1, diag + end + kmerSize);
                        start = 0;
                    };
                }
                prev_d = curr_d;
            }
            if (start != 0)  //End of hits list
            {
                frag.emplace_back(diag + start + 1, diag + end + kmerSize);
            }
        }
    }
}





