#include "Duster.h"

//-------------------------------------------------------------------------
void Duster::load(const SDGString& filenameS, unsigned kmer_size, unsigned kmask, unsigned bkmer_size, unsigned mkmer_size, double count_cutoff, double diversity_cutoff,
		unsigned min_count,bool & valid_idx_file, bool first_iter)
{

  nbseqS=0;

  clock_t begin, end;
  begin = clock();

  std::cout<<"Effective Kmer size (without kmer mask size):"<<hseq.getEffectiveKmerSize()<<std::endl;
  if(valid_idx_file)
	  valid_idx_file=read_idx(filenameS,count_cutoff,diversity_cutoff, min_count,kmask);
  if(!valid_idx_file)
    {
	  //Count kmers and filters
	  std::vector<unsigned> kmer_count((unsigned)pow(4,hseq.getEffectiveKmerSize()),0);
	  unsigned nb_kmer=0;
	  std::list< Info_kmer > list_infokmer;
	  Info_kmer kmer_threshold;
	  kmer_analysis(filenameS,kmer_size, kmask, bkmer_size, mkmer_size, count_cutoff, diversity_cutoff, kmer_count, nb_kmer, list_infokmer,kmer_threshold);
      kmer_filter(list_infokmer, kmer_threshold, min_count, kmer_count, first_iter);
      end = clock();
      std::cout<<" --> Time spent: "<<(double)(end-begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;


      //Prepare hash table of kmer positions
      begin = clock();
      std::cout<<"Prepare hash table"<<std::endl;
      kmer_pos.resize(nb_kmer);
      unsigned k=0;
      std::vector<KmerSpos>::iterator last_it=kmer_pos.begin();
      for(std::vector<unsigned>::iterator i=kmer_count.begin(); i!=kmer_count.end();i++)
		{
		  hash2wpos[k]=last_it;
		  last_it=last_it+(*i);
		  k++;
		}
      hash2wpos[k]=last_it;

      //Load hash table with kmer positions
      std::cout<<"Load kmer positions in hash table"<<std::endl;
      hash_ptr=hash2wpos;
      SDGFastaIstream inS2(filenameS);
      while(inS2)
		{
		  SDGBioSeq sS;
		  if(inS2)
			inS2>>sS;
		  hashSeqPos(sS,kmer_count);
		}
      inS2.close();
      hash_ptr.clear();
      save_idx(filenameS,count_cutoff,diversity_cutoff,min_count,kmask,kmer_count);
    }
  std::cout<<"end preparation of subjects";
  end = clock();
  std::cout<<" --> Time spent: "<<(double)(end-begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
}
//-------------------------------------------------------------------------
void Duster::kmer_analysis(const SDGString& filenameS, unsigned kmer_size, unsigned mask, unsigned bkmer_size, unsigned mkmer_size,
		double count_cutoff, double diversity_cutoff, std::vector<unsigned>& kmer_count, unsigned& nb_kmer,
		std::list< Info_kmer >& list_infokmer, Info_kmer& kmer_threshold)
{

	  std::vector<unsigned> background_count((unsigned)pow(4,bkmer_size),0);
	  std::vector<unsigned> model_count((unsigned)pow(4,mkmer_size),0);
	  std::vector<unsigned> nuc_count(4,0);

  unsigned nb_bkmer=0;
  unsigned nb_mkmer=0;
  unsigned nb_nuc=0;
  nb_kmer=0;

  kmer_counts(filenameS,kmer_size,bkmer_size,mkmer_size,
		  kmer_count, nb_kmer,
		  background_count, nb_bkmer,
		  model_count, nb_mkmer,
		  nuc_count, nb_nuc);



  kmer_count[0]=0; //remove kmers AAAAAA... NNNNNN.... XXXXX....
  kmer_prob(kmer_size, bkmer_size, mkmer_size,mask,
		  kmer_count,nb_kmer,
		  background_count, nb_bkmer,
		  model_count, nb_mkmer,
		  nuc_count, nb_nuc,
		  list_infokmer);

  kmer_count_percentiles(list_infokmer, count_cutoff, kmer_threshold);

  double cutoff_entropy=1.0;
  kmer_entropy_percentiles(list_infokmer, cutoff_entropy, kmer_threshold);

  kmer_diversity_percentiles(list_infokmer, diversity_cutoff, kmer_threshold);
  kmer_goodkmer_percentiles(list_infokmer);

}
//-------------------------------------------------------------------------
void Duster::kmer_counts(const SDGString& filenameS, unsigned kmer_size, unsigned bkmer_size, unsigned mkmer_size,
		std::vector<unsigned>& wcount, unsigned& nb_kmer,
		std::vector<unsigned>& bcount, unsigned& nb_bkmer,
		std::vector<unsigned>& mcount, unsigned& nb_mkmer,
		std::vector<unsigned>& ncount, unsigned& nb_nuc)
{
  std::cout<<"Counting kmer of size "<<kmer_size<<" ... "<<std::flush;

  SDGFastaIstream inS(filenameS);
  unsigned count_seq=0;
  while(inS)
	{
	  SDGBioSeq sS;
	  if(inS)
		inS>>sS;
	  std::cout<<"\n"<<++count_seq<<"->"<<sS.getDE()<<std::flush;
	  nb_kmer+=hashSeqCount(sS,kmer_size, wcount);
	  nb_bkmer+=hashSeqBackgroundCount(sS,bkmer_size, bcount);
	  nb_mkmer+=hashSeqModelCount(sS,mkmer_size, mcount);
	  nb_nuc+=hashSeqNucCount(sS, ncount);
	}
  inS.close();
  std::cout<<"\nNumber of kmers="<<nb_kmer<<std::endl;
}
//-------------------------------------------------------------------------
//
void Duster::kmer_count_percentiles(const std::list< Info_kmer >& list_infokmer, double cutoff,
		Info_kmer& kmer_threshold)
{
	std::cout<<"\nCount percentiles:"<<std::endl;

    std::vector<unsigned> distr;
	for(std::list< Info_kmer >::const_iterator i=list_infokmer.begin(); i!=list_infokmer.end(); i++)
    	distr.push_back(i->count);
    sort(distr.begin(),distr.end());

    std::vector<unsigned> clean_distr;
    copy(lower_bound(distr.begin(),distr.end(),(unsigned)1), distr.end(),back_inserter(clean_distr));
    std::cout<<" min="<<clean_distr.front()<<std::endl;
    std::cout<<" 1st quartile="
	       << clean_distr[(unsigned)floor(0.25*clean_distr.size())]<<std::endl;
    std::cout<<" median="
	       << clean_distr[(unsigned)floor(0.5*clean_distr.size())]<<std::endl;
    std::cout<<" 3st quartile="
	       << clean_distr[(unsigned)floor(0.75*clean_distr.size())]<<std::endl;
    std::cout<<" 95 percentile="
	       << clean_distr[(unsigned)floor(0.95*clean_distr.size())]<<std::endl;
    std::cout<<" 99 percentile="
	       << clean_distr[(unsigned)floor(0.99*clean_distr.size())]<<std::endl;
    std::cout<<" max="<< clean_distr.back()<<std::endl;
    long threshold_index=(long)(cutoff*distr.size())-1;
    if(threshold_index<0) threshold_index=0;
    kmer_threshold.count=distr[threshold_index];

	unsigned count=0;
	for(std::vector<unsigned>::iterator it=distr.begin(); it!=distr.end(); it++)
	{
		if(*it>kmer_threshold.count)
			break;
		count++;
	}
    std::cout<<"=>cut-off="<<cutoff<<". Kmers occuring more than "<<kmer_threshold.count
	       <<" will be removed! Will remove "<<count<<" kmers."<<std::endl;
}
//-------------------------------------------------------------------------
//
void Duster::kmer_entropy_percentiles(const std::list< Info_kmer >& list_infokmer, double cutoff_entropy,
		Info_kmer& kmer_threshold)
{

	std::cout<<"\nEntropy percentiles:"<<std::endl;

	std::vector<double> distr;
	for(std::list< Info_kmer >::const_iterator i=list_infokmer.begin(); i!=list_infokmer.end(); i++)
    	distr.push_back(i->entropy);
    sort(distr.begin(),distr.end());

    std::cout<<" min="<<distr.front()<<std::endl;
    std::cout<<" 1st quartile="
	       << distr[(unsigned)floor(0.25*distr.size())]<<std::endl;
    std::cout<<" median="
	       << distr[(unsigned)floor(0.5*distr.size())]<<std::endl;
    std::cout<<" 3st quartile="
	       << distr[(unsigned)floor(0.75*distr.size())]<<std::endl;
    std::cout<<" 95 percentile="
	       << distr[(unsigned)floor(0.95*distr.size())]<<std::endl;
    std::cout<<" 99 percentile="
	       << distr[(unsigned)floor(0.99*distr.size())]<<std::endl;
    std::cout<<" max="<< distr.back()<<std::endl;
    long threshold_index=(long)(cutoff_entropy*distr.size())-1;
    if(threshold_index<0) threshold_index=0;
    kmer_threshold.entropy=distr[threshold_index];

	unsigned count=0;
	for(std::vector<double>::iterator it=distr.begin(); it!=distr.end(); it++)
	{
		if(*it>kmer_threshold.entropy)
			break;
		count++;
	}
    std::cout<<"=>cut-off="<<cutoff_entropy<<". Kmers with entropy greater than "<<kmer_threshold.entropy
	       <<" will be removed! Will remove "<<count<<" kmers."<<std::endl;
}
//-------------------------------------------------------------------------
//
void Duster::kmer_diversity_percentiles(const std::list< Info_kmer >& list_infokmer, double diversity_cutoff,
		Info_kmer& kmer_threshold)
{

	std::cout<<"\nDiversity percentiles:"<<std::endl;

	std::vector<double> distr;
	for(std::list< Info_kmer >::const_iterator i=list_infokmer.begin(); i!=list_infokmer.end(); i++)
    	distr.push_back(i->diversity);
    sort(distr.begin(),distr.end());

    std::cout<<" min="<<distr.front()<<std::endl;
    std::cout<<" 1st quartile="
	       << distr[(unsigned)floor(0.25*distr.size())]<<std::endl;
    std::cout<<" median="
	       << distr[(unsigned)floor(0.5*distr.size())]<<std::endl;
    std::cout<<" 3st quartile="
	       << distr[(unsigned)floor(0.75*distr.size())]<<std::endl;
    std::cout<<" 95 percentile="
	       << distr[(unsigned)floor(0.95*distr.size())]<<std::endl;
    std::cout<<" 99 percentile="
	       << distr[(unsigned)floor(0.99*distr.size())]<<std::endl;
    std::cout<<" max="<< distr.back()<<std::endl;


	kmer_threshold.diversity=diversity_cutoff;
	unsigned count=0;
	for(std::vector<double>::iterator it=distr.begin(); it!=distr.end(); it++)
	{
		if(*it>kmer_threshold.diversity)
			break;
		count++;
	}
	std::cout<<"=>diversity_cutoff="<<kmer_threshold.diversity<<". Kmers with diversity less than "<<kmer_threshold.diversity
		   <<" will be removed! Will remove "<<count<<" kmers."<<std::endl;

}
//-------------------------------------------------------------------------
//
void Duster::kmer_goodkmer_percentiles(const std::list< Info_kmer >& list_infokmer)
{

	std::cout<<"\nGood kmer percentiles:"<<std::endl;

	std::vector<double> distr;
	for(std::list< Info_kmer >::const_iterator i=list_infokmer.begin(); i!=list_infokmer.end(); i++)
    	distr.push_back(i->good_kmer);
    sort(distr.begin(),distr.end());

    std::cout<<" min="<<distr.front()<<std::endl;
    std::cout<<" 1st quartile="
	       << distr[(unsigned)floor(0.25*distr.size())]<<std::endl;
    std::cout<<" median="
	       << distr[(unsigned)floor(0.5*distr.size())]<<std::endl;
    std::cout<<" 3st quartile="
	       << distr[(unsigned)floor(0.75*distr.size())]<<std::endl;
    std::cout<<" 95 percentile="
	       << distr[(unsigned)floor(0.95*distr.size())]<<std::endl;
    std::cout<<" 99 percentile="
	       << distr[(unsigned)floor(0.99*distr.size())]<<std::endl;
    std::cout<<" max="<< distr.back()<<std::endl;


	unsigned count=0;
	for(std::vector<double>::iterator it=distr.begin(); it!=distr.end(); it++)
	{
		if(*it>0.0)
			break;
		count++;
	}
	std::cout<<"Kmers with poor model probability will be removed (if not one iteration)! Will remove "<<count<<" kmers."<<std::endl;

}
//-------------------------------------------------------------------------
//
void Duster::kmer_prob(unsigned wsize, unsigned bwsize,unsigned mwsize, unsigned mask,
		const std::vector<unsigned>& wcount, unsigned nb_kmer,
		const std::vector<unsigned>& bcount, unsigned nb_bkmer,
		const std::vector<unsigned>& mcount, unsigned nb_mkmer,
		const std::vector<unsigned>& ncount, unsigned nb_nuc,
		std::list< Info_kmer >& list_infokmer)
{
	std::cout<<"\nCompute kmer stats ... "<<std::endl;
	Duster h(wsize,mask,bwsize, 0);
	std::cout<<"kmer size="<<wsize
			<<"\tmask size="<<mask
			<<"\teffective kmer size="<<h.hseq.getEffectiveKmerSize()<<std::endl;

	std::cout<<"\nSequence composition:\n";

	//Compute nucleotide probabilities
	unsigned max_key_bkmer=pow(4,bwsize);
	unsigned max_key_mkmer=pow(4,mwsize);

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
    	list_prob.push_back(std::make_pair(bprob[k],k));

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
    for(std::list< std::pair<double,unsigned> >::reverse_iterator i=list_prob.rbegin()  ; i!=list_prob.rend() && prob2display < maxProb2display; i++)
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
       	double diversity=list_skmer.size()/pow(4,bkmer_size);

       	s=kmer.c_str();

    	double prob_mkmer=mprob[h.mhseq(s++)];
    	for(unsigned i=1; i<=wsize-mwsize && *s!='\0'; i++)
    	{
    		prob_mkmer*=mprob_cond[h.mhseq(s++)];
    	}
    	double belong2model=log10(prob_mkmer)-log10(prob_bkmer);

    	list_infokmer.push_back(Info_kmer(k,wcount[k],ecount,entropy*-1,diversity,belong2model));
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
    for(std::list< Info_kmer >::reverse_iterator i=list_infokmer.rbegin(); i!=list_infokmer.rend() && count2display < maxCount2display; i++)
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
    for(std::list< Info_kmer >::iterator i=list_infokmer.begin(); i!=list_infokmer.end() && count2display < maxCount2display; i++)
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
void Duster::kmer_filter(const std::list< Info_kmer >& list_infokmer, const Info_kmer& kmer_threshold, unsigned min_count,
			std::vector<unsigned>& wcount, bool first_iter)
{
    std::cout<<"\nFilter kmers ... "<<std::flush;
    unsigned nb_kmer=0;
    unsigned nb_removed=1; // 1 for AAAAAA...
    for(std::list< Info_kmer >::const_iterator i=list_infokmer.begin();i!=list_infokmer.end();i++)
    {
    	unsigned k=i->hash_key;
    	if(wcount[k]==0) continue;
    	nb_kmer++;
    	if(i->count<=min_count*(i->expected_count))
    	{
    		wcount[k]=0;
    		nb_removed++;
    	}
    	if(i->count>=kmer_threshold.count)
    	{
    		wcount[k]=0;
    		nb_removed++;
    	}
    	if(i->entropy>=kmer_threshold.entropy)
    	{
    		wcount[k]=0;
    		nb_removed++;
    	}
    	if(i->diversity<=kmer_threshold.diversity)
    	{
    		wcount[k]=0;
    		nb_removed++;
    	}
    	if(!first_iter && i->good_kmer<=0.0)
    	{
    		wcount[k]=0;
    		nb_removed++;
    	}
    }

    std::cout<<"removed kmers:"<<nb_removed<<" over "<<nb_kmer<<" existing ("<<(double)nb_removed/nb_kmer<<") and "<<max_key+1<<" possible."<<std::endl;
}
//-------------------------------------------------------------------------
// Read a hash index
bool Duster::read_idx(const SDGString& filename, double count_cutoff , double diversity_cutoff, unsigned min_count, unsigned kmask)
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
  if(val!=min_count) return false;
  double fval;
  in>>fval;
  if(fval!=count_cutoff) return false;
  in>>fval;
  if(fval!=diversity_cutoff) return false;

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
  std::vector<KmerSpos>::iterator last_it=kmer_pos.begin();
  for(std::vector<unsigned>::iterator i=kmer_count.begin(); i!=kmer_count.end();i++)
    {
      hash2wpos[k]=last_it;
      last_it=last_it+(*i);
      k++;
    }
  hash2wpos[k]=last_it;

  return true;
}
//-------------------------------------------------------------------------
// Save a hash index
void Duster::save_idx(const SDGString& filename, double count_cutoff, double diversity_cutoff, unsigned min_count, unsigned kmask, const std::vector<unsigned>& wcount)
{
  std::ofstream fout(filename+".kidx");
  std::cout<<"Save idx"<<std::endl;
  fout<<kmer_size<<std::endl;
  fout<<kmask<<std::endl;
  fout<<min_count<<std::endl;
  fout<<count_cutoff<<std::endl;
  fout<<diversity_cutoff<<std::endl;
  fout<<kmer_pos.size()<<std::endl;
  for(std::vector<KmerSpos>::iterator i=kmer_pos.begin();i!=kmer_pos.end();i++)
    {
      fout<<i->pos<<" "<<i->numSeq<<"\t";
    }
  fout<<wcount.size()<<std::endl;
  for(std::vector<unsigned>::const_iterator i=wcount.begin();
      i!=wcount.end();i++)
    {
      fout<<*i<<"\t";
    }
}
//-------------------------------------------------------------------------
// Count kmers
unsigned Duster::hashSeqCount(const SDGBioSeq& seq, unsigned wsize, std::vector<unsigned>& wcount)
{
  unsigned len=seq.length();
  unsigned nb_kmer=0;
  if(len<=wsize) return 0;
  unsigned last_pos=len-wsize;
  const char* s=seq.toString().c_str();
  for(unsigned i=0;i<=last_pos;i++)
    {
	  nb_kmer++;
      wcount[hseq(s)]++;
      s++;
    }
  return nb_kmer;
}
//-------------------------------------------------------------------------
// Count background kmers
unsigned Duster::hashSeqBackgroundCount(const SDGBioSeq& seq, unsigned wsize, std::vector<unsigned>& wcount)
{
  unsigned len=seq.length();
  unsigned nb_kmer=0;
  if(len<=wsize) return 0;
  unsigned last_pos=len-wsize;
  const char* s=seq.toString().c_str();
  for(unsigned i=0;i<=last_pos;i++)
    {
	  nb_kmer++;
      wcount[bhseq(s)]++;
      s++;
    }
  return nb_kmer;
}
//-------------------------------------------------------------------------
// Count background kmers
unsigned Duster::hashSeqModelCount(const SDGBioSeq& seq, unsigned wsize, std::vector<unsigned>& wcount)
{
  unsigned len=seq.length();
  unsigned nb_kmer=0;
  if(len<=wsize) return 0;
  unsigned last_pos=len-wsize;
  const char* s=seq.toString().c_str();
  for(unsigned i=0;i<=last_pos;i++)
    {
	  nb_kmer++;
      wcount[mhseq(s)]++;
      s++;
    }
  return nb_kmer;
}
//-------------------------------------------------------------------------
// Count nucleotides
unsigned Duster::hashSeqNucCount(const SDGBioSeq& seq, std::vector<unsigned>& wcount)
{
  unsigned len=seq.length();
  unsigned nb_kmer=0;
  if(len<1) return 0;
  const char* s=seq.toString().c_str();
  for(unsigned i=0;i<len;i++)
    {
	  nb_kmer++;
      wcount[nhseq(s)]++;
      s++;
    }
  return nb_kmer;
}
//-------------------------------------------------------------------------
void Duster::hashSeqPos(const SDGBioSeq& seq, const std::vector<unsigned>& wcount)
{
  nbseqS++;
  unsigned len=seq.length();
  if(len<=kmer_size) return;
  unsigned last_pos=len-kmer_size;
  unsigned key;
  const char* s=seq.toString().c_str();
  for(unsigned i=0;i<=last_pos;i++)
    {
      key=hseq(s);

      if(wcount[key]!=0)
		{
		  *(hash_ptr[key])=KmerSpos(i,nbseqS);
		  hash_ptr[key]++;
		}
      s++;
    }
}
//-------------------------------------------------------------------------
void Duster::search(const SDGBioSeq& sequence, unsigned start, unsigned end, unsigned numseq, bool repeat, std::vector< std::pair<unsigned,unsigned> >& fmerged)
{
	clock_t clock_begin, clock_end;
	clock_begin = clock();
	std::cout<<"hashing query sequence..."<<std::flush;
	std::vector< Diag > diag_map;
	matchKmers(sequence, start, end, numseq, repeat, diag_map);
	std::cout<<"ok"<<std::endl;
	std::cout<<diag_map.size()<<" hits found";
	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;

	clock_begin = clock();
	std::cout<<"search fragments..."<<std::flush;
	std::vector< std::pair<unsigned,unsigned> > frag;
	diagSearch(diag_map,(wdist+1)*kmer_size,kmer_size,frag);
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
// Search for alignments with kmer matches
void Duster::matchKmers(const SDGBioSeq& sequence,
		unsigned start, unsigned end, unsigned numseq, bool repeat,
		std::vector< Diag >& diag_map)
{
  diag_map.clear();

  std::string str=sequence.toString().substr(start,end-start);

  const char* seq=str.c_str();

  unsigned last_pos=end-kmer_size;
  if(end<=kmer_size) return;

  unsigned key_d=0;
  for(unsigned i=start;i<=last_pos;i+=step_q)
    {
      key_d=hseq(seq);
      std::vector<KmerSpos>::iterator begin_d=hash2wpos[key_d];
      std::vector<KmerSpos>::iterator end_d=hash2wpos[key_d+1];
      for(std::vector<KmerSpos>::iterator j=begin_d;j!=end_d;j++)
      {
    	  // Attention cas j->numSeq==0 à traiter (no kmer hit!)
    	  if(!(repeat && numseq==j->numSeq && i==j->pos))
    			  diag_map.push_back(Diag(i-j->pos,j->pos,j->numSeq));
      }
      seq+=step_q;
    }
}
//-------------------------------------------------------------------------
// Search for diagonal of kmer matches
void Duster::diagSearch(std::vector< Diag >& diag_map,
		unsigned connect_dist, unsigned kmer_size,
		std::vector< std::pair<unsigned,unsigned> >& frag)
{
  unsigned size=diag_map.size();
  if(size>=2)
    {
      sort(diag_map.begin(),diag_map.end());

      unsigned start=0;
      unsigned end=0;
      int diag=0;

      Diag& prev_d=diag_map[0];
      //std::cout<<prev_d<<std::endl;
      for( unsigned i=1; i<size; ++i)
		{
		  Diag& curr_d=diag_map[i];
	      //std::cout<<curr_d<<std::endl;
		  if(prev_d.diag==curr_d.diag // Same diagonal
			 && prev_d.wpos.numSeq==curr_d.wpos.numSeq // Same sequence
			 && prev_d.wpos.pos+connect_dist>=curr_d.wpos.pos) // hits overlaps or close enough to be joined
				{
				  if(start!=0) // already extending a diagonal
					{
				      //std::cout<<"Extend "<<prev_d<<std::endl;
					  end=curr_d.wpos.pos;
					}
				  else //create a diagonal
					{
					  //std::cout<<"Create "<<prev_d<<".."<<curr_d<<std::endl;
					  diag=prev_d.diag;
					  start=prev_d.wpos.pos;
					  end=curr_d.wpos.pos;
					}
				}
			  else // too far to be joined
				  if(start!=0) // Create a diagonal
					{
					  //std::cout<<"Diagonal found is "<<diag+start+1<<"-"<<diag+end+kmer_size<<std::endl;
					  frag.push_back(std::pair<unsigned,unsigned>(diag+start+1,diag+end+kmer_size));
					  start=0;
					}
		  prev_d=curr_d;
		}
      if(start!=0)  //End of hits list
		{
		  	  //std::cout<<"Diagonal found is "<<diag+start+1<<"-"<<diag+end+kmer_size<<std::endl;
			  frag.push_back(std::pair<unsigned,unsigned>(diag+start+1,diag+end+kmer_size));
		}
    }
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
		unsigned curr_start=0;
		unsigned curr_end=0;
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
					fmerged.push_back(std::pair<unsigned,unsigned>(start,end));
					start=0;
				  }
				  else
				  {
					  if(prev_end-prev_start>=min_size && prev_start!=0)
						  fmerged.push_back(std::pair<unsigned,unsigned>(prev_start,prev_end));
					  start=0;
				  }
			  }
			prev_start=curr_start;
			if(curr_end>end) prev_end=curr_end;
			else if (start!=0) prev_end=end;
			else prev_end=curr_end;
		}
		if(start!=0) //end of sequence reached, must save the match range if one is in extension
		{
			if(end-start>=min_size)
				fmerged.push_back(std::pair<unsigned,unsigned>(start,end));
		}
		else
		{
			if(prev_end-prev_start>=min_size)
				fmerged.push_back(std::pair<unsigned,unsigned>(prev_start,prev_end));
		}
      }
}
//-------------------------------------------------------------------------
void Duster::writeHitBED(SDGString qname, const std::vector< std::pair<unsigned,unsigned> >& frag, std::ostream& out)
{
  unsigned size=frag.size();
  for(unsigned i=0; i<size; ++i)
    {
	  out<<qname<<"\t" // chromosome
		 <<frag[i].first  //chrom start
		 <<"\t"<<frag[i].second  //chrom end
		 <<"\t"<<"duster" // name
		 <<"\t"<<frag[i].second-frag[i].first //score
		 <<"\t"<<"+" //strand
//		 <<"\t"<<frag[i].first //thickStart
//		 <<"\t"<<frag[i].second //thickEnd
		 <<std::endl;
    }
}//-------------------------------------------------------------------------
void Duster::writeBED(SDGString qname, const std::vector< std::pair<unsigned,unsigned> >& frag, std::ostream& out)
{
  unsigned size=frag.size();
  for(unsigned i=0; i<size; ++i)
    {
	  out<<qname<<"\t" // chromosome
		 <<frag[i].first  //chrom start
		 <<"\t"<<frag[i].second  //chrom end
		 <<"\t"<<"duster" // name
		 <<"\t"<<frag[i].second-frag[i].first //score
		 <<"\t"<<"+" //strand
		 <<std::endl;
    }
}
//-------------------------------------------------------------------------
unsigned Duster::compute_coverage(const std::vector< std::pair<unsigned,unsigned> >& frag)
{
  unsigned size=frag.size();
  unsigned coverage=0;
  for(unsigned i=0; i<size; ++i)
    {
	  coverage+=frag[i].second-frag[i].first;
    }
  return coverage;
}
//-------------------------------------------------------------------------
void Duster::get_sequences(const std::vector< std::pair<unsigned,unsigned> >& frag, SDGBioSeq& seq, SDGFastaOstream& out)
{
  unsigned size=frag.size();
  for(unsigned i=0; i<size; ++i)
    {
	  SDGBioSeq subseq=seq.subseq(frag[i].first,frag[i].second-frag[i].first+1);
	  std::ostringstream name;
	  name<<seq.getDE()<<":"<<frag[i].first<<".."<<frag[i].second;
	  subseq.setDE((SDGString)name.str());
	  out<<subseq;
    }
}





