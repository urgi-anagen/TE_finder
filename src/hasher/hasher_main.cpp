#include <getopt.h>
#include <fstream>
#include <cstdlib>
#include <list>

#include "SDGString.h"
#include "SDGFastaIstream.h"
#include "SDGFastaOstream.h"
#include "SDGMemBioSeq.h"

#include "Hasher.h"


double filter_cutoff=0.0;
unsigned kmer_size=15, step_q=15, bkmer_size=1, kmer_dist=5,
min_size=20, min_frag_size, chunk_size_kb=0, min_count=0, kmask=4, verbosity=0, overlap=0;
double count_cutoff=1.0, diversity_cutoff=0.0;
bool repeat=false, stat_only=false;

SDGString outfilename="";

    
void help(void)
{
  std::cerr<<"usage: hasher"
      <<" [<options>] <fasta query sequence> [<fasta sequence model>]. Note: without sequence model, duster will perform its own repeat search."<<std::endl
      <<" options:"<<std::endl
      <<"   -h, --help:\n\t this help"<<std::endl
      <<"   -w, --kmer:\n\t kmer length, default: "<<kmer_size<<std::endl
      <<"   -S, --step_q:\n\t step on query sequence, default: "<<step_q<<std::endl
      <<"   -k, --kmask:\n\t period of k-mer mask, default: "<<kmask<<std::endl
      <<"   -d, --kmer_dist:\n\t max number of kmer between two matching kmer to connect, default: "
	<<kmer_dist<<std::endl
    <<"   -s, --min_size:\n\t min size range to report, default: "
	<<min_size<<std::endl
    <<"   -C, --filter_cutoff:\n\t filter kmer with counts over a percentile (Value [0-1]), default: "
	<<count_cutoff<<std::endl
    <<"   -D, --diversity_cutoff:\n\t filter kmer with diversity measure of kmer size used for background probability (Value [0-1]), default: "
    <<diversity_cutoff<<std::endl
    <<"   -m, --min_count:\n\t filter kmer with counts less than this value, default: "
	<<min_count<<std::endl
    <<"   -b, --background_kmer_size:\n\t kmer size to compute kmer background probability, default: "
	<<bkmer_size<<std::endl
      <<"   -o, --file_out:\n\t filename for output,"<<std::endl
      <<"   -c, --chunk_size:\n\t sequence chunk size in kb, default: None"<<std::endl
      <<"   -a, --analysis:\n\t compute kmer statistics only"<<std::endl
  	  <<"   -v, --verbosity:\n\t verbosity level, default:"<<verbosity<<std::endl;
};
void show_parameter(SDGString filename1,SDGString filename2)
{
  std::cout<<"\nRun with parameters:\n"
	  <<"Query sequences: "<<filename1<<std::endl
	  <<"Model sequences: "<<filename2<<std::endl
      <<"   -w, --kmer:\t kmer length: "<<kmer_size<<std::endl
      <<"   -S, --step_q:\t step on query sequence: "<<step_q<<std::endl
      <<"   -k, --kmask:\t length of k-mer mask: "<<kmask<<std::endl
      <<"   -d, --kmer_dist:\t max number of kmer between two matching kmer to connect: "<<kmer_dist<<std::endl
      <<"   -s, --min_size:\t min size range to report: "<<min_size<<std::endl
      <<"   -C, --filter_cutoff:\t filter kmer with counts in the last percentile: "<<count_cutoff<<std::endl
      <<"   -D, --diversity_cutoff:\n\t filter kmer with diversity measure of kmer size used for background probability: "<<diversity_cutoff<<std::endl
      <<"   -m, --min_count:\t filter kmer with counts less than this value: "<<min_count<<std::endl
      <<"   -b, --background_kmer_size:\t kmer size to compute kmer background probability: "<<bkmer_size<<std::endl
      <<"   -o, --file_out:\t filename for output:"<<outfilename<<std::endl
      <<"   -c, --chunk_size:\t sequence chunk size in kb: "<<chunk_size_kb<<std::endl
  	  <<"   -v, --verbosity:\t verbosity level: "<<verbosity<<std::endl;
};

int main(int argc, char* argv[])
{
  try{
		std::cout<<"Beginning Hasher (version "<<VERSION<<")"<<std::endl;
		clock_t begin, end;
              double time_spent;
              begin = clock();

    if(argc==1)
      {
		help();
		exit(EXIT_SUCCESS);
      } 
    int c;
    while (1)
      {
		static struct option long_options[] =
		{
		  {"help",no_argument, 0, 'h'},
		  {"kmer",required_argument, 0, 'w'},
		  {"step_q",required_argument, 0, 'S'},
		  {"kmask",required_argument, 0, 'k'},
		  {"kmer_dist",required_argument, 0, 'd'},
		  {"min_size",required_argument, 0, 's'},
		  {"filter_cutoff",required_argument, 0, 'C'},
		  {"min_count",required_argument, 0, 'm'},
		  {"diversity_cutoff",required_argument, 0, 'D'},
		  {"background_kmer_size",required_argument, 0, 'b'},
		  {"file_out",required_argument, 0, 'o'},
		  {"chunk_size",required_argument, 0, 'c'},
		  {"analysis",no_argument, 0, 'a'},
		  {"verbosity",no_argument, 0, 'v'},
		  {0, 0, 0, 0}
		};
		/* `getopt_long' stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "hd:f:w:S:k:d:s:C:D:m:b:o:c:av:",
				 long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
		  break;

		switch (c)
		  {
		  case 'h':
			{
			  help();
			  return 0;
			}
		  case 'w':
			{
			  kmer_size=atoi(optarg);
			  break;
			}
		  case 'S':
			{
			  step_q=atoi(optarg);
			  break;
			}
		  case 'k':
			{
			  kmask=atoi(optarg);
			  break;
			}
		  case 'd':
			{
				kmer_dist=atoi(optarg);
			  break;
			}
		  case 's':
			{
				min_size=atoi(optarg);
			  break;
			}
		  case 'C':
			{
				count_cutoff=atof(optarg);
			  break;
			}
		  case 'D':
			{
			  diversity_cutoff=atof(optarg);
			  break;
			}
		  case 'm':
			{
				min_count=atoi(optarg);
			  break;
			}
		  case 'b':
			{
				bkmer_size=atoi(optarg);
				break;
			}
		  case 'o':
			{
			  outfilename=optarg;
			  break;
			}
		  case 'c':
			{
			  chunk_size_kb=atoi(optarg);
			  break;
			}
		  case 'a':
			{
			  stat_only=true;
			  break;
			}

		  case 'v':
			{
			  verbosity=atoi(optarg);
			  break;
			}
		  case '?':
		    {
				help();
				break;
			}
		    return 1;
		  default:
		  	{
		  		abort();
		  		break;
		    }
		  }
      }

    SDGString filename1,filename2;

    filename1=argv[optind];
    if(++optind < argc)
    	filename2=argv[optind];

    /* Print any remaining command line arguments (not options). */
    if (++optind < argc)
      {
		help();
		std::cout<<"non-option ARGV-elements: "<<std::endl;
		while (optind < argc)
		  std::cout<<argv[optind++]<<std::endl;
		return 1;
      }
     
    if(filename2.empty())
    {
    	repeat=true;
    	filename2=filename1;
    	if(min_count==0)
    		min_count=1;
    	std::cout<<"De novo mode! min_count="<<min_count<<std::endl;
    }

    show_parameter(filename1,filename2);

    //Check parameters
    if(count_cutoff<0 || count_cutoff>1)
    {
    	std::cout<<"count_cutoff must be in interval [0-1]! Entered value is "<<count_cutoff<<std::endl;
    	exit( EXIT_FAILURE );
    }
    if(diversity_cutoff<0 || diversity_cutoff>1)
    {
    	std::cout<<"diversity_cutoff must be in interval [0-1]! Entered value is "<<diversity_cutoff<<std::endl;
    	exit( EXIT_FAILURE );
    }


    if(stat_only)
    {
    	std::cout<<"\nCompute kmer stat only!"<<std::endl;
    	for(unsigned bw=1; bw<=bkmer_size; bw++)
    	{
    	    HashDNASeq hsrch(kmer_size, kmask, bw, kmer_dist, 0, min_size, step_q );
        	std::vector<unsigned> kmer_count((unsigned)pow(4,hsrch.getEffectiveKmerSize()),0);
        	std::list< Info_kmer > list_infokmer;
        	Info_kmer kmer_threshold;
        	unsigned nb_kmer;

        	std::cout<<"\n======Compute kmer background probability for "<<bw-1<<" Markov's chain order======"<<std::endl;

    		hsrch.kmer_analysis(filename2,kmer_size,kmask, bw, kmer_size/2, count_cutoff, diversity_cutoff, kmer_count, nb_kmer, list_infokmer, kmer_threshold);
    	}
    	std::cout<<"\nEnd HashDNASeq (version "<<VERSION<<")"<<std::endl;
    	exit( EXIT_SUCCESS );
    }

    Hasher hsrch(kmer_size, kmask, bkmer_size,kmer_dist,0, min_size,step_q);
    bool valid_idx_file=true;

	hsrch.load(filename2,kmer_size, kmask, bkmer_size,kmer_size/2 , count_cutoff, diversity_cutoff, min_count,valid_idx_file);

	std::ofstream out;
	std::stringstream out_name;
	if(outfilename!="")
		out_name<<outfilename<<".align";
	else
		out_name<<filename1<<".hasher.align";

	out.open(out_name.str());


	SDGFastaIstream in(filename1);
	if(!in)
	{
		std::cerr<<"file:"<<filename1<<" does not exist!"<<std::endl;
		return 1;
	}

	min_frag_size=min_size;
	unsigned numseq=0;
	while(in)
	  {
		SDGBioSeq s;
		if(in)
		  in>>s;
		numseq++;
		std::cout<<s.getDE()<<" len:"<<s.length()<<" read!"<<std::endl;
		if(chunk_size_kb!=0)
		{
			unsigned start=1;
			unsigned chunk_size=((chunk_size_kb*1000)/kmer_size)*kmer_size;
			unsigned nb_chunk=s.length()/chunk_size;
			for(unsigned i=1;i<nb_chunk;i++)
			{
				std::cout<<"==>chunk #"<<i<<"/"<<nb_chunk<<":"<<start<<".."<<start+chunk_size-1<<std::endl;
				hsrch.search(s,start,start+chunk_size-1,numseq,(kmer_dist+1)*kmer_size, min_frag_size, repeat, out, verbosity);

				start=start+chunk_size;
			}
			std::cout<<"==>chunk #"<<nb_chunk<<"/"<<nb_chunk<<":"<<start<<".."<<s.length()<<std::endl;
			hsrch.search(s,start,s.length(),numseq,(kmer_dist+1)*kmer_size, min_frag_size, repeat, out, verbosity);
		}else
		{
			hsrch.search(s,1,s.length(),numseq,(kmer_dist+1)*kmer_size, min_frag_size, repeat, out, verbosity);
		}
		std::cout<<"ok!\n"<<std::endl;
	  }
	out.close();




	std::cout<<"\nEnd Hasher (version "<<VERSION<<")"<<std::endl;
	end = clock();
    time_spent=(double)(end-begin)/CLOCKS_PER_SEC;
    std::cout<<"====>Total time spent****: "<<time_spent<<std::endl;

	exit( EXIT_SUCCESS );
  }
	catch( SDGException e )
	{
		std::cout<<"******Exception catched: "<<e.message<<" ******"<<std::endl;
		exit( EXIT_FAILURE );
	}
	catch(...)
	{
		std::cout<<"****** unknown exception catch !!! ******"<<std::endl;
		exit( EXIT_FAILURE );
	}
}












