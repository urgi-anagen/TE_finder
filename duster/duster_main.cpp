#include <getopt.h>
#include <fstream>
#include <cstdlib>
#include <list>

#include "SDGString.h"
#include "SDGFastaIstream.h"
#include "SDGFastaOstream.h"
#include "SDGMemBioSeq.h"

#include "Duster.h"


unsigned kmer_size=16, step_q=1, bkmer_size=1, kmer_dist=5, frag_connect_dist=100, min_size=20, chunk_size_kb=0, nb_iter=1, min_count=0, kmask=2;
double count_cutoff=1.0, diversity_cutoff=0.0;
bool repeat=false, stat_only=false;

SDGString outfilename="";

    
void help(void)
{
  std::cerr<<"usage: duster"
      <<" [<options>] <fasta query sequence> [<fasta sequence model>]. Note: without sequence model, duster will perform its own repeat search."<<std::endl
      <<" options:"<<std::endl
      <<"   -h, --help:\n\t this help"<<std::endl
      <<"   -w, --kmer:\n\t kmer length (<16), default: "<<kmer_size<<std::endl
      <<"   -S, --step_q:\n\t step on query sequence, default: "<<step_q<<std::endl
      <<"   -k, --kmask:\n\t length of k-mer mask, default: "<<kmask<<std::endl
      <<"   -d, --kmer_dist:\n\t max number of kmer between two matching kmer to connect, default: "
	<<kmer_dist<<std::endl
    <<"   -f, --frag_connect_dist:\n\t max distance between two fragments to connect, default: "
	<<frag_connect_dist<<std::endl
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
      <<"   -c, --chunk_size:\n\t sequence chunk size in kb,"<<" default: None"<<std::endl
      <<"   -n, --nb_iter:\n\t number of iteration. A value of 0 make iteration stop if coverage variation is less than 1%."<<" default:"<<nb_iter<<std::endl
      <<"   -a, --analysis:\n\t compute kmer statistics only"<<std::endl;
};
void show_parameter(SDGString filename1,SDGString filename2)
{
  std::cout<<"\nRun with parameters:\n"
	  <<"Query sequences: "<<filename1<<std::endl
	  <<"Model sequences: "<<filename2<<std::endl
      <<"   -w, --kmer:\t kmer length (<16): "<<kmer_size<<std::endl
      <<"   -S, --step_q:\t step on query sequence: "<<step_q<<std::endl
      <<"   -k, --kmask:\t length of k-mer mask: "<<kmask<<std::endl
      <<"   -d, --kmer_dist:\t max number of kmer between two matching kmer to connect: "<<kmer_dist<<std::endl
      <<"   -f, --frag_connect_dist:\n\t max distance between two fragments to connect: "
   	<<frag_connect_dist<<std::endl
      <<"   -s, --min_size:\t min size range to report: "<<min_size<<std::endl
      <<"   -C, --filter_cutoff:\t filter kmer with counts in the last percentile: "<<count_cutoff<<std::endl
      <<"   -D, --diversity_cutoff:\n\t filter kmer with diversity measure of kmer size used for background probability: "<<diversity_cutoff<<std::endl
      <<"   -m, --min_count:\t filter kmer with counts less than this value: "<<min_count<<std::endl
      <<"   -b, --background_kmer_size:\t kmer size to compute kmer background probability: "<<bkmer_size<<std::endl
      <<"   -o, --file_out:\t filename for output:"<<outfilename<<std::endl
      <<"   -c, --chunk_size:\t sequence chunk size in kb: "<<chunk_size_kb<<std::endl
      <<"   -n, --nb_iter:\t number of iteration: "<<nb_iter<<std::endl;
};

int main(int argc, char* argv[])
{
  try{
		std::cout<<"Beginning Duster (version "<<VERSION<<")"<<std::endl;
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
		  {"frag_connect_dist",required_argument, 0, 'f'},
		  {"min_size",required_argument, 0, 's'},
		  {"filter_cutoff",required_argument, 0, 'C'},
		  {"min_count",required_argument, 0, 'm'},
		  {"diversity_cutoff",required_argument, 0, 'D'},
		  {"background_kmer_size",required_argument, 0, 'b'},
		  {"file_out",required_argument, 0, 'o'},
		  {"chunk_size",required_argument, 0, 'c'},
		  {"nb_iter",required_argument, 0, 'n'},
		  {"analysis",no_argument, 0, 'a'},
		  {0, 0, 0, 0}
		};
		/* `getopt_long' stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "hd:f:w:S:k:s:C:D:m:b:o:c:n:a",
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
		  case 'f':
			{
				frag_connect_dist=atoi(optarg);
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
		  case 'n':
			{
			  nb_iter=atoi(optarg);
			  break;
			}
		  case 'a':
			{
			  stat_only=true;
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
        	std::vector<unsigned> kmer_count((unsigned)pow(4,kmer_size-kmask),0);
        	std::list< Info_kmer > list_infokmer;
        	Info_kmer kmer_threshold;
        	unsigned nb_kmer;

        	std::cout<<"\n======Compute kmer background probability for "<<bw-1<<" Markov's chain order======"<<std::endl;
    	    Duster hsrch(kmer_size,kmask,bw,kmer_dist,frag_connect_dist, min_size, step_q );
    		hsrch.kmer_analysis(filename2,kmer_size,kmask, bw, kmer_size/2, count_cutoff, diversity_cutoff, kmer_count, nb_kmer, list_infokmer, kmer_threshold);
    	}
    	std::cout<<"\nEnd Duster (version "<<VERSION<<")"<<std::endl;
    	exit( EXIT_SUCCESS );
    }

    Duster hsrch(kmer_size, kmask, bkmer_size,kmer_dist,frag_connect_dist, min_size,step_q);
    bool valid_idx_file=true;
    bool first_iter=true;
	double prev_genome_perc_coverage=0.0;
    for(unsigned iter=1; iter<=nb_iter || nb_iter==0;iter++)
    {
		hsrch.load(filename2,kmer_size, kmask, bkmer_size,kmer_size/2 , count_cutoff, diversity_cutoff, min_count,valid_idx_file, first_iter);

		std::ofstream out;
		std::stringstream out_name;
		if(outfilename!="")
			out_name<<outfilename<<"."<<iter<<".bed";
		else
			out_name<<filename1<<"."<<iter<<".duster.bed";

		out.open(out_name.str());

		SDGFastaOstream seqout;
		std::stringstream seqout_name;
		if(outfilename!="")
			seqout_name<<outfilename<<"."<<iter<<".fa";
		else
			seqout_name<<filename1<<"."<<iter<<".duster.bed.fa";

		seqout.open(seqout_name.str());
	
		SDGFastaIstream in(filename1);
		if(!in)
		{
			std::cerr<<"file:"<<filename1<<" does not exist!"<<std::endl;
			return 1;
		}
		unsigned genome_size=0;
		unsigned genome_coverage=0;
		unsigned numseq=0;
		while(in)
		  {
			SDGBioSeq s;
			if(in)
			  in>>s;
			numseq++;
			genome_size+=s.length();
			std::cout<<s.getDE()<<" len:"<<s.length()<<" read!"<<std::endl;
			SDGBioSeq comp_s=s.complement();
			std::vector< std::pair<unsigned,unsigned> > frag,fmerged;
			if(chunk_size_kb!=0)
			{
				unsigned start=1;
				unsigned chunk_size=((chunk_size_kb*1000)/kmer_size)*kmer_size;
				unsigned nb_chunk=s.length()/chunk_size;
				for(unsigned i=1;i<nb_chunk;i++)
				{
					std::cout<<"==>chunk #"<<i<<"/"<<nb_chunk<<":"<<start<<".."<<start+chunk_size-1<<std::endl;
					std::cout<<"---direct strand---"<<std::endl;
					hsrch.search(s,start,start+chunk_size-1,numseq,repeat,frag);
					std::cout<<"---reverse strand---"<<std::endl;
					hsrch.search(comp_s,start,start+chunk_size-1,numseq,repeat,frag);
					start=start+chunk_size;
				}
				std::cout<<"==>chunk #"<<nb_chunk<<"/"<<nb_chunk<<":"<<start<<".."<<s.length()<<std::endl;
				std::cout<<"---direct strand---"<<std::endl;
				hsrch.search(s,start,s.length(),numseq,repeat,frag);
				std::cout<<"---reverse strand---"<<std::endl;
				hsrch.search(comp_s,start,s.length(),numseq,repeat,frag);

				hsrch.fragMerge(frag,(kmer_dist+1)*kmer_size,fmerged);
			}else
			{
				std::cout<<"---direct strand---"<<std::endl;
				hsrch.search(s,1,s.length(),numseq,repeat,frag);
				std::cout<<"---reverse strand---"<<std::endl;
				hsrch.search(comp_s,1,s.length(),numseq,repeat,frag);
				hsrch.fragMerge(frag,(kmer_dist+1)*kmer_size,fmerged);
			}
			genome_coverage+=hsrch.compute_coverage(fmerged);
			hsrch.writeBED(s.getDE(),fmerged,out);
			hsrch.get_sequences(fmerged,s,seqout);
		  }
		seqout.close();
		std::cout<<"Coverage="<<genome_coverage<<" ("<<(float)genome_coverage/genome_size<<")"
				<<" coverage % difference="<<fabs(((float)genome_coverage/genome_size)-prev_genome_perc_coverage)<<std::endl;
		if(genome_coverage==0) break;
		if(fabs(((float)genome_coverage/genome_size)-prev_genome_perc_coverage)<0.01 && nb_iter==0 && iter>1) break;
		prev_genome_perc_coverage=(float)genome_coverage/genome_size;
		filename2=seqout_name.str();
		//count_cutoff=1.0;
		min_count=0;
		repeat=false;
		first_iter=false;
    }
	std::cout<<"\nEnd Duster (version "<<VERSION<<")"<<std::endl;
	end = clock();
    time_spent=(double)(end-begin)/CLOCKS_PER_SEC;
    std::cout<<"====>Total time spent****: "<<time_spent<<std::endl;

	exit( EXIT_SUCCESS );
  }
	catch( SDGException* e )
	{
		std::cout<<"******Exception catched: "<<e->message<<" ******"<<std::endl;
		exit( EXIT_FAILURE );
	}
	catch(...)
	{
		std::cout<<"****** unknown exception catch !!! ******"<<std::endl;
		exit( EXIT_FAILURE );
	}
}












