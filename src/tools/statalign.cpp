#include <unistd.h> //getopt
#include <fstream>
#include <algorithm>

#include "Quantile.h"
#include <SDGString.h>
#include <SDGFastaIstream.h>
#include <SDGFastaOstream.h>
#include <SDGMemBioSeq.h>
#include "SDGBioSeqDB.h"

#include "Galign.h"

int match=10,mismh=8,gap_extend=4,gap_open=16,gap_len=10;

void help(void)
{
  std::cerr<<"usage: statalign"
      <<" [<options>] <fasta database sequence>"<<std::endl
      <<" options:"<<std::endl
      <<"   -h, --help:\n\t this help"<<std::endl
      <<"   -m, --match:\n\t match bonus (>0), default: "
      <<match<<std::endl
      <<"   -d, --mismatch:\n\t mismatch penalty (>0), default: "
      <<mismh<<std::endl
      <<"   -g, --gapopen:\n\t gap open penalty (>0), default: "
      <<gap_open<<std::endl
      <<"   -e, --gapextend:\n\t gap extend penalty (>0), default: "
      <<gap_extend<<std::endl
      <<"   -l, --lengap:\n\t max length of penalized gap, default: "
      <<gap_len<<std::endl
      <<"   -o, --outfile:\n\t prefix of output files, default: <fasta database sequence>"
      <<std::endl;
};

int main(int argc, char* argv[])
{
   try{

    int c;
    SDGString filename_out="";
    while (1)
      {
// 	static struct option long_options[] =
// 	{
// 	  {"help",no_argument, 0, 'h'},
// 	  {"match",required_argument, 0, 'm'},
// 	  {"mismatch",required_argument, 0, 'd'},
// 	  {"gapopen",required_argument, 0, 'g'},
// 	  {"gapextend",required_argument, 0, 'e'},
// 	  {"lengap",required_argument, 0, 'l'},
// 	  {"outfile",required_argument, 0, 'o'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;
	
// 	c = getopt_long (argc, argv, "hm:d:g:e:l:o:",
// 			 long_options, &option_index);

	c = getopt(argc, argv, "hm:d:g:e:l:o:");
	
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
	  case 'm':
	    {
	      match=atoi(optarg);
	      break;
	    }
	  case 'd':
	    {
	      mismh=atoi(optarg);
	      break;
	    }
	  case 'g':
	    {
	      gap_open=atoi(optarg);
	      break;
	    }
	  case 'e':
	    {
	      gap_extend=atoi(optarg);
	      break;
	    }
	  case 'l':
	    {
	      gap_len=atoi(optarg);
	      break;
	    }
	  case 'o':
	    {
	      filename_out=optarg;
	      break;
	    }
	  case '?':
	    help();
	    return 1;
	  default:
	    abort ();
	  }
      }
   
    if(argc==optind)
      {
	help();
	return 1;
      }

    SDGString filename=argv[optind];
    /* Print any remaining command line arguments (not options). */
    if (++optind < argc)
      {
	help();
	std::cout<<"non-option ARGV-elements: "<<std::endl;
	while (optind < argc)
	  std::cout<<argv[optind++]<<std::endl;
	return 1;
      }
     
       
  Galign al;
  al.setMismatch(match,mismh);
  al.setGap(gap_open, gap_extend, gap_len);

  Quantile quant_identity;
  Quantile quant_gaps;
  Quantile quant_scores;
  Quantile quant_lengths;

  if(filename_out=="")
      filename_out=filename.afterlast("/");
  filename_out+=".statalign.out";
  std::ofstream fout(filename_out);
  int x=0,y=0;

  SDGBioSeqDB db(filename);

  std::vector<double> vec_id;

  fout<<"seqname1\tseqname2\tscore\tlength\tidentity\tgaps"<<std::endl;
  for( SDGBioSeqDB::iterator i=db.begin(); i!=db.end();i++)
    {
      x++;
      y=x;
      for( SDGBioSeqDB::iterator j=i+1; j!=db.end();j++)
	{
	  y++;
	  al.setSeq(*i,*j);	  
	  std::cout<<" aligning ..."
		   <<i->getDE().before(" ")
		   <<" with "
		   <<j->getDE().before(" ")
		   <<std::flush;
	  fout<<i->getDE().before(" ")
		   <<"\t"
		   <<j->getDE().before(" ")
		   <<"\t";
	  al.align();
	  int score=al.getScore();
	  if(score!=0)
	    {
	      unsigned length=al.getAlignmentLength();
	      unsigned nbmatch=al.getNumberMatch();
	      unsigned nbmismatch=al.getNumberMismatch();

	      std::cout<<" score="<<score
		  <<" length="<<length
		  <<" identity="<<(double)nbmatch/(nbmatch+nbmismatch)
		  <<" gaps="<<(double)(length-nbmatch-nbmismatch)/length
		  <<std::endl<<std::flush;

	      fout<<score
		  <<"\t"<<length
		  <<"\t"<<(double)nbmatch/(nbmatch+nbmismatch)
		  <<"\t"<<(double)(length-nbmatch-nbmismatch)/length
		  <<std::endl<<std::flush;

	      quant_scores.add(score);
	      quant_lengths.add(length);
	      quant_identity.add((double)nbmatch/(nbmatch+nbmismatch));
	      quant_gaps.add((double)(length-nbmatch-nbmismatch)/length);
	    }
	  else //no alignment
	    {
	      std::cout<<" score="<<score
		  <<" length="<<0
		  <<" identity="<<(double)0.0
		  <<" gaps="<<(double)1.0
		  <<std::endl<<std::flush;

	      fout<<score
		  <<"\t"<<0
		  <<"\t"<<(double)0.0
		  <<"\t"<<(double)1.0
		  <<std::endl<<std::flush;

	    }
	}
    }

  std::cout<<"\n           \t"
      <<"min\t"
      <<"Q25\t"
      <<"median\t"
      <<"Q75\t"
      <<"max"<<std::endl;
  
  std::cout<<"Stat Identity\t"
      <<quant_identity.getMin()<<"\t"
      <<quant_identity.quantile(0.25)<<"\t"
      <<quant_identity.quantile(0.5)<<"\t"
      <<quant_identity.quantile(0.75)<<"\t"
	   <<quant_identity.getMax()<<std::endl;

  std::cout<<"Stat Gaps\t"
      <<quant_gaps.getMin()<<"\t"
      <<quant_gaps.quantile(0.25)<<"\t"
      <<quant_gaps.quantile(0.5)<<"\t"
      <<quant_gaps.quantile(0.75)<<"\t"
      <<quant_gaps.getMax()<<std::endl;

  std::cout<<"Stat Lengths\t"
      <<quant_lengths.getMin()<<"\t"
      <<quant_lengths.quantile(0.25)<<"\t"
      <<quant_lengths.quantile(0.5)<<"\t"
      <<quant_lengths.quantile(0.75)<<"\t"
      <<quant_lengths.getMax()<<std::endl;

  std::cout<<"Stat Scores\t"
      <<quant_scores.getMin()<<"\t"
      <<quant_scores.quantile(0.25)<<"\t"
      <<quant_scores.quantile(0.5)<<"\t"
      <<quant_scores.quantile(0.75)<<"\t"
      <<quant_scores.getMax()<<std::endl;



  fout<<"\n           \t"
      <<"min\t"
      <<"Q25\t"
      <<"median\t"
      <<"Q75\t"
      <<"max"<<std::endl;
  
  fout<<"Stat Identity\t"
      <<quant_identity.getMin()<<"\t"
      <<quant_identity.quantile(0.25)<<"\t"
      <<quant_identity.quantile(0.5)<<"\t"
      <<quant_identity.quantile(0.75)<<"\t"
	   <<quant_identity.getMax()<<std::endl;

  fout<<"Stat Gaps\t"
      <<quant_gaps.getMin()<<"\t"
      <<quant_gaps.quantile(0.25)<<"\t"
      <<quant_gaps.quantile(0.5)<<"\t"
      <<quant_gaps.quantile(0.75)<<"\t"
      <<quant_gaps.getMax()<<std::endl;
  
  fout<<"Stat Lengths\t"
      <<quant_lengths.getMin()<<"\t"
      <<quant_lengths.quantile(0.25)<<"\t"
      <<quant_lengths.quantile(0.5)<<"\t"
      <<quant_lengths.quantile(0.75)<<"\t"
      <<quant_lengths.getMax()<<std::endl;

  fout<<"Stat Scores\t"
      <<quant_scores.getMin()<<"\t"
      <<quant_scores.quantile(0.25)<<"\t"
      <<quant_scores.quantile(0.5)<<"\t"
      <<quant_scores.quantile(0.75)<<"\t"
      <<quant_scores.getMax()<<std::endl;

  }
  catch(SDGException e)
    {
      std::cerr << e.message << std::endl;
    }

  return 0;
}




