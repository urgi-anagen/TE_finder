#include <unistd.h> //getopt
#include <fstream>
#include <vector>
#include <algorithm>

#include "CountStat.h"
#include "Quantile.h"
#include <SDGString.h>
#include <SDGFastaIstream.h>
#include <SDGFastaOstream.h>
#include <SDGMemBioSeq.h>

#include "Galign.h"

int match=10,mismh=8,gap_extend=4,gap_open=16,gap_len=10;
void help(void)
{
  std::cerr<<"usage: refalign"
      <<" [<options>] <reference fasta sequence> <fasta database sequence>"<<std::endl
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

    if(argc==optind+1 || argc==1)
      {
	help();
	return 1;
      }

    SDGString filename1=argv[optind];
    SDGString filename2=argv[++optind];
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

  if(filename_out=="")
      filename_out=filename2.afterlast("/");

  SDGString aligner_filename=filename_out;
  aligner_filename+=".aligner";
  std::ofstream fout1(aligner_filename);

  SDGString fileout1(filename_out);
  fileout1+=".refalign.stat";
  std::ofstream fout(fileout1);

  SDGString fileout2(filename_out);
  fileout2+=".oriented";
  SDGFastaOstream seqout(fileout2);


  Quantile quant_identity,quant_gaps,quant_lengths,quant_scores;
  unsigned nbcomp=0;
  std::vector<int> cec_id;

  SDGBioSeq refSeq;
  SDGFastaIstream in1(filename1);
  in1>>refSeq;
  std::cout<<refSeq.getDE()<<std::endl;
  in1.close();

  SDGFastaIstream in2(filename2);
  int i=0;
  fout<<"seqname\tlength\tidentity\tgaps\tscore"<<std::endl;
  while(in2)
    {
      SDGBioSeq seq;

      std::cout<<"load sequence #"<<++i;
      in2>>seq;
      std::cout<<" "<<seq.getDE()<<std::flush;
      std::cout<<" length="<<seq.length()<<std::flush;
      fout<<seq.getDE()<<"\t"<<seq.length()<<std::flush;
     
      al.setSeq(refSeq,seq);	  
      al.align();
      SDGString s1,s2;
      al.getAlignedSeq(s1,s2);
      int start=0;
      while(s1.charAt(start)=='-')  start++;
      int end=s1.length()-1;
      while(s1.charAt(end)=='-')  end--;
      int score=al.getScore();
      unsigned length=al.getAlignmentLength();
      unsigned nbmatch=al.getNumberMatch();
      unsigned nbmismatch=al.getNumberMismatch();
      al.setSeq(refSeq,seq.complement());
      al.align();
      int scorecomp=al.getScore();
      if(scorecomp==0 && score==0)
	{
	  std::cout<<" no alignment!!"<<std::endl;
	  continue;
	}
      if(score<scorecomp)
	{
	  std::cout<<" best on complement -> complement!";
	  fout<<"comp";
	  seqout<<seq.complement();
	  score=al.getScore();
	  length=al.getAlignmentLength();
	  nbmatch=al.getNumberMatch();
	  nbmismatch=al.getNumberMismatch();
	  al.getAlignedSeq(s1,s2);
	  start=0;
	  while(s1.charAt(start)=='-')  start++;
	  end=s1.length()-1;
	  while(s1.charAt(end)=='-')  end--;
	  nbcomp++;
	}
      else
	{
	  seqout<<seq;
	  std::cout<<" OK!";
	}

      double identity=(double)nbmatch/(nbmatch+nbmismatch);
      double gaps=(double)(length-nbmatch-nbmismatch)/length;
      std::cout<<" identity="<<identity
	  <<" gaps="<<gaps
	  <<" Score="<<score
	  <<std::endl<<std::flush;

      fout<<"\t"<<identity
	  <<"\t"<<gaps
	  <<"\t"<<score
	  <<std::endl<<std::flush;

      fout1<<s1.substr(start,end-start+1)<<"\t"<<s2.substr(start,end-start+1)<<"\t"<<seq.getDE()<<std::endl;
      quant_identity.add(identity);
      quant_gaps.add(gaps);
      quant_lengths.add(length);
      quant_scores.add(score);
      in2.peek();
    }

  
  std::cout<<"\n           \t"
      <<"min\t"
      <<"Q25\t"
      <<"median\t"
      <<"mean\t"
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
      <<"mean\t"
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




