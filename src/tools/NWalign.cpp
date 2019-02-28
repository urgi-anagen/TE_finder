/*
 * \file NWalign.cpp
 */

#include <unistd.h> //getopt
#include <fstream>
#include <SDGString.h>
#include <SDGFastaIstream.h>
#include <SDGMemBioSeq.h>

#include "Galign.h"

int match=1,mismh=3,gap_extend=2,gap_open=5,gap_len=50;
unsigned sim_rep=0;
SDGString outfile="";
bool rev=true;
unsigned verbose=0;

void help(void)
{
  std::cerr<<"usage: NWalign"
      <<" [<options>] <fasta sequence> <fasta sequence>"<<std::endl
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
      <<"   -o, --out:\n\t out fasta aligned filename, default: no output"
      <<std::endl
      <<"   -D, --direct:\n\t only on direct strand, default: no"
      <<std::endl
      <<"   -v, --verbose:\n\t verbosity level, default: "
      <<verbose<<std::endl;
//       <<"   -s, --sim:\n\t simulation repeats for z_score, default: "
//       <<sim_rep<<std::endl;
};

int main(int argc, char* argv[])
{
  try
  {
    int c;
    while (1)
    {
// 	static struct option long_options[] =
// 	{
// 	  {"help",no_argument, 0, 'h'},
// 	  {"match",required_argument, 0, 'm'},
// 	  {"mismatch",required_argument, 0, 'd'},
// 	  {"gapopen",required_argument, 0, 'g'},
// 	  {"gapextend",required_argument, 0, 'e'},
// 	  {"sim",required_argument, 0, 's'},
// 	  {"lengap",required_argument, 0, 'l'},
// 	  {"direct",no_argument, 0, 'D'},
// 	  {"out",required_argument, 0, 'o'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;
	
// 	c = getopt_long (argc, argv, "hm:d:g:e:l:Ds:o:",
// 			 long_options, &option_index);
    	c = getopt(argc, argv, "hm:d:g:e:l:Do:v:");


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
    	case 'D':
	    {
	      rev=false;
	      break;
	    }
    	case 's':
	    {
	      sim_rep=atoi(optarg);
	      break;
	    }
    	case 'o':
	    {
	      outfile=optarg;
	      break;
	    }
    	case 'v':
    	{
    		verbose = atoi( optarg );
    		break;
    	}
    	case '?':
    		help();
    		return 1;
    	default:
    		abort ();
    	}
    }
    if(argc==optind || argc==optind+1)
    {
    	help();
    	return 1;
    }

    SDGString filename1 = argv[optind];
    SDGString filename2 = argv[++optind];
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

    SDGFastaIstream in1(filename1);
    SDGFastaIstream in2(filename2);

    SDGBioSeq s1;
    in1>>s1;
    if(verbose>0)
    	std::cout<<s1.getDE()<<": "<<s1.length()<<std::endl;

    SDGBioSeq s2;
    in2>>s2;
    if(verbose>0)
    	std::cout<<s2.getDE()<<": "<<s2.length()<<std::endl;

    SDGString als1,als2,nals1,nals2;

    al.setSeq(s1,s2);
    al.align();
    int score=al.getScore();
    unsigned match=al.getNumberMatch();
    unsigned mismatch=al.getNumberMismatch();
    double identity=((double)match)/(match+mismatch);
    if(verbose>0)
    	std::cout<<"score="<<score<<" identity="<<identity<<std::endl;

    if(outfile=="")
    	al.view();
    else
    {
    	al.getAlignedSeq(als1,als2);
    	nals1=al.getNameSeq1();
    	nals2=al.getNameSeq2();
    }
    //   if(sim_rep>0)
    //     std::cout<<"Z score="<<al.z_score(sim_rep)<<std::endl;

    if(rev)
    {
    	al.setSeq(s1,s2.complement());
    	al.align();
    	if(score<al.getScore())
    	{
    		score=al.getScore();
    		match=al.getNumberMatch();
    		mismatch=al.getNumberMismatch();
    		identity=((double)match)/(match+mismatch);
        if(verbose>0)
        {
        	std::cout<<"WARNING: best on complement!"<<std::endl;
        	std::cout<<"score="<<score<<" identity="<<identity<<std::endl;
        }
    		if(outfile=="")
    			al.view();
    		else
    		{
    			al.getAlignedSeq(als1,als2);
    			nals1=al.getNameSeq1();
    			nals2=al.getNameSeq2();
    		}
    		// 	  if(sim_rep>0)
    		// 	    std::cout<<"Z score="<<al.z_score(sim_rep)<<std::endl;
    	}
    	else
        if(verbose>0)
        	std::cout<<"WARNING: not best on complement!"<<std::endl;
    }

  if( outfile != "" )
  {
  	SDGBioSeq seqal1=newSDGMemBioSeq(als1);
  	seqal1.setDE(nals1);

  	SDGBioSeq seqal2=newSDGMemBioSeq(als2);
  	seqal2.setDE(nals2);

  	SDGFastaOstream out(outfile);
  	out<<seqal1<<seqal2;
  	out.close();
  }
  }
  catch(SDGException e)
  {
  	std::cerr << e.message << std::endl;
  }

  return 0;
}
