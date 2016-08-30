#include <unistd.h> //getopt
#include <fstream>
#include <SDGFastaIstream.h>
#include <SDGMemBioSeq.h>

#include "Lalign.h"

int match=1,mismh=3,gap_extend=2,gap_open=5,max_subalign=1;
unsigned sim_rep=0;
bool tab_output=false;

void help(void)
{
  std::cerr<<"usage: SWalign"
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
      <<"   -s, --subalign:\n\t number of subalignments, default: "
      <<max_subalign<<std::endl
//       <<"   -r, --sim:\n\t simulation repeats for z_score\n\t (not for suboptimal alignment), default: "
//       <<sim_rep<<std::endl
      <<"   -t, --tab:\n\t tabulated output, default: FALSE"<<std::endl;
};

int main(int argc, char* argv[])
{
  try{
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
// 	  {"sim",required_argument, 0, 'r'},
// 	  {"subalign",required_argument, 0, 's'},
// 	  {"tab",required_argument, 0, 't'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;
	
// 	c = getopt_long (argc, argv, "htm:d:g:e:s:r:",
// 			 long_options, &option_index);
	
	c = getopt (argc, argv, "htm:d:g:e:s:");
	
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
	  case 't':
	    {
	      tab_output=true;
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
	  case 's':
	    {
	      max_subalign=atoi(optarg);
	      break;
	    }
	  case 'r':
	    {
	      sim_rep=atoi(optarg);
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
    SDGString filename1=argv[optind];
    SDGString filename2;
    if(++optind<argc)
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
     
    Lalign al(max_subalign);
    al.setMismatch(match,mismh);
    al.setGap(gap_open, gap_extend);

    SDGBioSeq s1,s2;
    if(filename2=="" || filename1==filename2)
      {
	SDGFastaIstream in1(filename1);
	in1>>s1;
	std::cout<<argv[1]<<": "<<s1.length()<<std::endl;
	al.setSeq(s1);
      }
    else
      {
	SDGFastaIstream in1(filename1);
	SDGFastaIstream in2(filename2);
	in1>>s1;
	std::cout<<argv[1]<<": "<<s1.length()<<std::endl;
	in2>>s2;
	std::cout<<argv[2]<<": "<<s2.length()<<std::endl;

	al.setSeq(s1,s2);
      }

  al.align();
  int score=al.getScore();
  for(int i=0;i<max_subalign;i++)
    {
      if(tab_output)
	al.write_map();
      else
	{
	  al.view();
// 	  if(sim_rep>0)
// 	    if(i==0)
// 	      std::cout<<"Z score="<<al.z_score(sim_rep)<<std::endl;
// 	    else
// 	      std::cout<<"Z score= Not available !"<<std::endl;
	}
	if(max_subalign==1)
	{
	  al.setSeq(s1,s2.complement());
	  al.align();
	  if(score<al.getScore())
	    {
	      std::cout<<"*** WARNING best on complement! ***"<<std::endl;
	      std::cout<<" -> follows alignement with complement of sequence 2"<<std::endl;
	      if(tab_output)
		al.write_map();
	      else
		al.view();
	      std::cout<<"*** WARNING best on complement! ***"<<std::endl;
	    }
	  else
	    std::cout<<"*** not best on complement! ***"<<std::endl;
	}
      if(i!=max_subalign-1) al.findNext();
    }
  }
  catch(SDGException e)
    {
      std::cerr << e.message << std::endl;
    }

  return 0;
}
