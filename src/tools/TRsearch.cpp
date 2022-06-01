#include <unistd.h> //getopt
#include <fstream>
#include <SDGFastaIstream.h>
#include <SDGMemBioSeq.h>

#include "Lalign.h"
#include "ItrAlign.h"

int match=10,mismh=16,gap_extend=32,gap_open=32;
double id_thres=0.8, ext_len=0.5;
unsigned max_align_len=7000, nb_sub=3, min_size=10, min_len=10;
    
void help(void)
{
  std::cerr<<"usage: TRsearch"
      <<" [<options>] <fasta sequence database>"<<std::endl
      <<" options:"<<std::endl
      <<"   -h, --help:\n\t this help"<<std::endl
      <<"   -i, --identity:\n\t min identity threshold, default: "
      <<id_thres<<std::endl
      <<"   -l, --minlen:\n\t minimum length repeat, default: "
      <<min_len<<std::endl
      <<"   -m, --match:\n\t match bonus (>0), default: "
      <<match<<std::endl
      <<"   -d, --mismatch:\n\t mismatch penalty (>0), default: "
      <<mismh<<std::endl
      <<"   -g, --gapopen:\n\t gap open penalty (>0), default: "
      <<gap_open<<std::endl
      <<"   -e, --gapextend:\n\t gap extend penalty (>0), default: "
      <<gap_extend<<std::endl
      <<"   -x, --extlen:\n\t length from extremities in %, default: "
      <<ext_len<<std::endl
      <<"   -s, --minsize:\n\t min element size to report, default: "
      <<min_size<<std::endl;
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
// 	  {"id_thres",required_argument, 0, 'i'},
// 	  {"minlen",required_argument, 0, 'l'},
// 	  {"match",required_argument, 0, 'm'},
// 	  {"mismatch",required_argument, 0, 'd'},
// 	  {"gapopen",required_argument, 0, 'g'},
// 	  {"gapextend",required_argument, 0, 'e'},
// 	  {"minsize",required_argument, 0, 's'},
// 	  {"extlen",required_argument, 0, 'x'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;
	
// 	c = getopt_long (argc, argv, "hi:l:m:d:g:e:s:x:",
// 			 long_options, &option_index);

	c = getopt(argc, argv, "hi:l:m:d:g:e:s:x:");

	
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
	  case 'i':
	    {
	      id_thres=atof(optarg);
	      break;
	    }
	  case 'l':
	    {
	      min_len=atoi(optarg);
	      break;
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
	  case 's':
	    {
	      min_size=atoi(optarg);
	      break;
	    }
	  case 'x':
	    {
	      ext_len=atof(optarg);
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
     
    Lalign al(nb_sub);
    al.setMismatch(match,mismh);
    al.setGap(gap_open, gap_extend);
    
    ItrAlign altir(min_len,id_thres);
    altir.setMismatch(match,mismh);
    altir.setGap(gap_open, gap_extend);

    SDGFastaIstream in(filename);
    SDGBioSeq s;

    SDGString outfilename=filename+".TR.set";
    std::ofstream fout(outfilename);

    unsigned count=1;
    while(in)
      {
	in>>s;
	std::cout<<s.getDE()<<" "<<s.length()<<".."<<std::flush;
	unsigned len=(unsigned)floor(ext_len*s.length());
	unsigned slen=s.length();

	SDGBioSeq s1,s2,s2c;

	unsigned align_len=len<max_align_len?len:max_align_len;
	std::cout<<"Alignment length="<<align_len<<std::endl;
	if(max_align_len==align_len)
	  std::cout<<"WARNING: max alignment length reached !!"<<std::endl;
	s1=s.subseq(0,align_len);
	s2=s.subseq(slen-align_len,slen);

	//run_test_search_wSW LTR method1 as best local alignments between start and end
	//of the sequence
	al.setSeq(s1,s2);	  
	for(unsigned n=0;n<nb_sub;n++)
	  {
	    if(n==0) al.align();
	    else al.findNext();
	    if(al.getScore()==0) break;
	    unsigned r1s=al.getStartSeq1();
	    unsigned r1e=al.getEndSeq1();
	    unsigned r2s=al.getStartSeq2()+slen-align_len;
	    unsigned r2e=al.getEndSeq2()+slen-align_len;

	    unsigned lenRep=r1e-r1s+1;
	    if( lenRep>=100 && al.getIdentity()>id_thres )
	      {
		unsigned lengthTE=r2e-r1s+1;
		if(lengthTE < min_size)
		  continue;

		if(r1s==r2s && r1e==r2e &&  al.getIdentity()==1 
		   && r1s==slen-align_len+1 && r1e==align_len)
		  continue;

		SDGString name="";
		if(r1s<5 && slen-r2e<5)
		  name="term";
		name+="LTR|len="
		  +SDGString(lenRep)
		  +"|id="+SDGString(al.getIdentity())
		  +"|lenTE="+SDGString(lengthTE);
		std::cout<<name<<": "<<r1s<<".."<<r1e
			 <<"   "<<r2s<<".."<<r2e<<std::endl;
		SDGString str="Rep"+SDGString(count)+name;
		fout<<count<<"\t"<<str<<"\t"<<s.getDE()
		    <<"\t"<<r1s<<"\t"<<r1e<<std::endl;
		fout<<count<<"\t"<<str<<"\t"<<s.getDE()
		    <<"\t"<<r2s<<"\t"<<r2e<<std::endl;
		count++;
	      }
	  }

	//run_test_search_wSW ITR method1 as best local alignments between start and end
	//of the sequence
	s2c=s2.complement();
	al.setSeq(s1,s2c);	  
	for(unsigned n=0;n<nb_sub;n++)
	  {
	    if(n==0) al.align();
	    else al.findNext();
	    if(al.getScore()==0) break;
	    unsigned r1s=al.getStartSeq1();
	    unsigned r1e=al.getEndSeq1();
	    unsigned r2s=align_len-al.getStartSeq2()+slen-align_len+1;
	    unsigned r2e=align_len-al.getEndSeq2()+slen-align_len+1;
	    
	    unsigned lenRep=r1e-r1s+1;
	    if( lenRep>=min_len && al.getIdentity()>id_thres )
	      {
		unsigned lengthTE=r2s-r1s+1;
		if(lengthTE < min_size)
		  continue;

		SDGString name="non-termTIR|len="
		  +SDGString(lenRep)
		  +"|id="+SDGString(al.getIdentity())
		  +"|lenTE="+SDGString(lengthTE);
		std::cout<<name<<": "<<r1s<<".."<<r1e
			 <<"   "<<r2s<<".."<<r2e<<std::endl;
		SDGString str="Rep"+SDGString(count)+name;
		fout<<count<<"\t"<<str<<"\t"<<s.getDE()
		    <<"\t"<<r1e<<"\t"<<r1s<<std::endl;
		fout<<count<<"\t"<<str<<"\t"<<s.getDE()
		    <<"\t"<<r2s<<"\t"<<r2e<<std::endl;
		count++;
	      }
	  }
	    
	//run_test_search_wSW TIR methode 2 at the boundaries
	altir.alignItr(s);
	if(altir.getScore()!=0)
	  {
	    unsigned r1s=altir.getStartSeq1();
	    unsigned r1e=altir.getEndSeq1();
	    unsigned r2s=altir.getEndSeq2();
	    unsigned r2e=altir.getStartSeq2();
	    
	    unsigned lenRep=r1e-r1s+1;
	    if( lenRep>=min_len && altir.getIdentity()>id_thres )
	      {
		unsigned lengthTE=r2s-r1s+1;
		if(lengthTE < min_size)
		  continue;
		SDGString name="termTIR|len="
		  +SDGString(lenRep)
		  +"|id="+SDGString(altir.getIdentity())
		  +"|lenTE="+SDGString(lengthTE);
		std::cout<<name<<": "<<r1s<<".."<<r1e
			 <<"   "<<r2s<<".."<<r2e<<std::endl;
		
		SDGString str="Rep"+SDGString(count)+name;
		fout<<count<<"\t"<<str<<"\t"<<s.getDE()
		    <<"\t"<<r1e<<"\t"<<r1s<<std::endl;
		fout<<count<<"\t"<<str<<"\t"<<s.getDE()
		    <<"\t"<<r2s<<"\t"<<r2e<<std::endl;
		count++;
	      }
	  }
      }

    in.close();
  }
  catch(SDGException e)
    {
      std::cerr << e.message << std::endl;
    }
  return 0;
}






