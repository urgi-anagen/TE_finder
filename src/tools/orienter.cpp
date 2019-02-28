#include <stdio.h>
#include <unistd.h> //getopt
#include <string>
#include <fstream>
#include <algorithm>

#include <SDGFastaIstream.h>
#include <SDGFastaOstream.h>
#include <SDGMemBioSeq.h>
#include "SDGBioSeqDB.h"
#include <Alea.h>

#include "HashAlignClone.h"

int match=1, mismh=3, gap_extend=2, gap_open=5, len_thres=30, keepseq=false, verbose=0;

void help( void )
{
	std::cerr<<"usage: orienter"
	<<" [<options>] <fasta bank name>"<<std::endl
	<<"options:"<<std::endl
	<<"   -h, --help:\n\t this help"<<std::endl
	<<"   -m, --match:\n\t match bonus (>0), default: "
	<<match<<std::endl
	<<"   -d, --mismatch:\n\t mismatch penalty (>0), default: "
	<<mismh<<std::endl
	<<"   -g, --gapopen:\n\t gap open penalty (>0), default: "
	<<gap_open<<std::endl
	<<"   -e, --gapextend:\n\t gap extend penalty (>0), default: "
	<<gap_extend<<std::endl
	<<"   -l, --lenthres:\n\t significant alignment minimum len , default: "
	<<len_thres<<std::endl
	<<"   -k, --keepseq:\n\t return all sequences even if non-orientable"
	<<std::endl
	<<"   -v, --verbose:\n\t verbosity level, default: "
	<<verbose<<std::endl;
};

int main( int argc, char* argv[] )
{
	try{

		int c;
		while( 1 )
		{
// 	static struct option long_options[] =
// 	{
// 	  {"help",no_argument, 0, 'h'},
// 	  {"match",required_argument, 0, 'm'},
// 	  {"mismatch",required_argument, 0, 'd'},
// 	  {"gapopen",required_argument, 0, 'g'},
// 	  {"gapextend",required_argument, 0, 'e'},
// 	  {"lenthres",required_argument, 0, 'l'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;

// 	c = getopt_long (argc, argv, "hm:d:g:e:l:",
// 			 long_options, &option_index);

			c = getopt(argc, argv, "hm:d:g:e:l:kv:");
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
				len_thres=atoi(optarg);
				break;
			}
			case 'k':
			{
				keepseq=true;
				break;
			}
			case 'v':
			{
				verbose=atoi(optarg);
				break;
			}
			case '?':
				help();
				return 1;
			default:
				abort ();
			}
		}

		if( argc == optind )
		{
			help();
			return 1;
		}

		std::cout<<"START orienter"<<std::endl;

		SDGString filename = argv[optind];
		/* Print any remaining command line arguments (not options). */
		if (++optind < argc)
		{
			help();
			std::cout<<"non-option ARGV-elements: "<<std::endl;
			while (optind < argc)
				std::cout<<argv[optind++]<<std::endl;
			return 1;
		}

		HashAlignClone al;

		int x=0,y=0;

		SDGBioSeqDB db( filename );
		SDGFastaOstream fout( filename+".oriented" );
		int lendb=db.size();
		if(lendb==1)
		{
			std::cout<<" *** only one sequence, no job to do !"<<std::endl;
			fout<<db[0];
			fout.close();
			return 0;
		}

		std::ofstream logfile( "orienter_error.log", std::ios::app );

		std::vector<std::vector<int> > matches;
		std::vector<int> vec_strand;

		for( int i=0; i<lendb; i++ )
		{
			std::vector<int> seq_match( lendb );
			matches.push_back( seq_match );
		}

		std::cout<<"aligning input sequences..."<<std::endl;

		for( SDGBioSeqDB::iterator i=db.begin(); i!=db.end();i++ )
		{
			vec_strand.push_back(+1);
			matches[x][x]=+1;
			y=x;
			for( SDGBioSeqDB::iterator j=i+1; j!=db.end();j++)
			{
				y++;
				if( verbose > 0 )
					std::cout<<"aligning "<<x+1<<" with "<<y+1<<std::flush;
				
				
				al.setSeq(*i,*j);
				al.search();
				
				al.fragAlign(match,mismh,gap_open, gap_extend,0,0);
				
				int score=al.score();
				int len_direct=al.getMatches().back().getRangeQ().getLength();
				if( verbose > 0 )
					std::cout<<" "<<score;

				al.setSeq(*i,j->complement());
				al.search();
				al.fragAlign(match,mismh,gap_open, gap_extend,0,0);

				int score_comp=al.score();
				int len_comp=al.getMatches().back().getRangeQ().getLength();
				if( verbose > 0 )
					std::cout<<" "<<score_comp;

				if( score >= score_comp )
				{
					if( len_direct < len_thres )
					{
						matches[x][y]=0;
						matches[y][x]=0;
					}
					else
					{
						matches[x][y]=score;
						matches[y][x]=score;
					}
					if( verbose > 0 )
					{
						std::cout<<" score="<<score
						<<" strand=+"
						<<" len="<<len_direct<<std::endl;
					}
				}
				else
				{
					if( len_comp < len_thres )
					{
						matches[x][y]=0;
						matches[y][x]=0;
					}
					else
					{
						matches[x][y]=-1*score_comp;
						matches[y][x]=-1*score_comp;
					}
					if( verbose > 0 )
					{
						std::cout<<" score="<<score_comp
						<<" strand=-"
						<<" len="<<len_comp<<std::endl;
					}
				}
			}
			x++;
		}

		std::vector<int> count_comp;
		count_comp.resize( lendb );

		if( verbose > 1 )
		{
			std::cout<<std::endl;
			std::cout<<std::setw(6)<<"\\";
			for(int i=0;i<lendb;i++)
				std::cout<<std::setw(6)<<i+1;
			std::cout<<std::endl;
			for(int i=0;i<lendb;i++)
			{
				std::cout<<std::setw(6)<<i+1;
				for(int j=0;j<lendb;j++)
				{
					std::vector<int> seq_match(lendb);
					std::cout<<std::setw(6)<<matches[i][j];
				}
				std::cout<<std::endl;
			}
		}

		RandomIndex candidats;
		double max=0.0,count_try=0;
		while(1)
		{
			max=0;
			int numseq=-1;
			for(int i=0;i<lendb;i++)
			{
				int count=0,nb=0;
				for(int j=0;j<lendb;j++)
				{
					if(matches[i][j]<0) count++;
					if(matches[i][j]!=0 && matches[i][j]!=1) nb++;
				}
				double prop=(double)count/nb;
				if(prop==max)
				{
					candidats.add(i);
				}
				else
					if(prop>max)
					{
						max=prop;
						candidats.clear();
						candidats.add(i);
					}
				count_comp[i]=count;
			}
			count_try++;
			if(max==0 || count_try==10000) break;
			if( verbose > 1 )
			{
				std::cout<<"max="<<max<<" seq:";
				for(RandomIndex::iterator i=candidats.begin(); i!=candidats.end();i++)
					std::cout<<*i<<" ";
				std::cout<<std::endl;
			}
			numseq=candidats.draw();
			vec_strand[numseq]*=-1;
			for(int i=0;i<lendb;i++)
			{
				if(i!=numseq)
				{
					matches[numseq][i]*=-1;
					matches[i][numseq]*=-1;
				}
			}
		}

		if( verbose > 1 )
		{
			std::cout<<"Results:"<<std::endl;
			std::cout<<std::setw(6)<<"\\";
			for(int i=0;i<lendb;i++)
				std::cout<<std::setw(6)<<i+1;
			std::cout<<std::endl;
			for(int i=0;i<lendb;i++)
			{
				std::cout<<std::setw(6)<<i+1;
				for(int j=0;j<lendb;j++)
				{
					std::vector<int> seq_match(lendb);
					std::cout<<std::setw(6)<<matches[i][j];
				}
				std::cout<<std::endl;
			}
		}

		std::cout<<"saving oriented sequences..."<<std::endl;

		for(int i=0;i<lendb;i++)
		{
			if(count_comp[i]>0 && ! keepseq)
			{
				logfile<<"seq #"<<i+1<<":"<<db[i].getDE()<<" not oriented!"<<std::endl;
				if( verbose > 0 )
					std::cout<<"seq #"<<i+1<<":"<<db[i].getDE()<<" not oriented!"<<std::endl;
			}
			else
			{
				if(vec_strand[i]>0)
					fout<<db[i];
				else
				{
					SDGBioSeq seq=db[i];
					SDGString header=db[i].getDE()+" re-oriented";
					seq=seq.complement();
					seq.setDE(header);
					if( verbose > 0 )
						std::cout<<"seq #"<<i+1<<":"<<seq.getDE()<<"!"<<std::endl;
					fout<<seq;
				}
			}
		}
		fout.close();
	}
	catch(SDGException e)
	{
		std::cerr << e.message << std::endl;
	}

	std::cout<<"END orienter"<<std::endl;

	return 0;
}
