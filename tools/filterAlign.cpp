#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
using namespace std;

//-----------------------------------------------------------------------------

void help ( char *programName )
{
  cout << "" << endl;
  cout << "usage: " << programName << " [ options ]" << endl;
  cout << "options:" << endl;
  cout << "     -h: this help" << endl;
  cout << "     -i: name of the input file (format='align')" << endl;
  cout << "     -E: maximum E-value (default=100)" << endl;
  cout << "     -S: minimum score (default=0)" << endl;
  cout << "     -I: minimum identity (default=0)" << endl;
  cout << "     -l: minimum length (default=0)" << endl;
  cout << "     -L: maximum length (default=1000000000)" << endl;
  cout << "     -o: name of the output file (default=inFileName+'.filtered')" << endl;
  cout << "     -v: verbose (default=0/1)" << endl;
  cout << "" << endl;
}

//-----------------------------------------------------------------------------

void parseArgs ( int argc, char **argv, string &inFileName, float &maxEvalue, int &minScore, int &minIdentity, int &minLength, int &maxLength, string &outFileName, int &verbose )
{
  extern char *optarg;
  char c;
  while ( (c = getopt(argc,argv,"hi:E:S:I:l:L:o:v:")) != EOF )
  {
	  switch ( c )
	  {
	  case 'h':
		  help( argv[0] );
		  exit( EXIT_SUCCESS );
	  case 'i':
		  inFileName = optarg;
		  break;
	  case 'E':
		  maxEvalue = atof(optarg);
		  break;
	  case 'S':
		  minScore = atoi(optarg);
		  break;
	  case 'I':
		  minIdentity = atoi(optarg);
		  break;
	  case 'l':
		  minLength = atoi(optarg);
		  break;
	  case 'L':
		  maxLength = atoi(optarg);
		  break;
	  case 'o':
		  outFileName = optarg;
		  break;
	  case 'v':
		  verbose = atoi(optarg);
		  break;
	  }
  }
}

//-----------------------------------------------------------------------------

void tokenize( const string& str, vector<string>& tokens, const string& delimiters = " " )
{
  // skip delimiters at beginning
  string::size_type lastPos = str.find_first_not_of( delimiters, 0 );

  // find first "non-delimiter"
  string::size_type pos = str.find_first_of( delimiters, lastPos );

  while ( string::npos != pos || string::npos != lastPos )
    {
      // found a token, add it to the vector.
      tokens.push_back( str.substr( lastPos, pos - lastPos ) );

      // skip delimiters.  Note the "not_of".
      lastPos = str.find_first_not_of( delimiters, pos );

      // find next "non-delimiter"
      pos = str.find_first_of( delimiters, lastPos );
    }
}

//-----------------------------------------------------------------------------

void filterAlign( ifstream &inFile, float &maxEvalue, int &minScore, int &minIdentity, int &minLength, int &maxLength, ofstream &outFile, int verbose )
{
  int nbMatches = 0;
  int nbFiltered = 0;

  while ( ! inFile.eof() )
  {
	  string line;
	  getline( inFile, line );
	  if ( line == "" )
		  break;
	  ++ nbMatches;

      vector<string> tokens;
      tokenize( line, tokens, "\t" );

      string qryName = tokens[0];
      int qryStart = atoi(tokens[1].c_str());
      int qryEnd = atoi(tokens[2].c_str());
      string sbjName = tokens[3];
      int sbjStart = atoi(tokens[4].c_str());
      int sbjEnd = atoi(tokens[5].c_str());
      float evalue = atof( tokens[6].c_str() );
      int score = atoi( tokens[7].c_str() );
      float identity = atof(tokens[8].c_str());

      int length;
      if ( qryStart  < qryEnd )
    	  length = qryEnd - qryStart + 1;
      else
    	  length = qryStart - qryEnd + 1;

      if ( evalue <= maxEvalue && score >= minScore && identity >= minIdentity && length >= minLength && length <= maxLength )
      {
    	  stringstream outLine;
    	  outLine << qryName << "\t" << qryStart << "\t" << qryEnd << "\t" << sbjName << "\t" << sbjStart << "\t" << sbjEnd << "\t" << evalue << "\t" << score << "\t" << identity << "\n";
    	  if ( verbose > 1 )
    		  cout << outLine.str();
    	  outFile << outLine.str();
    	  if( nbMatches % 500 == 0 )
    	  {
    		  outFile << flush;
    		  if( verbose > 1 )
    			  cout << nbMatches << endl;
    	  }
      }
      else
    	  ++ nbFiltered;
  }
  if ( verbose > 0 )
  {
	  cout << "total number of matches: " << nbMatches << endl;
	  cout << "number of filtered matches: " << nbFiltered << endl;
  }
}

//-----------------------------------------------------------------------------

int main ( int argc, char **argv )
{
  string inFileName = "";
  float maxEvalue = 100;
  int minScore = 0;
  int minIdentity = 0;
  int minLength = 0;
  int maxLength = 1000000000;
  string outFileName = "";
  int verbose = 0;

  parseArgs( argc, argv, inFileName, maxEvalue, minScore, minIdentity, minLength, maxLength, outFileName, verbose );
  if ( inFileName == "" )
    {
      cerr << "*** Error: missing name of the input file" << endl;
      help( argv[0] );
      exit( EXIT_FAILURE );
    }
  if ( verbose > 0 )
    cout << "input file: " << inFileName << endl;

  ifstream inFile;
  inFile.open( inFileName.c_str() );
  if ( ! inFile.is_open() )
    {
      cerr << "*** Error: unable to open input file " << inFileName << endl;
      exit( EXIT_FAILURE );
    }

  if ( outFileName == "" )
    outFileName = inFileName + ".filtered";
  if ( verbose > 0 )
    cout << "output file: " << outFileName << endl;

  ofstream outFile;
  outFile.open( outFileName.c_str() );
  if ( ! outFile.is_open() ) {
    cerr << "*** Error: unable to open output file" << endl;
    exit( EXIT_FAILURE ); }

  filterAlign( inFile, maxEvalue, minScore, minIdentity, minLength, maxLength, outFile, verbose );

  inFile.close();
  if ( outFileName != "" )
    outFile.close();

  return EXIT_SUCCESS;
}
