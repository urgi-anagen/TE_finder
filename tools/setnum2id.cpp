#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
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
  cout << "     -i: name of the input file (format='set')" << endl;
  cout << "     -o: name of the output file (default=inFileName+'.out')" << endl;
  cout << "     -v: verbose (default=0/1/2)" << endl;
  cout << "" << endl;
}

//-----------------------------------------------------------------------------

void parseArgs ( int argc, char **argv, string &inFileName, string &outFileName, int &verbose )
{
  extern char *optarg;
  char c;
  while ( (c = getopt(argc,argv,"hi:o:v:")) != EOF )
    {
      switch ( c )
        {
        case 'h':
          help( argv[0] );
	  exit( EXIT_SUCCESS );
        case 'i':
	  inFileName = optarg;
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

void setnum2id( ifstream &inFile, ofstream &outFile, int verbose )
{
  map<string,int> dID2count;
  map<string,int>::iterator iter;
  int countID = 1;
  int countLine = 0;
  string line;

  while ( ! inFile.eof() )
    {

      getline( inFile, line );
      if ( line == "" )
	break;
      countLine += 1;
       if ( verbose > 2 )
 	cout << line << endl;

      vector<string> tokens;
      tokenize( line, tokens, "\t" );
      string ID = tokens[0] + "-" + tokens[2] + "-" + tokens[1];
      //cout << ID << endl;
      //dID2count[ ID ] = atoi( tokens[0].c_str() );

      iter = dID2count.find( ID );
      int newPath;
      if( iter == dID2count.end() )
	{
	  newPath = countID;
	  countID += 1;
	  dID2count[ ID ] = newPath;
	}
      else
	newPath = dID2count[ ID ];

      ostringstream newPath_str;
      newPath_str << newPath;
      string outLine = newPath_str.str() + "\t" + tokens[1] + "\t" + tokens[2] + "\t" + tokens[3] + "\t" + tokens[4] + "\n";

      if ( ! outFile.is_open() )
	cout << outLine;
      else
	outFile << outLine;
      if( countLine % 500 == 0 )
	{
	  outFile << flush;
	  if( verbose > 1 )
	    cout << countLine << endl;
	}
    }
  if ( verbose > 0 ) {
    cout << "number of lines: " << countLine << endl;
    cout << "distinct IDs: " << dID2count.size() << endl; }
}

//-----------------------------------------------------------------------------

int main ( int argc, char **argv )
{
  string inFileName = "";
  string outFileName = "";
  int verbose = 0;

  parseArgs( argc, argv, inFileName, outFileName, verbose );
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
    outFileName = inFileName + ".out";
  if ( verbose > 0 )
    cout << "output file: " << outFileName << endl;

  ofstream outFile;
  outFile.open( outFileName.c_str() );
  if ( ! outFile.is_open() ) {
    cerr << "*** Error: unable to open output file" << endl;
    exit( EXIT_FAILURE ); }

  setnum2id( inFile, outFile, verbose );

  inFile.close();
  if ( outFileName != "" )
    outFile.close();

  return EXIT_SUCCESS;
}
