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
  cout << "     -o: name of the ouput file (default=inFileName+'.gff')" << endl;
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

void align2piler( ifstream &inFile, ofstream &outFile, int verbose )
{
  int countLine = 0;

  while ( ! inFile.eof() )
    {
      string line;
      getline( inFile, line );
      if ( line == "" )
	break;
      countLine += 1;
      if ( verbose > 2 )
 	cout << line << endl;

      vector<string> tokens;
      tokenize( line, tokens, "\t" );

      string qryName = tokens[0];
      int qryStart = atoi(tokens[1].c_str());
      int qryEnd = atoi(tokens[2].c_str());
      string sbjName = tokens[3];
      int sbjStart = atoi(tokens[4].c_str());
      int sbjEnd = atoi(tokens[5].c_str());
      string evalue = tokens[6];
      string score = tokens[7];
      float identity = atof(tokens[8].c_str());

      string strand;
      if ( qryStart  < qryEnd )
	{
	  if ( sbjStart < sbjEnd )
	    strand = "+";
	  else {
	    strand = "-";
	    int tmp = sbjStart;
	    sbjStart = sbjEnd;
	    sbjEnd = tmp; }
	}
      else if ( qryStart > qryEnd )
	{
	  int tmp = qryStart;
	  qryStart = qryEnd;
	  qryEnd = tmp;
	  if ( sbjStart < sbjEnd )
	    strand = "-";
	  else {
	    strand = "+";
	    int tmp = sbjStart;
	    sbjStart = sbjEnd;
	    sbjEnd = tmp; }
	}

      float divergence = (100.0 - identity) / 100.0;

      stringstream outLine;
      outLine << qryName << "\tblaster\thit\t" << qryStart << "\t" << qryEnd << "\t" << score << "\t" << strand << "\t.\tTarget " << sbjName << " " << sbjStart << " " << sbjEnd << "; maxe ";
      outLine.precision(4);
      outLine << fixed << divergence << "\n";

      if ( verbose > 1 )
	cout << outLine.str();

      outFile << outLine.str();
      if( countLine % 500 == 0 )
	{
	  outFile << flush;
	  if( verbose > 1 )
	    cout << countLine << endl;
	}
    }
  if ( verbose > 0 )
    cout << "number of lines: " << countLine << endl;
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
    outFileName = inFileName + ".gff";
  if ( verbose > 0 )
    cout << "output file: " << outFileName << endl;

  ofstream outFile;
  outFile.open( outFileName.c_str() );
  if ( ! outFile.is_open() ) {
    cerr << "*** Error: unable to open output file" << endl;
    exit( EXIT_FAILURE ); }

  align2piler( inFile, outFile, verbose );

  inFile.close();
  if ( outFileName != "" )
    outFile.close();

  return EXIT_SUCCESS;
}
