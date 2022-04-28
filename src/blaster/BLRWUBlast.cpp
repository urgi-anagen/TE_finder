/*
 * \file BLRWUBlast.cpp
 */

#include <sstream>
#include "BLRWUBlast.h"

void BLRWUBlast::blast( int verbose )
{
    int sys_return = 0;
    std::string blast_command;
    if (para.getType() == "blastn") {
        blast_command = para.getType()
                        + " " + para.getBankCut() + " " + query_filename
                        + auto_blastparam + auto_blastnparam + " " + para.getOption()
                        + " O=" + result_filename + " > " + query_filename + "_blast.log 2>&1";
    }
    if (para.getType() == "blastx"
        || para.getType() == "tblastx"
        || para.getType() == "tblastn"
        || para.getType() == "blastp") {
        blast_command = para.getType()
                        + " " + para.getBankCut() + " " + query_filename
                        + auto_blastparam + " " + para.getOption() + " O=" + result_filename + " > " + query_filename +
                        "_blast.log 2>&1";
    }
    if (verbose > 0)
        std::cout << blast_command << std::endl;

    sys_return = system( blast_command.c_str() );

    if( verbose > 0 )
        std::cout<<"system call to blast return:"<<sys_return<<std::endl;

    if (sys_return == -1) {
        std::cout << std::ifstream(query_filename + "_blast.log").rdbuf();
        std::ostringstream ostr;
        ostr << " call to 'system()' function return error -1 !";
        throw SDGException(NULL, ostr.str(), -1);
    }

	if( sys_return != 0 )
  	{
        std::cout << std::ifstream(query_filename+"_blast.log").rdbuf();
  		std::ostringstream ostr;
  		ostr<<" program "<<para.getType()
  		<<" return with error value: "<<sys_return
  		<<", stopping blaster";
  		throw SDGException(NULL,ostr.str(),-1);
  	}

    SDGString rm_command="rm -f "+query_filename;
    system( rm_command );
    if(para.getCleanTmpFiles()){
        rm_command="rm -f "+query_filename+ "_blast.log";
        system( rm_command );
    }
}

void BLRWUBlast::pressdb( int verbose )
{
	SDGString command;
	if(para.getType()=="blastx"
			||para.getType()=="blastp")
	{
		std::ifstream in(para.getBankCut()+".xpd");
		if(!in)
		{
			command="xdformat -p "+para.getBankCut();
		  	if( verbose > 0 )
		  		std::cout<<command<<std::endl;
			system(command.start());
		}
		else
			if( verbose > 0 )
				std::cout<<"file '"<<para.getBankCut()<<"' already formatted"<<std::endl;
	}
	else
	{
		std::ifstream in(para.getBankCut()+".xnd");
		if(!in)
		{
			command="xdformat -n "+para.getBankCut();
		  	if( verbose > 0 )
		  		std::cout<<command<<std::endl;
			system(command.start());
		}
		else
			if( verbose > 0 )
				std::cout<<"file '"<<para.getBankCut()<<"' already formatted"<<std::endl;
	}
}
