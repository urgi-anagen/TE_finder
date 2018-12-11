/*
 * \file BLRWUBlast.cpp
 */

#include <sstream>
#include "BLRWUBlast.h"

void BLRWUBlast::blast( int verbose )
{
	int sys_return = 0;
	if( para.getType() == "blastn" )
	{
		SDGString blast_command=para.getType()
		+" "+para.getBankCut()+" "+query_filename
		+auto_blastparam+auto_blastnparam+" "+para.getOption()
		+" O="+result_filename;
		if( verbose > 0 )
			std::cout<<blast_command<<std::endl;
		sys_return = system( blast_command );
	}
	if( para.getType() == "blastx"
			|| para.getType() == "tblastx"
			|| para.getType() == "tblastn"
			||para.getType() == "blastp" )
	{
		SDGString blast_command=para.getType()
		+" "+para.getBankCut()+" "+query_filename
		+auto_blastparam+" "+para.getOption()+" O="+result_filename;
		if( verbose > 0 )
			std::cout<<blast_command<<std::endl;
		sys_return = system( blast_command );
	}

	SDGString rm_command = "rm -f " + query_filename;
	system( rm_command );

	if( sys_return == -1 )
		throw SDGException(NULL," fork error!, stopping program.",-1);

	if( sys_return != 0 )
	{
		if( ( para.getType()!="tblastx" && para.getType()!="blastx" )
				|| ( para.getType()=="tblastx" && sys_return!=4096 )
				|| ( para.getType()=="blastx" && sys_return!=5888 ) )
  	{
  		std::ostringstream ostr;
  		ostr<<" program "<<para.getType()
  		<<" return with error value: "<<sys_return
  		<<", stopping blaster";
  		throw SDGException(NULL,ostr.str(),-1);
  	}
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
