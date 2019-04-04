/**
 *
 * BLRMatchMapLoader.cpp
 *
 **/
#include <stdlib.h>
#include <regex.h>
#include <fstream>
#include "FragAlign.h"
#include "BLRMatchMapLoader.h"
#include <SDGError.h>

//---------------------------------------------------------------------------
void BLRMatchMapLoader::readAlign(BLRMatchMap& blrmm, std::istream& input_align, int verbose)
{
   unsigned countseqS=0,countseqQ=0;
 
   blrmm.getName2NumQ().clear();
   blrmm.getName2NumS().clear();
   blrmm.getNum2NameQ().clear();	
   blrmm.getNum2NameS().clear();	
   
   std::map<std::string,long> name2numQ = blrmm.getName2NumQ();
   std::map<std::string,long> name2numS = blrmm.getName2NumS();
   std::map<long,std::string> num2nameQ = blrmm.getNum2NameQ();	
   std::map<long,std::string> num2nameS = blrmm.getNum2NameS();

    if (blrmm.getParameter()->getBank() == blrmm.getParameter()->getQuery() &&
        blrmm.getParameter()->getBank() != "<not set>")
        blrmm.setSameDb(true);

    //Check format
    std::string str;
    char buff[1024];
    input_align.getline(buff, 1023, '\n');
    str = buff;
    size_t n = std::count(str.begin(), str.end(), '\t');
    if (n != 8) {
        std::cout << "BLRMatchMap ERROR: Number of columns is:" << n + 1 << std::endl
                  << "first line:\n" << str << std::endl
                  << "! not an *align* formated file" << std::endl;
        exit(0);
    }

    //Read the first used to test the file format
    std::istringstream is(str);
    RangePair rp;
    rp.readtxt(is);
    if (blrmm.getParameter()->getEvalFilter() > rp.getE_value()
        && blrmm.getParameter()->getIdFilter() < rp.getIdentity()
        && blrmm.getParameter()->getLenFilter() < rp.getLength()) {
        std::map<std::string, long>::iterator it
                = name2numQ.find(rp.getRangeQ().getNameSeq());
        if (it == name2numQ.end()) {
            name2numQ[rp.getRangeQ().getNameSeq()] = ++countseqQ;
            num2nameQ[countseqQ] = rp.getRangeQ().getNameSeq();
            rp.getRangeQ().setNumChr(countseqQ);
        } else
            rp.getRangeQ().setNumChr(it->second);

        if (blrmm.isSameDb()) {
            it = name2numQ.find(rp.getRangeS().getNameSeq());
            if (it == name2numQ.end()) {
                name2numQ[rp.getRangeS().getNameSeq()] = ++countseqQ;
                num2nameQ[countseqQ] = rp.getRangeS().getNameSeq();
                rp.getRangeS().setNumChr(countseqQ);
            } else
                rp.getRangeS().setNumChr(it->second);
        } else {
            it = name2numS.find(rp.getRangeS().getNameSeq());
            if (it == name2numS.end()) {
                name2numS[rp.getRangeS().getNameSeq()] = ++countseqS;
                num2nameS[countseqS] = rp.getRangeS().getNameSeq();
                rp.getRangeS().setNumChr(countseqS);
            } else
                rp.getRangeS().setNumChr(it->second);
        }
        blrmm.insert(rp);
    }

    // read the others
    while (input_align) {
        RangePair rp;
        rp.readtxt(input_align);
        if (input_align) {
            if (blrmm.getParameter()->getEvalFilter() < rp.getE_value()
                || blrmm.getParameter()->getIdFilter() > rp.getIdentity()
                || blrmm.getParameter()->getLenFilter() > rp.getLength())
                continue;

            std::map<std::string, long>::iterator it
                    = name2numQ.find(rp.getRangeQ().getNameSeq());
            if (it == name2numQ.end()) {
                name2numQ[rp.getRangeQ().getNameSeq()] = ++countseqQ;
                num2nameQ[countseqQ] = rp.getRangeQ().getNameSeq();
                rp.getRangeQ().setNumChr(countseqQ);
            } else
                rp.getRangeQ().setNumChr(it->second);

            if (blrmm.isSameDb()) {
                it = name2numQ.find(rp.getRangeS().getNameSeq());
                if (it == name2numQ.end()) {
                    name2numQ[rp.getRangeS().getNameSeq()] = ++countseqQ;
                    num2nameQ[countseqQ] = rp.getRangeS().getNameSeq();
                    rp.getRangeS().setNumChr(countseqQ);
                } else
                    rp.getRangeS().setNumChr(it->second);
            } else {
                it = name2numS.find(rp.getRangeS().getNameSeq());
                if (it == name2numS.end()) {
                    name2numS[rp.getRangeS().getNameSeq()] = ++countseqS;
                    num2nameS[countseqS] = rp.getRangeS().getNameSeq();
                    rp.getRangeS().setNumChr(countseqS);
                } else
                    rp.getRangeS().setNumChr(it->second);
            }
            blrmm.insert(rp);
        }
    }

    blrmm.setName2NumQ(name2numQ);
    blrmm.setName2NumS(name2numS);
    blrmm.setNum2NameQ(num2nameQ);
    blrmm.setNum2NameS(num2nameS);

    if (verbose > 0) {
        std::cout << "nb of matches: " << blrmm.getNbMatchesInMapAlign() << std::endl;
        if (blrmm.isSameDb())
            std::cout << "nb of distinct queries/subjects: " << blrmm.getNbQseq() << std::endl;
        else {
            std::cout << "nb of distinct queries: " << blrmm.getNbQseq() << std::endl;
            std::cout << "nb of distinct subjects: " << blrmm.getNbSseq() << std::endl;
        }
    }
}
//---------------------------------------------------------------------------
void BLRMatchMapLoader::readPath(BLRMatchMap& blrmm, std::istream& input_path, int verbose)
{
	unsigned countseqS=0,countseqQ=0;

	blrmm.getName2NumQ().clear();
	blrmm.getName2NumS().clear();
	blrmm.getNum2NameQ().clear();
	blrmm.getNum2NameS().clear();

	std::map<std::string,long> name2numQ = blrmm.getName2NumQ();
	std::map<std::string,long> name2numS = blrmm.getName2NumS();
	std::map<long,std::string> num2nameQ = blrmm.getNum2NameQ();
	std::map<long,std::string> num2nameS = blrmm.getNum2NameS();

	if(blrmm.getParameter()->getBank()==blrmm.getParameter()->getQuery() &&
		blrmm.getParameter()->getBank()!="<not set>")
		blrmm.setSameDb(true);

	unsigned long currentId=0, previousId=0;
	RangePairSet rps;

	while(input_path)
	{
		RangePairSet rp;
		rp.readtxt(input_path);
		if(input_path)
		{
			if(blrmm.getParameter()->getEvalFilter()<rp.getE_value()
				|| blrmm.getParameter()->getIdFilter()>rp.getIdentity()
				|| blrmm.getParameter()->getLenFilter()>rp.getLength())
				continue;

			std::map<std::string,long>::iterator it
			=name2numQ.find(rp.getRangeQ().getNameSeq());
			if(it==name2numQ.end()) //unknown name -> add to name list to associate a number
			{
				name2numQ[rp.getRangeQ().getNameSeq()]=++countseqQ;
				num2nameQ[countseqQ]=rp.getRangeQ().getNameSeq();
				rp.getRangeQ().setNumChr(countseqQ);
			}
			else
			  rp.getRangeQ().setNumChr(it->second); // name already found -> add number

			if(blrmm.isSameDb()) // if QueryDB and Subject DB are the same
			{
				it=name2numQ.find(rp.getRangeS().getNameSeq());
				if(it==name2numQ.end())
				{
					name2numQ[rp.getRangeS().getNameSeq()]=++countseqQ;
					num2nameQ[countseqQ]=rp.getRangeS().getNameSeq();
					rp.getRangeS().setNumChr(countseqQ);
				}
				else
					rp.getRangeS().setNumChr(it->second);
			}
			else // QueryDB and Subject DB are different
			{
				it=name2numS.find(rp.getRangeS().getNameSeq());
				if(it==name2numS.end())
				{
					name2numS[rp.getRangeS().getNameSeq()]=++countseqS;
					num2nameS[countseqS]=rp.getRangeS().getNameSeq();
					rp.getRangeS().setNumChr(countseqS);
				}
				else
					rp.getRangeS().setNumChr(it->second);
			}

			//remove sequence name to gain space
			rp.getRangeQ().setNameSeq("");
			rp.getRangeS().setNameSeq("");

			currentId = rp.getId();
			if (previousId != currentId && previousId!=0) //new path
			{
			//save previous path
				rps.setPath(blrmm.para->getDist_pen(),0.0,blrmm.para->getGap_pen());
				blrmm.rpsList.push_back(rps);

			//clear for new path
				rps.clear();
			}
			rps.addPath(rp);
			previousId=currentId;
		}	
	}
	rps.setPath(blrmm.para->getDist_pen(),0.0,blrmm.para->getGap_pen());
	blrmm.rpsList.push_back(rps);
	blrmm.setName2NumQ(name2numQ);
	blrmm.setName2NumS(name2numS);
	blrmm.setNum2NameQ(num2nameQ);
	blrmm.setNum2NameS(num2nameS);

	if(verbose>0)
	{
		std::cout<<"nb of matches: "<<blrmm.rpsList.size()<<std::endl;
		if( blrmm.isSameDb() )
			std::cout<<"nb of distinct queries/subjects: "<<blrmm.getNbQseq()<<std::endl;
		else
		{
			std::cout<<"nb of distinct queries: "<<blrmm.getNbQseq()<<std::endl;
			std::cout<<"nb of distinct subjects: "<<blrmm.getNbSseq()<<std::endl;
		}
	}
}

