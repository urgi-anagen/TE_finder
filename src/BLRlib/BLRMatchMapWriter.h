/*
 * BLRMatchMapWriter.h
 *
 *  Created on: 30 sept. 2015
 *      Author: hquesnev
 */

#ifndef BLRMATCHMAPWRITER_H_
#define BLRMATCHMAPWRITER_H_

#include <iostream>
#include <stdlib.h>
#include <SDGString.h>
#include <utility>
#include <map>
#include "BLRMatchMap.h"

class BLRMatchMapWriter
{
	public:
  		void writePath(const RpsList& rpsList, const std::map<long,std::string>& num2nameQ, const std::map<long,std::string>& num2nameS,std::ostream& out);
  		void writeGFF3(const RpsList& rpsList, const std::map<long,std::string>& num2nameQ, const std::map<long,std::string>& num2nameS,std::ostream& out, const std::string& source);

};


#endif /* BLRMATCHMAPWRITER_H_ */
