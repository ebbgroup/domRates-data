/*
 * utility.cpp
 *
 *  Created on: May 10, 2015
 *      Author: Carsten Kemena
 *	 Copyright: 2015
 *
 *  This file is part of BioSeqDataLib.
 *
 *  BioSeqDataLib is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  BioSeqDataLib is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with BioSeqDataLib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "utility.hpp"

namespace BioSeqDataLib {

/*
void
openInFile(const std::string &inFile, std::ifstream &inS)
{
	inS.open(inFile);
	if (inS.fail())
		throw std::ios_base::failure("Error opening file '" + inFile + "': " +strerror(errno));
}*/

/*
void
openOutFile(const std::string &outFile, std::ofstream &outS, std::ios_base::openmode mode)
{
	outS.open(outFile, mode);
	if (outS.fail())
		throw std::ios_base::failure("Error opening file '" + outFile + "': " +strerror(errno));
}*/

fs::path
getEnv(const std::string &name)
{
	  char* path = nullptr;
	  std::string pathS;
	  path = getenv (name.c_str());
	  if (path != nullptr)
		  pathS = path;
	  else
		  throw std::runtime_error("Environment variable " + name + " not set!");
	  return pathS;
}

} /* namespace BioSeqDataLib */
