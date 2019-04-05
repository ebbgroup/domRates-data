/*
 * utility.hpp
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

#include <ios>
#include <fstream>
#include <string>
#include <errno.h>
#include <string.h>
#include <ios>
#include <iostream>

#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

#ifndef UTILITY_HPP_
#define UTILITY_HPP_

/**
 * \file utility.hpp
 * \brief File containing a few helper functions.
 */

namespace BioSeqDataLib {

/**
 * \brief Opening a file for writing.
 * Throws an exception in case of failure.
 * @param inFile The file to open.
 * @return An ofstream object.
 */
void
openOutFile(const std::string &outFile, std::ofstream &outS, std::ios_base::openmode mode = std::ios_base::out);



fs::path
getEnv(const std::string &name);

} /* namespace BioSeqDataLib */

#endif /* UTILITY_HPP_ */
