/*
 * runDBsearch.hpp
 *
 *  Created on: May 11, 2015
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

/**
 * @file runDBsearch.hpp
 * \brief This progam contains functions to handle database searches
 *
 */

#ifndef RUNDBSEARCH_HPP_
#define RUNDBSEARCH_HPP_

#include <string>
#include <stdexcept>

#include <boost/filesystem.hpp>

#include "../annotation/BlastHit.hpp"

namespace fs=boost::filesystem;

namespace BioSeqDataLib
{

/**
 * \brief Makes a blast database
 * @param file The sequence file
 * @param protein if true a protein database is created
 * @param logFile The logfile for the run. If empty, output will be redirected to /dev/null.
 */
void
makeBlastDB(const fs::path &file, bool protein=false, const fs::path &logFile="");

/**
 * \brief Runs a blast
 * @param program The blast program to use
 * @param query The query file
 * @param db The database path
 * @param outFile The output file
 * @param log The logfile for the run. If empty, output will be redirected to /dev/null.
 */
void
runblast(const std::string &program, const fs::path &query, const fs::path &db, const fs::path &outFile, const fs::path &log);

// void
// runblast(const std::string &program, const std::string &db, BlastHit &hits, const std::string &outFile);
//

}



#endif /* RUNDBSEARCH_HPP_ */
