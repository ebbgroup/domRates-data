/*
 * rundDBsearch.cpp
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


#include "runDBsearch.hpp"

namespace BioSeqDataLib
{

void
makeBlastDB(const fs::path &file, bool protein, const fs::path &logFile)
{
	std::string logF;
	if (logFile.empty())
		logF =  " > /dev/null";
	else
		logF = " >> " + logFile.string() +" 2>> " + logFile.string();

	std::string dbtype = protein ? "prot" : "nucl";
	std::string command = "makeblastdb -in " + file.string() + " -dbtype " + dbtype +logF;
	if (system(command.c_str()))
	{
		throw std::runtime_error("Error: makeblastdb on " + file.string() + " failed!");
	}
}

void
runblast(const std::string &program, const fs::path &query, const fs::path &db, const fs::path &outFile, const fs::path &logFile)
{
	if (outFile.empty())
	{
		// todo
	}

	std::string logF;
	if (logFile.empty())
		logF =  " > /dev/null";
	else
		logF = " >> " + logFile.string() +" 2>> " + logFile.string();

	std::string command = program + " -query " + query.string() + " -db " + db.string() + " -outfmt 6 -out " + outFile.string() + logF;
	if (system(command.c_str()))
	{
		throw std::runtime_error("Error: Blast of sequence " + query.string() + " failed!");
	}

}

// void
// runblast(const std::string &program, const std::string &db, BlastHit &hits, const std::string &outFile)
// {
//
//
// }




}
