/*
 * BlastHit.cpp
 *
 *  Created on: 29 Jul 2014
 *      Author: ckemena
 *	 Copyright: 2014
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
#include "BlastHit.hpp"

namespace BioSeqDataLib
{

BlastHit::BlastHit()
{
	// TODO Auto-generated constructor stub

}


BlastHit::BlastHit(std::string line)
{
	std::vector<std::string> tokens;
	split(line, "\t", tokens);
	queryID = tokens[0];
	subjectID = tokens[1];
	percentIdentity = std::stof(tokens[2]);
	alignmentLength = std::stoi(tokens[3]);
	nMismatches = std::stoi(tokens[4]);
	nGapOpenings = std::stoi(tokens[5]);
	queryStart = std::stoi(tokens[6])-1;
	queryEnd = std::stoi(tokens[7])-1;
	subjectStart = std::stoi(tokens[8])-1;
	subjectEnd = std::stoi(tokens[9])-1;
	eValue = std::stod(tokens[10]);
	HSPBitScore = std::stod(tokens[11]);
	if (queryStart > queryEnd)
	{
		queryStrand = '-';
		std::swap(queryStart, queryEnd);
	}
	else
		queryStrand = '+';
	if (subjectStart > subjectEnd)
	{
		subjectStrand = '-';
		std::swap(subjectStart, subjectEnd);
	}
	else
		subjectStrand = '+';
}

BlastHit::~BlastHit()
{
	// TODO Auto-generated destructor stub
}


void
readBlast(const std::string &inFile, std::map<std::string, std::vector<BlastHit> > &hits)
{

	AlgorithmPack::Input inS(inFile);

	std::string line;
	size_t lineN=0;
	try
	{
		size_t pos;
		while (getline(inS, line))
		{
			++lineN;
			if ((line.empty()) || (line[0] == '#'))
				continue;


			if ( (pos = line.find_first_of('\t')) != std::string::npos)
				hits[line.substr(0,pos)].emplace_back(line);
			else
				throw std::runtime_error("Error parsing line: " + std::to_string(lineN) + " of file: " + inFile);
		}
	}
	catch (const std::exception &e)
	{
		throw std::runtime_error("Error parsing line: " + std::to_string(lineN) + " of file: " + inFile);
	}
}


} /* namespace BioSeqDataLib */
