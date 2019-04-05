/*
 * BlastHit.h
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
#ifndef BLASTHIT_H_
#define BLASTHIT_H_

#include <string>
#include <vector>
#include <map>

#include "../utility/stringHelpers.hpp"
#include "../external/Input.hpp"
/**
 * \file BlastHit.hpp
 * \brief File containing the class definition for a BLAST hit.
 */
namespace BioSeqDataLib
{

/**
 * \brief Class representing a BLAST hit.
 */
class BlastHit
{
public:
	std::string queryID;
	std::string subjectID;
	float percentIdentity;
	size_t alignmentLength;
	size_t nMismatches;
	size_t nGapOpenings;
	size_t queryStart;
	size_t queryEnd;
	char queryStrand;
	size_t subjectStart;
	size_t subjectEnd;
	char subjectStrand;
	double eValue;
	float HSPBitScore;

	/**
	 * Standard constructor
	 */
	BlastHit();

	/**
	 *
	 * @param line
	 */
	explicit BlastHit(std::string line);
	virtual ~BlastHit();
};


/**
 * \brief Reading the blast output.
 * @param inFile The input file.
 * @param hits Stores the hits. key=queryName, value=vector of hits
 */
void
readBlast(const std::string &inFile, std::map<std::string, std::vector<BlastHit> > &hits);



} /* namespace BioSeqDataLib */

#endif /* BLASTHIT_H_ */
