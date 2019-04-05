/*
 * OrthologySet.hpp
 *
 *  Created on: 18 Aug 2014
 *      Author: Carsten Kemena
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

/**
 * \file OrthologySet.hpp
 * \brief File containing the OrthologySet class.
 */
#ifndef ORTHOLOGYSET_HPP_
#define ORTHOLOGYSET_HPP_

// C++ header
#include <fstream>
#include <map>
#include <string>
#include <set>
#include <vector>
#include <iostream>

// BioSeqDataLib headers
#include "../external/Input.hpp"
#include "../utility/stringHelpers.hpp"

namespace AP=AlgorithmPack;

namespace BioSeqDataLib
{

/**
 * key: species name
 * value: vector of all genes
 */
typedef std::map<std::string, std::vector<std::string> > OrthoGroup;

/**
 * \brief Class to store orthoglogy groups.
 */
class OrthologySet
{
private:
	std::set<std::string> speciesSet_;
	std::map<std::string, OrthoGroup> groups_;
	std::map<std::string, std::map<std::string, std::string> > id2group_;

	void
	readOrthoMCL_(AP::Input &orthoS);

	void
	readProteinortho_(AP::Input &orthoS);
public:

	/**
	 * \brief Standard constructor
	 */
	OrthologySet();

	/**
	 * \brief Standard destructor
	 */
	virtual ~OrthologySet();

	/**
	 * \brief Reads an orthology file.
	 * @param orthoFile The file containing the orthology data.
	 */
	void
	read(const std::string &orthoFile);

	/**
	 *
	 * @param groupID
	 * @return
	 */
	OrthoGroup &operator[](const std::string &groupID)
	{
		return groups_[groupID];
	}

	/**
	 *
	 * @param groupID
	 * @return
	 */
	const OrthoGroup &operator[](const std::string &groupID) const
	{
		return groups_.at(groupID);
	}

	/**
	 *
	 * @param species
	 * @param seqID
	 * @return
	 */
	const OrthoGroup &group(const std::string &species, const std::string &seqID) const
	{
		 return groups_.at(id2group_.at(species).at(seqID));
	}

	const OrthoGroup *find(const std::string &species, const std::string &seqID) const;

	std::set<std::string>
	speciesSet()
	{
		return speciesSet_;
	}
};

} /* namespace BioSeqDataLib */

#endif /* ORTHOLOGYSET_HPP_ */
