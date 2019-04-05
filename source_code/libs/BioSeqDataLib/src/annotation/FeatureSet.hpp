/*
 * FeatureSet.hpp
 *
 *  Created on: 29 Nov 2013
 *      Author: Carsten Kemena
 *		 Email: c.kemena[@]uni-muenster.de
 *	 Copyright: 2013
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
 * \file FeatureSet.hpp
 * \brief Header for the FeatureSet class.
 */
#ifndef FEATURESET_HPP_
#define FEATURESET_HPP_

// C++ header
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <iostream>

// BioSeqDataLib header
#include "Feature.hpp"
#include "../utility/Exceptions.hpp"


namespace BioSeqDataLib
{



/**
 * \brief Class to store all features of a GFF file
 */
class FeatureSet
{

private:

	// Key = ID
	// Value = pair(seqID, featureType)
	std::map<std::string, std::string> id2position_;

	//Key = seqID
	//Value = map(key=type, key=parentId, value=vector<feature>)
	std::map<std::string, std::vector<Feature> > features_;

	std::map<std::string, size_t> typeCounter;

public:
	/**
	 * \brief Standard Constructor
	 */
	FeatureSet();

	/**
	 * \brief Standard Destructor
	 */
	virtual ~FeatureSet();

	/**
	 * \brief Reads a GFF file.
	 * @param infile The filename.
	 */
	void
	readGff(const std::string &infile);

	/**
	 * \brief Accessing a feature.
	 * @param id The feature id to search for.
	 * @return Iterator to the feature.
	 */
	std::vector<Feature>::iterator
	operator[](const std::string &id)
	{
		const std::string &seqId = id2position_.at(id);
		std::vector<Feature> &vec = features_[seqId];
		for (size_t i=0; i<vec.size(); ++i)
			if (vec[i].attribute("ID") == id)
				return vec.begin()+i;
		return vec.end();
	}


	/**
	 * \brief Accessing a feature.
	 * @param id The feature id to search for.
	 * @return Iterator to the feature.
	 */
	std::vector<Feature>::const_iterator
	operator[](const std::string &id) const
	{
		const std::string &seqId = id2position_.at(id);
		const std::vector<Feature> &vec = features_.at(seqId);
		for (size_t i=0; i<vec.size(); ++i)
			if (vec[i].attribute("ID") == id)
				return vec.begin()+i;
		return vec.end();
	}

	/**
	 * \brief Returns the number of times a type occurs.
	 * @param type The type to count.
	 * @return The number of times the type occurs.
	 */
	size_t
	counts(const std::string &type) const
	{
		return typeCounter.at(type);
	}


};

} /* namespace BioSeqDataLib */

#endif /* FEATURESET_HPP_ */



