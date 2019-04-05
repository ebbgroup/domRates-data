/*
 * Feature.h
 *
 *  Created on: 3 Jul 2014
 *      Author: Carsten Kemena
 *		 Email: c.kemena[@]uni-muenster.de
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
 * \file Feature.hpp
 * \brief A class to store a Feature.
 */
#ifndef FEATURE_H_
#define FEATURE_H_

#include<vector>
#include<string>
#include<map>
#include<limits>
#include<fstream>
#include<iostream>

#include "../utility/stringHelpers.hpp"

namespace BioSeqDataLib {

/**
 * \brief Stores a single line of a GFF file.
 *
 */
class Feature
{

private:


public:
	std::string seqId_;
	std::string source_;
	std::string type_;
	size_t start_;
	size_t end_;
	double score_;
	char strand_;
	short phase_;
	std::vector<std::string> parents_;
	std::map<std::string, std::string> attributes_;

	/**
	 * \brief Standard constructor
	 */
	Feature() = default;

	/**
	 * \brief Standard copy constructor
	 * @param obj Object to copy
	 */
	Feature(const Feature&) = default;

	/**
	 * \brief Standard move constructor
	 * @param obj Object to move
	 */
	Feature(Feature&&) = default;

	/**
	 * \brief Constructor
	 * @param line A single line of a GFF file
	 */
	explicit Feature(const std::string &line);

    /**
     * \brief constructor
     * @param seq The sequence id
     * @param source The source
     * @param type The type
     * @param start The start
     * @param end The end
     * @param score The score
     * @param strand The strand
     * @param phase The phase
     */
	Feature(std::string &&seq, std::string &&source, std::string &&type, size_t start, size_t end, double score, char strand, short phase) : seqId_(std::forward<std::string>(seq)), source_(std::forward<std::string>(source)), type_(std::forward<std::string>(type)), start_(start),	end_(end), score_(score),strand_(strand), phase_(phase)
	{}

	/**
	 * \brief Assignment operator
	 * @param obj The object to assign
	 * @return Reference to the object
	 */
    Feature& operator=(const Feature&) = default;

    /**
     * \brief Move assignment operator
     * @param obj The obj to assign
     * @return Reference to the object
     */
    Feature& operator=(Feature&&) = default;

    /**
     * Standard destructor
     */
	virtual ~Feature() = default;

	/**
	 * Returns the name of an attribute
	 * @param key The name of the attribute
	 * @return The value
	 */
	const std::string&
	attribute(const std::string &key) const
	{
		return attributes_.at(key);
	}

	/**
	 * Sets the value of an attribute
	 * @param key The name of the attribute
	 * @param value The value of the attribute
	 */
	void
	attribute(const std::string &key, const std::string &value)
	{
		attributes_[key] = value;
	}


/*
    bool operator<=(const Feature &other) const
    {
    	if (other.seqId_ == seqId_)
    	{
			if ((!other.parents_.empty()) && (other.parents_.at(0) == attributes_.at("ID")))
				return true;
    		if (other.start_ == start_)
    			return other.end_ >= end_;
    		else
    			return other.start_ < start_;
    	}
    	else
    		return (other.seqId_ < seqId_);
    }




    bool operator>=(const Feature &other) const
    {
    	if (other.seqId_ == seqId_)
    	{
    		if (other.start_ == start_)
    			return other.end_ <= end_;
    		else
    			return other.start_ > start_;
    	}
    	else
    		return other.seqId_ > seqId_;
    }
*/
    /**
     * \brief equal comparison operator
     * @param other The object to compare to
     * @return true if equal else false
     * \warning Only the fields seqId, start_, end_ and type_ are compared
     */
	bool
	operator==(const Feature &other) const
	{

		if (other.seqId_ != seqId_)
			return false;
		if (other.start_ != start_)
			return false;
		if (other.end_ != end_)
			return false;
		if (other.type_ != type_)
			return false;
		return true;
	}

	/**
	 * \brief unqual comparison operator
	 * @param other The object to compare to
	 * @return true if unequal else false
	 * \warning Only the fields seqId, start_, end_ and type_ are compared
	 */
	bool
	operator!=(const Feature &other) const
	{

		if (other.seqId_ != seqId_)
			return true;
		if (other.start_ != start_)
			return true;
		if (other.end_ != end_)
			return true;
		if (other.type_ != type_)
			return true;
		return false;
	}

	friend bool operator<(const Feature &a, const Feature &b);
	friend bool operator>(const Feature &a, const Feature &b);

};

/**
 * \brief Simple function to read one type of feature.
 * If ignoreStrand is set the strand of a feature is ignored. This increases the likelihood of mixing up
 * features of different units.
 * @param inFile The file to read from.
 * @param featType The feature type to read.
 * @param ignoreStrand If true strands are ignored, otherwise not
 * @return A map from seqname -> features. Features are sorted by start position and size. When the strand is important,
 * '+' or '-' is added to the sequence name.
 */
std::map<std::string, std::vector<Feature> >
simpleGffRead(const std::string &inFile, const std::string &featType, bool ignoreStrand=false);

bool
operator<(const Feature &a, const Feature &b);

} /* namespace BioSeqDataLib */

#endif /* FEATURE_H_ */
