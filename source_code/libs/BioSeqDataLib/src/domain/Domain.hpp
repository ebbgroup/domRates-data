/*
 * Domain.hpp
 *
 *  Created on: Mar 14, 2016
 *      Author: ckeme_01
 * 	 Copyright: 2016
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
 * \file Domain.hpp
 * \brief Header containing the Domain class.
 */
#ifndef SRC_DOMAIN_HPP_
#define SRC_DOMAIN_HPP_

#include <memory>
#include <string>
#include <algorithm> 

#include "../utility/stringHelpers.hpp"


namespace BioSeqDataLib
{

/**
 * \brief A class to encode the databases of a domain.
 */
enum class DomainDB {unknown, pfam, superfamily, gene3d, smart, tigrfams, prosite, panther, hamap, pirsf, prints, prodom};


/** \addtogroup DomainGroup
 *  @{
 */

/**
 * \brief The most basic Domain class.
 */
class Domain
{
private:
	std::string accession_;
	unsigned long start_;
	unsigned long end_;
	double evalue_;
	DomainDB domDB_;

public:

	/**
	 * \brief Standard constructor.
	 */
	Domain();

	/**
	 * \brief Constructor setting values.
	 * @param accession The accession name.
	 * @param start The start position.
	 * @param end The end position.
	 * @param evalue The evalue.
	 */
	Domain(const std::string &accession, unsigned long start, unsigned long end, double evalue, DomainDB db = DomainDB::unknown);

	/**
	 * \overload
	 */
	Domain(std::string &&accession, unsigned long start, unsigned long end, double evalue, DomainDB db = DomainDB::unknown);

	/**
	 * \brief Copy constructor
	 * @param  Domain to copy
	 */
	Domain(const Domain &) = default;


	/**
	 * \brief Move constructor
	 * @param  Domain to move
	 */
	Domain(Domain &&) = default;

	/**
	 * \brief Standard destructor.
	 */
	virtual ~Domain();

	/**
	 * \brief Move assignment operator.
	 * @param The PfamDomain to move.
	 * @return Reference to the PfamDomain.
	 */
	Domain& operator=(Domain&&) = default;

	/**
	 * \brief Copy assignment operator.
	 * @param The PfamDomain to copy.
	 * @return Reference to the PfamDomain.
	 */
	Domain& operator=(const Domain&) = default;


	/**
	 * \brief comparison operator
	 * @param b Domain to compare to
	 * @result true if domain has a smaller start value than b
	 */
	bool operator<(const Domain& b)  const
	{
	    return (this->start_ < b.start_);
	}

	/**
	 * \brief comparison operator
	 * @param b Domain to compare to
	 * @result true if domain has a larger start value than b
	 */
	bool operator>(const Domain& b)  const
	{
	    return (this->start_ > b.start_);
	}

	/**
	 * \brief Returns the accession number
	 * \return The accession number
	 */
	std::string
	accession() const
	{
		return accession_;
	}

	/**
	 * \brief Sets a new value for the accession number.
	 * \param new_val The accession number.
	 */
	void
	accession(std::string &&new_val)
	{
		accession_ = std::forward<std::string>(new_val);
	}

	void
	accession(const std::string &new_val)
	{
		accession_ = new_val;
	}

	/**
	 * \brief Returns the domain start in the sequence.
	 * \return The first position of the domain.
	 */
	unsigned long
	start() const
	{
		return start_;
	}

	/**
	 * \brief Sets the PfamDomain start in the sequence
	 * \param new_val The first position of the domain.
	 */
	void
	start(size_t new_val)
	{
		start_ = new_val;
	}

	/**
	 * \brief Returns the last position of the Domain in the sequence.
	 * \return The last position of the Domain.
	 */
	unsigned long
	end() const
	{
		return end_;
	}

	/**
	 * \brief Sets a new value for the last position of the domain.
	 * \param new_val The last position of the domain.
	 */
	void
	end(size_t new_val)
	{
		end_ = new_val;
	}


	/**
	 * \brief Returns the length of the domain.
	 * @return The length of the domain.
	 */
	unsigned long
	length() const
	{
		return end_-start_+1;
	}

	/**
	 * \brief Returns the evalue.
	 * \return The evalue of the domain.
	 */
	double
	evalue() const
	{
		return evalue_;
	}

	/**
	 * \brief Sets a new value for evalue.
	 * \param new_val The new evalue.
	 */
	void
	evalue(double new_val)
	{
		evalue_ = new_val;
	}

	/**
	 * \brief Returns the database the domain belongs to.
	 * @return The database
	 */
	DomainDB
	db() const
	{
		return domDB_;
	}

	/**
	 * \brief Sets the database field.
	 * @param db The databaes. Has to be of type DomainDB
	 */
	void
	db(DomainDB db)
	{
		domDB_ = db;
	}


	/**
	 * \brief Calculates the distance (gap) between two domains
	 * return if positive the length of the gap between two domains if negative the number of overlapping positions.
	 */
	int
	distance_overlap(const Domain &dom)
	{
		auto start = std::max(dom.start(), this->start());
		auto end = std::min(dom.end(), this->end());
		int value = end - start + 1;
		return value *-1;
	}

};


/** @} */ // Domain module

} /* namespace BioSeqDataLib */

#endif /* SRC_DOMAIN_HPP_ */
