/*
 * SFDomain.hpp
 *
 *  Created on: 17 Aug 2015
 *      Author: Carsten Kemena
 *		 Email: c.kemena[@]uni-muenster.de
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



#ifndef SRC_ANNOTATION_SFDOMAIN_HPP_
#define SRC_ANNOTATION_SFDOMAIN_HPP_

#include "DomainExt.hpp"

namespace BioSeqDataLib
{


/**
 * \file SFDomain.hpp
 * \brief Header containing the Superfamily Domain class.
 */


/** \addtogroup DomainGroup
 *  @{
 */


/**
 * \brief The class to represent a SuperFamily domain.
 */
class SFDomain : public DomainExt {

private:
	std::string scopID_;

public:
	SFDomain();

	SFDomain(std::string accession, std::string name, size_t seqStart, size_t seqEnd, size_t envStart, size_t envEnd, size_t hmmStart, size_t hmmEnd, size_t hmm_length, double bit_score, double evalue) : DomainExt(accession, name, seqStart, seqEnd, envStart, envEnd, hmmStart, hmmEnd, hmm_length, bit_score, evalue)
	{}

	SFDomain(std::string accession, std::string name, size_t seqStart, size_t seqEnd, size_t envStart, size_t envEnd, size_t hmmStart, size_t hmmEnd, size_t hmm_length, double bit_score, double evalue, const std::string &scopID) : DomainExt(accession, name, seqStart, seqEnd, envStart, envEnd, hmmStart, hmmEnd, hmm_length, bit_score, evalue), scopID_(scopID)
	{}

	SFDomain(std::string accession, size_t start, size_t end, double evalue) : DomainExt(accession, start , end, evalue), scopID_("")
	{}


	virtual ~SFDomain();

	/**
	 * The model
	 * @return
	 */
	std::string
	modelID() const
	{
		return this->name();
	}

	/**
	 *
	 * @return
	 */
	std::string
	scopID() const
	{
		return scopID_;
	}

	/**
	 *
	 * @return
	 */
	std::string
	scopDomID() const
	{
		return this->accession();
	}
};

/** @} */ // Domain module




} // namespace BioSeqDataLib

#endif /* SRC_ANNOTATION_SFDOMAIN_HPP_ */
