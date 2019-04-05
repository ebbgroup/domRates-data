/*
 * PfamDomain.hpp
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




/**
 * \file PfamDomain.hpp
 * \brief Header containing the PfamDomain class.
 */
#ifndef SRC_ANNOTATION_PFAMDOMAIN_HPP_
#define SRC_ANNOTATION_PFAMDOMAIN_HPP_

#include <map>

#include "DomainExt.hpp"


namespace BioSeqDataLib
{

/** \addtogroup DomainGroup
 *  @{
 */

/**
 * \brief Class representing the Pfam domain.
 *
 * This domain inherits all fields from the Domain class and has additional fields specific for Pfam domains (e.g. clan).
 */
class PfamDomain : public DomainExt
{
private:

	// PfamDomain
	double significance_;	/// The significance
	std::string clan_;		/// The clan of the domain
	std::string type_;		/// The domain type

public:
	/**
	 * \brief Standard constructor
	 */
	PfamDomain()=default;

	/**
	 * \brief PfamDomain constructor.
	 * @param name The PfamDomain name.
	 * @param accession The accession number.
	 * @param seqStart The first position in the sequence.
	 * @param seqEnd The last position in the sequence.
	 * @param envStart The first position of the envelope.
	 * @param envEnd The last position of the envelope.
	 * @param hmmStart The first position in the HMM.
	 * @param hmmEnd The last position in the HMM.
	 * @param hmm_length The length of the HMM.
	 * @param bit_score The bitscore.
	 * @param evalue The e-value of the PfamDomain.
	 * @param significance The significance value.
	 * @param clan The clanName.
	 * @param type The type of the PfamDomain.
	 * @param version The version of the PfamDomain.
	 * @param id The id of the PfamDomain.
	 */
	PfamDomain(std::string accession, std::string name, size_t seqStart, size_t seqEnd, size_t envStart, size_t envEnd, size_t hmmStart, size_t hmmEnd, size_t hmm_length, double bit_score, double evalue, double significance, std::string clan, std::string type) : DomainExt(accession, name, seqStart, seqEnd, envStart, envEnd, hmmStart, hmmEnd, hmm_length, bit_score, evalue), significance_(significance), clan_(clan), type_(type)
	{}


	/**
	 * \brief Constructor initializing all values of the parent lass domain.
	 * All PfamDomain specific fields are set to either 0 or the empty string.
	 * @param accession The accession number
	 * @param name The name of the domain
	 * @param seqStart The start of the domain in the sequence (start at 0)
	 * @param seqEnd The last position of the domain in the sequence (start at 0)
	 * @param envStart Start of the envelope
	 * @param envEnd End of the envelope
	 * @param hmmStart Start of the match in the HMM
	 * @param hmmEnd End of the match in the HMM
	 * @param hmm_length Length of the HMM
	 * @param bit_score The bit score
	 * @param evalue The e-value of the domain
	 * @param id ID - default value: -1
	 */
	PfamDomain(const std::string &accession,
			const std::string &name,
			unsigned long seqStart,
			unsigned long seqEnd,
			unsigned long envStart,
			unsigned long envEnd,
			unsigned long hmmStart,
			unsigned long hmmEnd,
			unsigned long hmm_length,
			double bit_score,
			double evalue)
			: DomainExt(accession, name, seqStart, seqEnd, envStart, envEnd, hmmStart, hmmEnd, hmm_length, bit_score, evalue), significance_(0), clan_(""), type_("")
	{}



	/**
	 * \brief PfamDomain constructor
	 * @param accession The accession number.
	 */
	//explicit PfamDomain(const std::string &accession) : DomainExt(accession, 0, 0), significance_(0), clan_("NA"), type_("")
	//{}

	PfamDomain(std::string accession, unsigned long start, unsigned long end, double evalue) : DomainExt(accession, start , end, evalue), significance_(0), clan_(""), type_("")
	{}


	/**
	 * \brief Move constructor.
	 * @param PfamDomain to move.
	 */
	PfamDomain(PfamDomain &&rval) : DomainExt(std::move(rval))
	{
		significance_ = std::move(rval.significance_);
		clan_ = std::move(rval.clan_);
		type_ = std::move(rval.type_);
	}

	/**
	 * \brief Copy constructor.
	 * @param The PfamDomain to copy.
	 */
	PfamDomain(const PfamDomain &) = default;

	/**
	 * \brief Standard destructor.
	 */
	virtual ~PfamDomain();

	/**
	 * \brief Move assignment operator.
	 * @param The PfamDomain to move.
	 * @return Reference to the PfamDomain.
	 */
	PfamDomain& operator=(PfamDomain&&) = default;

	/**
	 * \brief Copy assignment operator.
	 * @param The PfamDomain to copy.
	 * @return Reference to the PfamDomain.
	 */
	PfamDomain& operator=(const PfamDomain&) = default;


	/**
	 * \brief Returns the clan name.
	 * \return The sequence clan name.
	 */
	std::string
	clan() const
	{
		return clan_;
	}

	/**
	 * \brief Sets a new value for the clan.
	 * \param new_val The new clan name.
	 */
	void
	clan(std::string &&new_val)
	{
		clan_ = std::forward<std::string>(new_val);
	}

	/**
	 * \brief Returns the significance.
	 * \return The significance of the domain.
	 */
	double
	significance() const
	{
		return significance_;
	}

	/**
	 * \brief Sets a new value for significance.
	 * \param new_val The new significance.
	 */
	void
	significance(double new_val)
	{
		significance_ = new_val;
	}

	/**
	 * \brief Returns the type of the domain
	 * @return The domain type
	 */
	std::string
	type() const
	{
		return type_;
	}

	/**
	 * \brief Sets the type of the domain
	 * @param The domain type
	 */
	void
	type(std::string &&new_val)
	{
		type_ = std::forward<std::string>(new_val);
	}

};









/** @} */ // Domain module



} //namespace BioSeqDataLib


#endif /* SRC_ANNOTATION_PFAMDOMAIN_HPP_ */
