/*
 * Domain.hpp
 *
 *  Created on: 15 Nov 2013
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
 * \file DomainExt.hpp
 * \brief Header containing the DomainExt class.
 */
#ifndef DomainExt_HPP_
#define DomainExt_HPP_

#include <string>

#include "Domain.hpp"

namespace BioSeqDataLib
{

/** \addtogroup DomainGroup
 *  @{
 */

/**
 * \brief Class to represent a single domain. It contains the most important fields available from a hmmscan output.
 */

class DomainExt : public Domain
{

private:
	std::string name_; // The name of the domain
	unsigned long envStart_; // Start of the envelope
	unsigned long envEnd_;	// End of the envelope
	unsigned long hmmStart_;	// Start of the match in the HMM
	unsigned long hmmEnd_;	// End of the match in the HMM
	unsigned long hmmLength_; // Length of the HMM
	double bit_score_; // The bit score


public:
	/**
	 * \brief Standard constructor
	 */
	DomainExt() = default;

	/**
	 * \brief Constructor initializing all values of a domain
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
	DomainExt(const std::string &accession,
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
			: Domain(accession, seqStart, seqEnd, evalue), name_(name), envStart_(envStart), envEnd_(envEnd), hmmStart_(hmmStart), hmmEnd_(hmmEnd), hmmLength_(hmm_length), bit_score_(bit_score)
	{}

	/**
	 * \brief Simple constructor only setting the domain id.
	 * @param accession The accession number.
	 */
	//DomainExt(std::string accession) : accession_(accession), name_(""), seqStart_(0), seqEnd_(0), envStart_(0), envEnd_(0), hmmStart_(0), hmmEnd_(0), hmmLength_(0), bit_score_(0), evalue_(0), id_(-1)
	//{}




	DomainExt(std::string accession, unsigned long start, unsigned long end, double evalue) : Domain(accession, start, end, evalue), name_(""), envStart_(start), envEnd_(end), hmmStart_(0), hmmEnd_(0), hmmLength_(0), bit_score_(0)
	{}

	//DomainExt(std::string accession) : Domain(accession, 0, 0, -1.0), name_(""), envStart_(0), envEnd_(0), hmmStart_(0), hmmEnd_(0), hmmLength_(0), bit_score_(0)
	//{}

	/**
	 * \brief Copy constructor
	 * @param  Domain to copy
	 */
	DomainExt(const DomainExt &) = default;

	/**
	 * \brief Move constructor
	 * @param  Domain to move
	 */
	DomainExt(DomainExt &&) = default;

	/**
	 * \brief Destructor
	 */
	virtual ~DomainExt();

	/**
	 * \brief Move assignment operator.
	 * @param The PfamDomain to move.
	 * @return Reference to the PfamDomain.
	 */
	DomainExt& operator=(DomainExt&&) = default;

	/**
	 * \brief Copy assignment operator.
	 * @param The PfamDomain to copy.
	 * @return Reference to the PfamDomain.
	 */
	DomainExt& operator=(const DomainExt&) = default;




	/**
	 * \brief Returns the sequence name
	 * \return The sequence name.
	 */
	std::string
	name() const
	{
		return name_;
	}

	/**
	 * \brief Sets a new value for name
	 * \param new_val The new name.
	 */
	void
	name(std::string &&new_val)
	{
		name_ = std::forward<std::string>(new_val);
	}



	/**
	 * \brief Returns the envelope start in the sequence.
	 * \return The first position of the envelope.
	 */
	size_t
	envStart() const
	{
		return envStart_;
	}

	/**
	 * \brief Sets a new value for the envelope start.
	 * \param new_val The first position of the envelope.
	 */
	void
	envStart(size_t new_val)
	{
		envStart_ = new_val;
	}

	/**
	 * \brief Returns the last position of the envelope.
	 * \return The last position of the envelope.
	 */
	size_t
	envEnd() const
	{
		return envEnd_;
	}

	/**
	 * \brief Sets a new value for the envelope end.
	 * \param new_val The new envelope end.
	 */
	void
	envEnd(size_t new_val)
	{
		envEnd_ = new_val;
	}

	/**
	 * \brief Returns the length of the envelope.
	 * @return The length of the envelope.
	 */
	size_t
	envLength() const
	{
		return envEnd_-envStart_+1;
	}

	/**
	 * \brief Returns the first position of the hmm match.
	 * \return The hmm start position.
	 */
	size_t
	hmmStart() const
	{
		return hmmStart_;
	}

	/**
	 * \brief Sets a new value for hmmStart
	 * \param new_val The first position of the hmm.
	 */
	void
	hmmStart(size_t new_val)
	{
		hmmStart_ = new_val;
	}

	/**
	 * \brief Returns the hmm end position.
	 * \return The hmm end position.
	 */
	size_t
	hmmEnd() const
	{
		return hmmEnd_;
	}

	/**
	 * \brief Sets a new value for HMM end.
	 * \param new_val The new HMM end position.
	 */
	void
	hmmEnd(size_t new_val)
	{
		hmmEnd_ = new_val;
	}

	/**
	 * \brief Returns the length of the HMM.
	 * @return The length of the HMM.
	 */
	size_t
	hmmLength() const
	{
		return hmmLength_;
	}

	/**
	 * \brief Returns the bit score of the domain
	 * @return The bits score of the domain
	 */
	double
	bitscore() const
	{
		return bit_score_;
	}
};


/** @} */ // Domain module




} /* namespace BioSeqDataLib */

#endif /* DomainExt_HPP_ */
