/*
 * DomainArrangement.hpp
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
 * \file DomainArrangement.hpp
 * \brief Header containing the DomainArrangment class.
 */
#ifndef DOMAINARRANGEMENT_HPP_
#define DOMAINARRANGEMENT_HPP_

// C++ header
#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <vector>
#include <iostream>
#include <math.h>

#include "PfamDomain.hpp"

// BioSeqDataLib header
#include "../utility/DSM.hpp"
#include "DomainExt.hpp"

namespace BioSeqDataLib
{

/** \addtogroup DomainGroup
 *  @{
 */


/**
 * \brief Class to store a domain arrangement.
 *
 * \tparam DomainType The domain type to use (e.g. PfamDomain)
 * \details A domain Arrangement is a set of domains in a certain order. A domain can occur several times in a row.
 */
template<typename DomainType>
class DomainArrangement
{
private:
	std::vector<DomainType> domains_;
	std::map<unsigned int, DomainType> collapsedDoms_;
	using Domain_t = std::vector<DomainType>;

public:
	using iterator = typename Domain_t::iterator;
	using const_iterator = typename Domain_t::const_iterator;
	using reverse_iterator = typename Domain_t::reverse_iterator;
	using const_reverse_iterator = typename Domain_t::const_reverse_iterator;

	typedef DomainType value_type;

	/**
	 * \brief Standard Constructor
	 */
	DomainArrangement():domains_()
	{}

	/**
	 * \brief Copy constructor
	 * @param Object to copy
	 */
	DomainArrangement(const DomainArrangement<DomainType>&) = default;


	/**
	 * \brief Move constructor
	 * @param rvalue of a DomamainArrangment
	 */
	DomainArrangement(DomainArrangement<DomainType> &&) = default;

	/**
	 * \brief Standard destructor
	 */
	virtual ~DomainArrangement() = default;

	/***********************************************************
	 *                       Operators                         *
	 ***********************************************************/

	/**
	 * \brief Copy assignment operator
	 * @param Object to copy
	 * @return Reference to the new object
	 */
	DomainArrangement& operator=(const DomainArrangement<DomainType>&) = default;

	/**
	 * \brief Move assignment operator
	 * @param The value to add.
	 * @return Reference to object.
	 */
	DomainArrangement& operator=(DomainArrangement<DomainType>&&) = default;


	/**
	 * \brief Gives access to a specific domain.
	 * @param index The index of the domain to access.
	 * @return A reference to the domain.
	 */
	DomainType
	&operator[](unsigned int index)
	{
		return domains_[index];
	}

	/**
	 * \brief Gives access to a specific domain.
	 * @param index The index of the domain to access.
	 * @return A const reference to the domain.
	 */
	const DomainType
	&operator[](unsigned int index) const
	{
		return domains_[index];
	}

	/**
	 * \brief Equality operator according to accession number.
	 * @param a1 First arrangement
	 * @param a2 second arrangement
	 * @return true if, equal else false
	 */
	friend bool operator== (const DomainArrangement<DomainType> &a1, const DomainArrangement<DomainType> &a2)
	{
		if (a1.size() != a2.size())
			return false;
		size_t len = a1.size();
		for (size_t i=0; i<len; ++i)
		{
			if (a1[i].accession() != a2[i].accession())
				return false;
		}
		return true;
	}

	/**
	 * \brief Inequality operator according to accession number.
	 * @param a1 First arrangement
	 * @param a2 second arrangement
	 * @return true if unequal, else false
	 */
	friend bool operator!= (const DomainArrangement<DomainType> &a1, const DomainArrangement<DomainType> &a2)
	{
		return !(a1==a2);
	}

	/**
	 * \brief Less operator using accession numbers.
	 * @param a1 First arrangement
	 * @param a2 Second arrangement
	 * @return true if a1 is smaller than a2 else false
	 */
	friend bool operator< (const DomainArrangement<DomainType> &a1, const DomainArrangement<DomainType> &a2)
	{

		size_t len = (a1.size() < a2.size()) ? a1.size() : a2.size();
		for (size_t i=0; i<len; ++i)
		{
			if (a1[i].accession() > a2[i].accession())
				return false;
			if (a1[i].accession() < a2[i].accession())
				return true;
		}
		if (a1.size() < a2.size())
			return true;
		else
			return false;
	}

	/**
	 * \brief Greater or equal operator using accession numbers.
	 * @param a1 First arrangement
	 * @param a2 Second arrangement
	 * @return true if a1 larger or equal to a2 else false
	 */
	friend bool operator>= (const DomainArrangement<DomainType> &a1, const DomainArrangement<DomainType> &a2)
	{
		return !(a1<a2);
	}


	/**
	 * \brief Greater operator using accession numbers.
	 * @param a1 First arrangement
	 * @param a2 Second arrangement
	 * @return true if a1 is smaller than a2 else false
	 */
	friend bool operator> (const DomainArrangement<DomainType> &a1, const DomainArrangement<DomainType> &a2)
	{

		size_t len = (a1.size() < a2.size()) ? a1.size() : a2.size();
		for (size_t i=0; i<len; ++i)
		{
			if (a1[i].accession() < a2[i].accession())
				return false;
			if (a1[i].accession() > a2[i].accession())
				return true;
		}
		if (a1.size() > a2.size())
			return true;
		else
			return false;
	}

	/**
	 * \brief Smaller or equal operator using accession numbers.
	 * @param a1 First arrangement
	 * @param a2 Second arrangement
	 * @return true if a1 smaller or equal to a2 else false
	 */
	friend bool operator<= (const DomainArrangement<DomainType> &a1, const DomainArrangement<DomainType> &a2)
	{
		return !(a1>a2);
	}

	template<class... Args>
	void
	erase(Args&&... args)
	{
		domains_.erase(std::forward<Args>(args)...);
	}


	/**
	 * \brief Returns an iterator the first domain of the arrangement.
	 * @return The iterator.
	 */
	iterator
	begin()
	{
		return domains_.begin();
	}

	/**
	 * \brief Returns an iterator to the end of the arrangement.
	 * @return The iterator.
	 */
	iterator
	end()
	{
		return domains_.end();
	}


	/**
	 * \brief Returns an iterator the first domain of the arrangement.
	 * @return The iterator.
	 */
	const_iterator
	begin() const
	{
		return domains_.begin();
	}

	/**
	 * \brief Returns an iterator to the end of the arrangement.
	 * @return The iterator.
	 */
	const_iterator
	end() const
	{
		return domains_.end();
	}

	reverse_iterator
	rbegin()
	{
		return domains_.rbegin();
	}


	reverse_iterator
	rend()
	{
		return domains_.rend();
	}


	const_reverse_iterator
	rbegin() const
	{
		return domains_.rbegin();
	}


	const_reverse_iterator
	rend() const
	{
		return domains_.rend();
	}

	/***********************************************************
	 *                       Functions                         *
	 ***********************************************************/
	/**
	 * \brief Adds a domain to the end of the arrangement.
	 * @param domain The domain to add.
	 */
	void
	push_back(DomainType &&domain)
	{
		domains_.push_back(std::forward<DomainType>(domain));
	}

	void
	push_back(DomainType &domain)
	{
		domains_.push_back(domain);
	}

	/**
	 * \brief Directly constructs a domain at the end of the arrangement.
	 * @param args Arguments needed by the constructor
	 */
	template<class... Args>
	void
	emplace_back(Args&&... args)
	{
		domains_.emplace_back(std::forward<Args>(args)...);
	}


	/**
	 * Inserts a Domain at the given postion
	 * @param  args The domain to Insert
	 * @return      Iterator to the inserted element
	 */
	template<class... Args>
	auto
	insert(Args&&... args) -> decltype (domains_.insert(std::forward<Args>(args)...))
	{
		return domains_.insert(std::forward<Args>(args)...);
	}


	/**
	 * Inserts a Domain at the given postion
	 * @param postion The position where the doamin should be inserted
	 * @param  args The domain to Insert
	 * @return      Iterator to the inserted element
	 */
	template <typename I, class... Args>
	auto emplace (I position, Args&&... args) -> decltype (domains_.emplace(position, std::forward<Args>(args)...))
	{
		return domains_.emplace(position, std::forward<Args>(args)...);
	}


	/**
	 * \brief Returns the number of domains in the arrangement.
	 * @return Number of domains.
	 */
	size_t
	size() const
	{
		return domains_.size();
	}


	/**
	 * \brief Returns a string of the domains.
	 * \param sep The separator to use.
	 * @return The string.
	 */
	std::string
	str(const char sep='-') const
	{
		std::string str = "";
		auto it = domains_.begin();
		str += it->accession();
		++it;
		for (auto itEnd=domains_.end(); it!=itEnd; ++it)
		{
			str.push_back(sep);
			str += it->accession();
		}
		return str;
	}


	/**
	 * \brief Removes all consecutive repeats of the same domain and only keeps the first one.
	 * \pre Domains are sorted by starting position
	 * @param reversable If set to true, domains are stored and can be reconstructed using the reconstruct function.
	 */
	void
	collapse(bool reversable = false)
	{
		if (domains_.empty())
			return;
		size_t nDomains = domains_.size();
		size_t pos=0;
		if (!reversable)
		{
			for (size_t i=1; i<nDomains; ++i)
			{
				if (domains_[pos].accession() != domains_[i].accession())
				{
					if (++pos != i)
						domains_[pos] = std::move(domains_[i]);
				}
			}
			domains_.resize(++pos);
		}
		else
		{
			for (size_t i=1; i<nDomains; ++i)
			{
				if (domains_[pos].accession() != domains_[i].accession())
					domains_[++pos] = domains_[i];
				else
					collapsedDoms_.emplace(i, std::move(domains_[i]));
			}
			domains_.resize(++pos);
		}
	}

	/**
	 * \brief Reconstructs the original domain arrangement.
	 * This works only if the collapse function has been called with reversable set to true.
	 */
	void
	reconstruct()
	{
		if (collapsedDoms_.empty())
			return;
		size_t oldSize = domains_.size();
		domains_.resize(domains_.size() + collapsedDoms_.size());
		size_t nDomains = domains_.size();
		auto it = collapsedDoms_.rbegin();
		auto itEnd = collapsedDoms_.rend();
		for (size_t i=nDomains; i>0; --i)
		{
			if (it != itEnd)
			{
				if ((i-1) != it->first)
					domains_[i-1] = std::move(domains_[--oldSize]);
				else
				{
					domains_[i-1] = std::move(it->second);
					++it;
				}
			}
			else
				break;
		}

	}

	/**
	 * \brief Removes all domains from the arrangement.
	 */
	void
	clear()
	{
		domains_.clear();
	}

	/**
	 * \brief Checks if the DA is empty.
	 * @return true if empty else false
	 */
	bool
	empty() const
	{
		return domains_.empty();
	}

	/**
	 * \brief Clear arrangement from all overlaps that are larger than the given threshold.
	 * @param dbImportance The order of importance for the domain databases. Domains belonging to none of the databases listed will be removed.
	 * @param maxAbsOverlap The maximum overlap in terms of number of amino acids between two domains that is allowed.
	 * @param maxFracOverlap The maximum fraction of a domain that can be overlapping another one. Since an overlap can happen on both sides the maximal overlap can be up to twice as much.
	 */
	void
	solveDbOverlaps(std::vector<DomainDB> dbImportance, size_t maxAbsOverlap, float maxFracOverlap)
	{
		std::set<DomainDB> usedDbs;
		for (auto db: dbImportance)
			usedDbs.insert(db);

		size_t nDomains = domains_.size();
		std::set<size_t> toDelete;
		for (size_t i=0; i<nDomains; ++i)
		{
			if (usedDbs.count(domains_[i].db()) == 0)
				toDelete.insert(i);
		}

		for (auto it = toDelete.rbegin(); it != toDelete.rend(); ++it)
			domains_.erase(domains_.begin() + (*it));
		toDelete.clear();

		for (DomainDB db : dbImportance)
		{

			bool changed = true;
			while (changed)
			{
				nDomains = domains_.size();
				for (size_t i=1; i<nDomains; ++i)
				{
					DomainType &dom1 = domains_[i-1];
					DomainType &dom2 = domains_[i];
					if (dom1.end() > dom2.start())
					{
						float overlap = dom1.end() - dom2.start() + 1;
						if (((dom1.db() == db) || (dom2.db() == db)) && ((overlap > maxAbsOverlap) || ((overlap/dom1.length()) > maxFracOverlap ) || ((overlap/dom2.length()) > maxFracOverlap ) ))
						{
							if (dom1.db() == dom2.db())
							{
								if (dom1.evalue() < dom2.evalue())
									toDelete.insert(i);
								else
									toDelete.insert(i-1);
							}
							else
								toDelete.insert(dom1.db() == db ? i : i-1);
						}
					}
				}
				for (auto it = toDelete.rbegin(); it != toDelete.rend(); ++it)
					domains_.erase(domains_.begin() + (*it));
				if (toDelete.empty())
					changed = false;
				toDelete.clear();
			}
		}
	}



};


/**
 * \brief Wrapper to be used in a Sequence
 * \details This wrapper is made to be used together with a Sequence.
 */
template<typename DomainType>
class DomainArrangementWrapper
{
private:
	DomainArrangement<DomainType> arrangement_;

public:

	typedef DomainType value_type;
	/**
	 * \brief Default constructor
	 */
	DomainArrangementWrapper():arrangement_()
	{}

	// No copy/move constructors
	/// @cond HIDDEN
	DomainArrangementWrapper(const DomainArrangementWrapper<DomainType>&) = delete;
	/// @endcond

	/**
	 * \brief The move constructor
	 * @param The object to move.
	 */
	DomainArrangementWrapper(DomainArrangementWrapper<DomainType>&&) = default;

	/**
	 * \brief Default destructor
	 */
	virtual
	~DomainArrangementWrapper()
	{}

	// No assignment operators
	/// @cond HIDDEN
	DomainArrangementWrapper& operator=(const DomainArrangementWrapper<DomainType>&)  = delete;
	DomainArrangementWrapper& operator=(DomainArrangementWrapper<DomainType> &&other)  = delete;
	/// @endcond

	/**
	 * \brief Returns the stored domain arrangement.
	 * @return The domain arrangement.
	 */
	DomainArrangement<DomainType> &
	domainArrangement()
	{
		return arrangement_;
	}

	/**
	 * \brief Sets the domain arrangement.
	 * @param arrangement The arrangement to store.
	 */
	void
	domainArrangement(DomainArrangement<DomainType> &&arrangement)
	{
		arrangement_=std::forward<DomainArrangement<DomainType> >(arrangement);
	}

};



/**
 *
 * \brief Counts the occurrence of each clan an a domain set.
 * \param useDomain If useDomain is true, the domain accession number is used if the domain is not assigned to any clan.
 * @return The domain and the number of their occurrence.
 * \relates DomainArrangement
 */
std::map<std::string, int>
clanCounts(const DomainArrangement<PfamDomain> &daSet, bool useDomain=true);


/**
 * \brief Counts the occurrence of each domain an a domain set.
 * @return The domain and the number of their occurrence.
 * \tparam DomainType The domain type
 * \relates DomainArrangement
 */
template<typename DomainType>
std::map<std::string, int>
domainCounts(const DomainArrangement<DomainType> &daSet)
{
	std::map<std::string, int> summary;
	for (const DomainType &domain : daSet)
		summary[domain.accession()] += 1;
	return summary;
}

/**
 * \brief Calculates the cos distance between two domain arrangements.
 * \param da1 The first domain arrangement.
 * \param da2 The second domain arrangement.
 * @return The cos score.
 * \tparam DomainType The domain type
 *
 * \relates DomainArrangement
 */
template<typename DomainType>
float
cos(const DomainArrangement<DomainType> &da1, const DomainArrangement<DomainType> &da2, const DSM &dsm)
{
	// create common universe
	std::set<std::string> universe;
	for (const DomainType &domain : da1)
		universe.insert(domain.accession());
	for (const DomainType &domain : da2)
		universe.insert(domain.accession());

	// create vectors
	std::vector<float> vec1, vec2;
	for (const std::string &dom : universe)
	{
		float score = 0;
		float val;
		for (const DomainType &domain : da1)
		{
			val = dsm.val(domain.accession(), dom);
			if (val > score)
				score =val;
		}
		vec1.push_back(score);
		score = 0;
		for (const DomainType &domain : da2)
		{
			val = dsm.val(domain.accession(), dom);
			if (val > score)
				score =val;
		}
		vec2.push_back(score);
	}

	// calculate norm
	size_t len = universe.size();
	float score= 0;
	float nVec1 = 0;
	float nVec2 = 0;
	for (size_t i=0; i<len; ++i)
	{
		score += vec1[i] * vec2[i];
		nVec1 += vec1[i] * vec1[i];
		nVec2 += vec2[i] * vec2[i];
	}

	return score/(sqrt(nVec1) * sqrt(nVec2));

}





/** @} */ // Domain module




/**
 *
 *\brief Function to reconstruct Split Domains into whole Domains
 *@param DomainArrangement The DomainArrangement
 *@param ot The overlap treshhold for domains
 *@param maxdist The maximal distance between two parts of a split domains.
 */
template<typename DomainType>
void
splitDomRec(DomainArrangement<DomainType> &domainArrangement, unsigned int ot, unsigned int maxdist)
{
	unsigned int k,l,n,q;
	std::string domaID;
	std::string nextdomaID;
	unsigned int maxhmmot=10;
	auto& splitdoma = domainArrangement;
	for (k=0;k<splitdoma.size(); ++k)
	{
		//Code to prevent segmentatation faults
		unsigned int i=splitdoma.size();


		if (i>1)
		{
			domaID = splitdoma[k].accession();
			if ((k+1)<i)
			{
				unsigned int m=0;
				nextdomaID = splitdoma[k+1].accession();
			//Checks if the next Domain has the same name as the current
				if (domaID == nextdomaID)
					m=1;
				//Continues to check how many domains of the same name are in a row
				while ((m>0) & (domaID == nextdomaID) & ((k+m+1)<splitdoma.size()))
				{
					nextdomaID = splitdoma[k+m+1].accession();
					if (domaID == nextdomaID)
						++m;
				}
				if (m>0)
				{
					size_t splitdomaHMM = splitdoma[k].hmmLength();
					size_t splitdomaHMMStart = splitdoma[k].hmmStart();
					size_t splitdomaHMMEnd = splitdoma[k].hmmEnd();
					size_t splitdomaHMMref = splitdomaHMMEnd - splitdomaHMMStart;
					size_t splitdomaHMMLength = splitdomaHMMEnd - splitdomaHMMStart;
					double crit =splitdomaHMM;
					//While the Domain length is less than the total HMMAllignment*minimal Differece adds the length of the next Domain to the total length
					l=0;
					unsigned int v=0;
					while ((splitdomaHMMLength <= crit) & (l<m) & (v!=1))
					{
						for (l=0; l<m; ++l)
						{
							size_t nextsplitdomaHMMStart = splitdoma[k+l+1].hmmStart();
							size_t nextsplitdomaHMMEnd = splitdoma[k+l+1].hmmEnd();
							size_t nextsplitdomaHMMLength = nextsplitdomaHMMEnd - nextsplitdomaHMMStart;
							size_t currentsplitdomaHMMEnd=splitdoma[k+l].hmmEnd();
							double hmmpercentoverlap=0.2*nextsplitdomaHMMLength;
							if (((currentsplitdomaHMMEnd-hmmpercentoverlap)<nextsplitdomaHMMStart) && (currentsplitdomaHMMEnd-maxhmmot<nextsplitdomaHMMStart))
								splitdomaHMMLength += nextsplitdomaHMMLength;
							else
							{
								v=1;
								break;
							}
						}
					}
					//If the Domain is now longer than it used to be but no more than 10% longer than the HMM allignment checks for continuity, overlap and max distance
					if ((splitdomaHMMref<splitdomaHMMLength) & (splitdomaHMMLength<=crit*1.1))
					{
						//If all criteria are met (Domain is longer, overlap of start and end no larger than specified and distance smaller than max distance, changes
						for (n=0; n<l; ++n)
						{
							unsigned long splitdomaEnd=splitdoma[k+n].end();
							unsigned long nextsplitdomaStart=splitdoma[k+n+1].start();
							if (((splitdomaEnd-ot)<nextsplitdomaStart) & ((splitdomaEnd+maxdist)>nextsplitdomaStart))
							{

								unsigned long ENDENV = splitdoma[k+n+1].envEnd();
								size_t ENDHMM = splitdoma[k+n+1].hmmEnd();
								size_t ENDSEQ = splitdoma[k+n+1].end();
								splitdoma[k].envEnd(ENDENV);
								splitdoma[k].hmmEnd(ENDHMM);
								splitdoma[k].end(ENDSEQ);
								double splitdomaEvalue=splitdoma[k+n].evalue();
								double nextsplitdomaEvalue=splitdoma[k+n+1].evalue();
								if (nextsplitdomaEvalue>splitdomaEvalue)
									splitdoma[k].evalue(nextsplitdomaEvalue);
							}
							else
								break;
						}
						//Deletes the surplus domains
						for (q=0; q<n; ++q)
							splitdoma.erase(splitdoma.begin()+(k+q+1));
					}

				}
			}
		}
	}

}

} /* namespace BioSeqDataLib */

#endif /* DOMAINARRANGEMENT_HPP_ */
