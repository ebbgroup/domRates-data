/*
 * AlignmentFunctions.hpp
 *
 *  Created on: 22 Oct 2013
 *      Author: CarstenK
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

#include "../utility/TwoValues.hpp"

/**
 * \file AlignmentFunctions.hpp
 * \brief Functions related to alignments analysis and modification
 */

#ifndef ALIGNMENTFUNCTIONS_HPP_
#define ALIGNMENTFUNCTIONS_HPP_

namespace BioSeqDataLib
{

/** @defgroup alignmentFunction Alignment functions
 *  Group with basic alignment functions
 *  @{
 */


/**
 * \brief Checks if a set of sequences is an alignment.
 *
 * This function checks if all sequences have the same length and if no columns consist of gaps only.
 * \tparam SequenceSet A SequenceSet or an alignment.
 * \param set The set to check.
 * \relates Alignment
 * \relates SequenceSet
 * \return A pair giving an error code and the position.
 *
 * Error codes:
 *  0 = everything is fine
 *  1 = length differ (the first sequence having a different length is reported)
 *  2 = column consists of gaps only (first column of gaps only is reported)
 *  3 = no sequences
 */
template<typename SequenceSetType>
TwoValues
alnCheck(const SequenceSetType &set)
{
	size_t nSeqs = set.size();
	if (nSeqs == 0)
		return TwoValues(3,0);
	size_t length = set[0].size();
	std::vector<size_t>gapCount(length, 0);
	size_t j;
	for (size_t i=0; i<nSeqs; ++i)
	{
		if (set[i].size() != length)
			return TwoValues(1,i);
		const typename SequenceSetType::value_type &seq = set[i];
		for (j=0; j<length; ++j)
			if (seq[j] == '-')
				++gapCount[j];
	}

	for (j=0; j<length; ++j)
	{
		if (gapCount[j] == nSeqs)
			return TwoValues(2,j);
	}


	return TwoValues(0,0);
}


/**
 * \brief Calculates the identity of a given alignment.
 * \param set The Alignment to analyse.
 */
template<typename SequenceType>
TwoValues
identity(Alignment<SequenceType> &aln)
{
	size_t nSeqs = aln.size();
	size_t length = aln[0].size();
	size_t same=0, total=0;
	size_t j;
	std::vector <std::map<char, size_t> >profile(length, std::map<char, size_t>());
	std::map<char, size_t>::iterator it, itEnd;
	for (size_t i=0; i<nSeqs; ++i)
	{
		Sequence<> &seq = aln[i];
		for (j=0; j<length; ++j)
		{
			if (seq[j] != '-')
				++profile[j][tolower(seq[j])];
		}
	}


	for (j=0; j<length; ++j)
	{
		itEnd = profile[j].end();
		size_t tmpTotal=0;
		for (it=profile[j].begin(); it != itEnd; ++it)
		{
			same += (it->second * (it->second-1))/2;
			tmpTotal += it->second;
		}
		total += (tmpTotal * (tmpTotal-1))/2;
	}

	return TwoValues(same, total);
}



/**
 * \brief Shrinks the alignment to a sub segment.
 * @param aln The sub alignment
 * @param first The first position
 * @param last The last position
 * @param seqName If given the sequence name identifies the sequence that should be used for the given positions.
 */
template<typename SequenceType>
void
shrinkToSubAln(Alignment<SequenceType> &aln, size_t first, size_t last, const std::string &seqName="")
{
	size_t nSeqs = aln.size();
	if (!seqName.empty())
	{
		size_t alnLength = aln.length();
		SequenceType &seq = aln[seqName];
		size_t counter=0;
		++first;
		++last;
		size_t tmpFirst=-1, tmpLast=-1;
		for (size_t i=0; i<alnLength; ++i)
		{
			if (seq[i] != '-')
			{
				++counter;
				if (counter == first)
					tmpFirst = i;
				if (counter == last)
				{
					tmpLast=i;
					break;
				}
			}
		}
		first=tmpFirst;
		last=tmpLast;
	}


	if (first != 0)
	{
		for (size_t i=0; i<nSeqs; ++i)
		{
			size_t pos=0;
			SequenceType &seq = aln[i];
			for (size_t j=first; j<=last; ++j)
				seq[pos++] = seq[j];
		}
	}
	for (size_t i=0; i<nSeqs; ++i)
		aln[i].resize(last-first+1);
}

/** @*/ //  Alignment functions

} //BioSeqDataLib





#endif /* ALIGNMENTFUNCTIONS_HPP_ */
