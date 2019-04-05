/*
 * SeqSetFunctions.hpp
 *
 *  Created on: 14 Oct 2013
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

/**
 * \file SeqSetFunctions.hpp
 * \brief Collection of several function related to SequenceSets.
 */


// C header
#include <cstdlib>

// C++ header
#include <algorithm>
#include <string>
#include <utility>

// BioSeqDataLib header
#include "../align/AlignmentMatrix.hpp"

#ifndef SEQSETFUNCTIONS_HPP_
#define SEQSETFUNCTIONS_HPP_


namespace BioSeqDataLib
{


/**
 * \brief Returns the length of the shortest sequence.
 * @param set The SequenceSet to analyse.
 * @return The length of the shortest sequence.
 * \tparam SequenceSetType The type of the sequence to use.
 * \relates SequenceSet
 */
template<typename SequenceSetType>
size_t
minLength(const SequenceSetType &set)
{
	size_t nSeqs = set.size();
	size_t min = std::numeric_limits<size_t>::max();
	for (size_t i=0; i<nSeqs; ++i)
	{
		if (set[i].size() < min)
			min=set[i].size();
	}
	return min;
}

/**
 * \brief Returns the length of the longest sequence.
 * @param set The SequenceSet to analyse.
 * @return The length of the longest sequence.
 * \tparam SequenceSetType The type of the sequence to use.
 * \relates SequenceSet
 */
template<typename SequenceSetType>
size_t
maxLength(const SequenceSetType &set)
{
	size_t nSeqs = set.size();
	size_t max = 0;
	for (size_t i=0; i<nSeqs; ++i)
	{
		if (set[i].size() > max)
			max=set[i].size();
	}
	return max;
}

/**
 * \brief Returns the average sequence length of the set.
 * @param set The SequenceSet to analyse.
 * @return The average length.
 * \tparam SequenceSetType The type of the sequence to use.
 * \relates SequenceSet
 */
template<typename SequenceSetType>
double
avgLength(const SequenceSetType &set)
{
	size_t nSeqs = set.size();
	size_t avgLength = 0;
	for (size_t i=0; i<nSeqs; ++i)
		avgLength += set[i].size();
	return (1.0*avgLength/nSeqs);
}


/**
 * \brief Calculates the average identity of the sequence set.
 * \details The identity is calculated by producing alignments
 * between each pair of sequences and calculating the identity for each of them.
 * \tparam The SequenceSet type.
 * @param set The sequence set to analyse.
 * @param simMat The similarity matrix to use for the alignments.
 * @return The average identity of the sequence set.
 */
template<typename SequenceSetType>
float
avgId(const SequenceSetType &set, const SimilarityMatrix<float> &simMat)
{
	float id=0;
	size_t nSeqs = set.size();
	//MatrixStack<3, std::pair<float, char> > dynMat(1,1);
	std::string editString1, editString2;
	std::string::size_type len;
	size_t j,k;
	size_t pos1, pos2;
	size_t same, total;
	AlignmentMatrix<float, SimilarityMatrix<float> > dynMat(-11,-1, simMat);
	for (size_t i=0; i<nSeqs; ++i)
	{
		const typename SequenceSetType::value_type &seq1 = set[i];
		for (j=i+1; j<nSeqs; ++j)
		{
			editString1.clear();
			editString2.clear();
			const typename SequenceSetType::value_type &seq2 = set[j];
			dynMat.gotoh(seq1, seq2);
			//initGotohMatrix(seq1, seq2, simMat, dynMat);
			//runGotoh(dynMat, seq1.size(), seq2.size(), -11, -1);
			//gotohTraceback(dynMat, seq1.size(), seq2.size(), editString1, editString2);
			const auto &result = dynMat.result();
			len=result.size();
			//pos1=seq1.size();
			//pos2=seq2.size();
			same=0;
			total=0;
			const auto &editString1 = result.eS1;
			const auto &editString2 = result.eS2;
			for (k=0; k<len; ++k)
			{
				if ((editString1[k] != -1) && (editString2[k] != -1))
				{
					if (std::tolower(seq1[editString1[k]]) == std::tolower(seq2[editString2[k]]))
						++same;
					++total;
				}
				else
				{
					if (editString1[k] != -1)
						--pos1;
					if (editString2[k] != -1)
						--pos2;
				}
			}
			id+=(1.0*same/total);
		}
	}
	return 1.0*id/(nSeqs*(nSeqs-1)/2);
}




} // BioSeqDataLib

#endif /* SEQSETFUNCTIONS_HPP_ */
