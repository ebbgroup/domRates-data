/*
 * sw.hpp
 *
 *  Created on: Aug 12, 2016
 *      Author: ckeme_01
 *	 Copyright: 2016
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
#ifndef SRC_ALIGN_SW_HPP_
#define SRC_ALIGN_SW_HPP_



#include <utility>
#include <tuple>

#include "../utility/Matrix.hpp"
#include "../utility/MatrixStack.hpp"
#include "../utility/SimilarityMatrix.hpp"
#include "../utility/LineMatrix.hpp"

#include "nw_gotoh.hpp"

namespace BioSeqDataLib
{

	/******************************************************************************
	 * 						Smith-Waterman Algorithm							  *
	 ******************************************************************************/

	/**
	 * \Fills the Matrix with matches.
	 * \tparam SeqType The type of sequences to align.
	 * \tparam DataType The type of the matrix.
	 * @param seq1 The first sequence.
	 * @param seq2 The second sequence.
	 * @param simMat The similarity matrix.
	 * @param matrix The dynamic programming matrix.
	 */
	template<typename SeqType, typename SimMat, typename DataType >
	void
	initSWMatrix(const SeqType &seq1, const SeqType &seq2, const SimMat &simMat, Matrix<std::pair<DataType, char> > &matrix)
	{
		size_t len1=seq1.size();
		size_t len2=seq2.size();
		size_t j;
		matrix.resize(len1+1, len2+1);
		for (size_t i=0; i<len1; ++i)
		{
			for (j=0; j<len2; ++j)
				matrix[i+1][j+1].first = simMat.val(seq1[i],seq2[j]);
		}
	}

	/**
	 * \Fills the Matrix with matches.
	 * \tparam SeqType The type of sequences to align.
	 * \tparam DataType The type of the matrix.
	 * @param seq1 The first sequence.
	 * @param seq2 The second sequence.
	 * @param simMat The similarity matrix.
	 * @param matrix The dynamic programming matrix.
	 */
	template<typename SeqType, typename SimMat, typename DataType>
	void
	initSWMatrix(const SeqType &seq1, const SeqType &seq2, const SimMat &simMat, LineMatrix<std::pair<DataType, char> > &matrix)
	{
		size_t len1=seq1.size();
		size_t len2=seq2.size();
		size_t j;
		matrix.resize(len1+1, len2+1);
		for (size_t i=0; i<len1; ++i)
		{
			for (j=0; j<len2; ++j)
				matrix.val(i+1,j+1).first = simMat.val(seq1[i],seq2[j]);
		}
	}


	/**
	 * \brief Fills the Smith-Waterman dynamic programming matrix.
	 * \tparam DataType The data type to be used for the matrix values.
	 * \pre The matrix object contains the match score at each position shifted by one (e.g. position 1,1 in the matrix should contain
	 * the score for a match of positions 0,0 in the sequences).
	 * @param matrix[fin|out] The matrix to fill.
	 * @param dim1[in] The size of the first dimension.
	 * @param dim2[in] The size of the second dimension.
	 * @param gapPenalty[in] The gap penalty to use.
	 * \see initNwMatrix
	 */
	template<typename DataType>
	void
	sw(Matrix<std::pair<DataType, char> > &matrix, size_t dim1, size_t dim2, DataType gapPenalty)
	{
		++dim2;
		auto itPrev = matrix[0].begin();
		auto itVert = matrix[0].begin();
		auto it = matrix[0].begin();
		auto itEnd = matrix[0].begin()+dim2;
		char path;
		size_t i;//, j=0;
		DataType score;
		for (; it!=itEnd; ++it)
		{
			it->first=0;
			it->second='o';
		}

		for (i=1; i <= dim1; ++i)
		{
			itPrev = it = matrix[i].begin();
			itEnd = matrix[i].begin()+dim2;
			itVert = matrix[i-1].begin();
			it->first = 0;
			it->second = 'o';
			++it;
			for (; it!=itEnd; ++it, ++itPrev)
			{
				score=it->first+itVert->first;
				path=itVert->second;
				++itVert;
				// check gaps;
				if (itVert->first > itPrev->first)
				{
					it->first = itVert->first;
					it->second = 'i';
				}
				else
				{
					it->first = itPrev->first;
					it->second = 'j';
				}
				it->first += gapPenalty;
				// check match
				if ((score > it->first) || ((score == it->first) && (path=='m')))
				{
					it->first = score;
					it->second = 'm';
				}

				if (it->first <= 0)
				{
					it->first = 0;
					it->second = 'o';
				}
			}
		}
	}


	/**
	 * \brief Fills the Smith-Waterman dynamic programming matrix.
	 * \tparam DataType The data type to be used for the matrix values.
	 * \pre The matrix object contains the match score at each position shifted by one (e.g. position 1,1 in the matrix should contain
	 * the score for a match of positions 0,0 in the sequences).
	 * @param matrix[fin|out] The matrix to fill.
	 * @param dim1[in] The size of the first dimension.
	 * @param dim2[in] The size of the second dimension.
	 * @param gapPenalty[in] The gap penalty to use.
	 * \see initNwMatrix
	 */
/*	template<typename DataType>
	void
	sw(LineMatrix<std::pair<DataType, char> > &matrix, size_t dim1, size_t dim2, DataType gapPenalty)
	{
		size_t i;
		DataType score;
		++dim2;
		auto itPrev = matrix.begin();
		auto itVert = matrix.begin();
		auto it = matrix.begin();
		auto itEnd = matrix.it(dim2);
		char path;
		size_t j=0;
		for (; it!=itEnd; ++it)
		{
			it->first=0;
			it->second = 'o';
		}
		j=0;
		for (i=1; i <= dim1; ++i)
		{
			itPrev = it = matrix.begin(i);
			itEnd+=dim2;
			it->first = 0;
			it->second = 'o';
			++it;
			for (; it!=itEnd; ++it, ++itPrev)
			{
				score=it->first+itVert->first;
				path=itVert->second;
				++itVert;
				// check gaps;
				if (itVert->first > itPrev->first)
				{
					it->first = itVert->first;
					it->second = 'i';
				}
				else
				{
					it->first = itPrev->first;
					it->second = 'j';
				}
				it->first += gapPenalty;
				// check match
				if ((score > it->first) || ((score == it->first) && (path=='m')))
				{
					it->first = score;
					it->second = 'm';
				}

				if (it->first < 0)
				{
					it->first = 0;
					it->second = 'o';
				}
			}
			++itVert;
		}
	}*/

	/**
	 * \brief Performs the traceback of the Smith-Waterman matrix.
	 * @param matrix The filled dynamic programming matrix.
	 * @param dim1 The size of the first dimension.
	 * @param dim2 The size of the second dimension.
	 * @param editString1 The gap string of sequence one.
	 * @param editString2 The gap string of sequence two.
	 */
	template<typename DataType>
	void
	swTraceback(const Matrix<std::pair<DataType, char> > &matrix, size_t dim1, size_t dim2, EditString &editString)
	{
		// find max
		DataType maximumVal = 0;
		size_t d1=0, d2=0;
		for (size_t i=0; i<dim1; ++i)
		{
			for (size_t j=0; j<dim2; ++j)
			{
				if (matrix[i][j].first > maximumVal)
				{
					maximumVal = matrix[i][j].first;
					d1 = i;
					d2 = j;
				}
			}
		}

		dim1 = d1;
		dim2 = d2;
		editString.s1.clear();
		editString.s2.clear();
		editString.end1 = d1-1;
		editString.end2 = d2-1;
		std::string &editString1 = editString.s1;
		std::string &editString2 = editString.s2;


		while (matrix[dim1][dim2].second != 'o')
		{
			if (matrix[dim1][dim2].second == 'm')
			{
				editString1.push_back('m');
				editString2.push_back('m');
				--dim1;
				--dim2;
			}
			else
			{
				if (matrix[dim1][dim2].second == 'i')
				{
					editString1.push_back('m');
					editString2.push_back('-');
					--dim1;
				}
				else
				{
					editString1.push_back('-');
					editString2.push_back('m');
					--dim2;
				}
			}
		}

		std::reverse(editString1.begin(), editString1.end());
		std::reverse(editString2.begin(), editString2.end());
		editString.start1 = dim1;
		editString.start2 = dim2;
	}


	/**
	 * \brief Performs the traceback of the Smith-Waterman matrix.
	 * @param matrix The filled dynamic programming matrix.
	 * @param dim1 The size of the first dimension.
	 * @param dim2 The size of the second dimension.
	 * @param editString1 The gap string of sequence one.
	 * @param editString2 The gap string of sequence two.
	 */
	/*template<typename DataType>
	void
	swTraceback(const LineMatrix<std::pair<DataType, char> > &matrix, size_t dim1, size_t dim2, std::string &editString1, std::string &editString2)
	{
		while ((dim1 != 0) && (dim2 != 0))
		{
			if (matrix.val(dim1,dim2).second == 'm')
			{
				editString1.push_back('m');
				editString2.push_back('m');
				--dim1;
				--dim2;
			}
			else
			{
				if (matrix.val(dim1,dim2).second == 'i')
				{
					editString1.push_back('m');
					editString2.push_back('-');
					--dim1;
				}
				else
				{
					editString1.push_back('-');
					editString2.push_back('m');
					--dim2;
				}
			}
		}

		// reached first position of a sequence, fill remaining length with gaps
		if (dim1 != 0)
		{
			editString1.append(dim1, 'm');
			editString2.append(dim1, '-');
		}
		if (dim2 != 0)
		{
			editString1.append(dim2, '-');
			editString2.append(dim2, 'm');
		}
	}

*/
}




#endif /* SRC_ALIGN_SW_HPP_ */
