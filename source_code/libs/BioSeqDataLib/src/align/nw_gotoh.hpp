/*
 * nw_gotoh.hpp
 *
 *  Created on: 28 Oct 2013
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
 * \file nw_gotoh.hpp
 * \brief Contains implementations of the Needleman-Wunsch and the Gotoh algorithms.
 */
#ifndef NW_GOTOH_HPP_
#define NW_GOTOH_HPP_

#include <utility>
#include <tuple>

#include "../utility/Matrix.hpp"
#include "../utility/MatrixStack.hpp"
#include "../utility/SimilarityMatrix.hpp"
#include "../utility/LineMatrix.hpp"

namespace BioSeqDataLib
{

	struct EditString
	{
		size_t start1;
		size_t end1;
		std::string s1;

		size_t start2;
		size_t end2;
		std::string s2;
	};

	/******************************************************************************
	 * 						Needleman-Wunsch Algorithm							  *
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
	initNwMatrix(const SeqType &seq1, const SeqType &seq2, const SimMat &simMat, Matrix<std::pair<DataType, char> > &matrix)
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
	initNwMatrix(const SeqType &seq1, const SeqType &seq2, const SimMat &simMat, LineMatrix<std::pair<DataType, char> > &matrix)
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
	 * \brief Fills the Needleman-Wunsch dynamic programming matrix.
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
	nw(Matrix<std::pair<DataType, char> > &matrix, size_t dim1, size_t dim2, DataType gapPenalty)
	{
		++dim2;
		auto itPrev = matrix[0].begin();
		auto itVert = matrix[0].begin();
		auto it = matrix[0].begin();
		auto itEnd = matrix[0].begin()+dim2;
		char path;
		size_t i, j=0;
		DataType score;
		for (; it!=itEnd; ++it)
			it->first=(j++)*gapPenalty;

		for (i=1; i <= dim1; ++i)
		{
			itPrev = it = matrix[i].begin();
			itEnd = matrix[i].begin()+dim2;
			itVert = matrix[i-1].begin();
			it->first = i*gapPenalty;
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
			}
		}
	}


	/**
	 * \brief Fills the Needleman-Wunsch dynamic programming matrix.
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
	nw(LineMatrix<std::pair<DataType, char> > &matrix, size_t dim1, size_t dim2, DataType gapPenalty)
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
			it->first=j*gapPenalty;

		j=0;
		for (i=1; i <= dim1; ++i)
		{
			itPrev = it = matrix.begin(i);
			itEnd+=dim2;
			it->first = i*gapPenalty;
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
			}
			++itVert;
		}
	}

	/**
	 * \brief Performs the traceback of the Needleman-Wunsch matrix.
	 * @param matrix The filled dynamic programming matrix.
	 * @param dim1 The size of the first dimension.
	 * @param dim2 The size of the second dimension.
	 * @param editString1 The gap string of sequence one.
	 * @param editString2 The gap string of sequence two.
	 */
	template<typename DataType>
	void
	nwTraceback(const Matrix<std::pair<DataType, char> > &matrix, size_t dim1, size_t dim2, std::string &editString1, std::string &editString2)
	{
		while ((dim1 != 0) && (dim2 != 0))
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

	/**
	 * \brief Performs the traceback of the Needleman-Wunsch matrix.
	 * @param matrix The filled dynamic programming matrix.
	 * @param dim1 The size of the first dimension.
	 * @param dim2 The size of the second dimension.
	 * @param editString1 The gap string of sequence one.
	 * @param editString2 The gap string of sequence two.
	 */
	template<typename DataType>
	void
	nwTraceback(const LineMatrix<std::pair<DataType, char> > &matrix, size_t dim1, size_t dim2, std::string &editString1, std::string &editString2)
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

	/******************************************************************************
	 * 								Gotoh Algorithm								  *
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
	template<typename SeqType, typename SimMatrix, typename DataType>
	void
	initGotohMatrix(const SeqType &seq1, const SeqType &seq2, SimMatrix &simMat, MatrixStack<3, std::pair<DataType, char> > &matrices)
	{
		size_t len1=seq1.size();
		size_t len2=seq2.size();
		size_t j;
		matrices.resize(len1+1, len2+1);
		auto &matrix = matrices[0];
		for (size_t i=0; i<len1; ++i)
		{
			for (j=0; j<len2; ++j)
				matrix[i+1][j+1].first = simMat.val(seq1[i],seq2[j]);
		}
	}

	/**
	 * \brief Fills the Gotoh dynamic programming matrix using a preinitialized matrix.
	 * \tparam DataType The data type to be used for the matrix values.
	 * \pre The matrix object contains the match score at each position shifted by one (e.g. position 1,1 in the matrix should contain
	 * the score for a match of positions 0,0 in the sequences).
	 * @param matrices[in|out] The matrix to fill.
	 * @param dim1[in] The size of the first dimension.
	 * @param dim2[in] The size of the second dimension.
	 * @param gop[in] The gap opening penalty to use.
	 * @param gep[in] The gap extension penalty to use.
	 * \see initGotohMatrix
	 */
	template<typename DataType, typename GapDataType>
	void
	runGotoh(MatrixStack<3, std::pair<DataType, char> > &matrices, size_t dim1, size_t dim2, GapDataType gop, GapDataType gep)
	{
		Matrix<std::pair<DataType, char> > &matM = matrices[0];
		Matrix<std::pair<DataType, char> > &matH = matrices[1];
		Matrix<std::pair<DataType, char> > &matV = matrices[2];
		typedef typename std::vector<std::pair<DataType, char> >::iterator MatrixIterator;
		MatrixIterator itM=matM[0].begin();
		MatrixIterator itH=matH[0].begin();
		MatrixIterator itV=matV[0].begin();
		MatrixIterator itVertM, itPrevH, itVertV;

		++dim1;
		++dim2;
		GapDataType use_gop;
		DataType match_score;
		const DataType MINIMUM = (-1)*(std::numeric_limits<DataType>::max()-1000);
		// 0=match, 1=insert, 2=deletion

		itM->first = 0;
		itH->first = 0;
		itV->first = 0;

		for (size_t j=1; j<dim2; ++j)
		{
			(++itM)->first = (++itH)->first=(long int)j*gep;
			(++itV)->first = MINIMUM;
		}


		for (size_t i=1; i<dim1; ++i)
		{
			itM=matM[i].begin();
			itPrevH=itH=matH[i].begin();
			itV = matV[i].begin();
			itH->first = MINIMUM;
			itM->first = itV->first=matV[i-1][0].first+gep;
			itVertM=matM[i-1].begin();
			itVertV=matV[i-1].begin();
			++itH;
			++itV;
			for (size_t j=1; j<dim2; ++j)
			{
				//calculate insert value
				match_score=itVertM->first;
				++itVertM;
				++itVertV;
				use_gop = (j==(dim2-1))? 0 : gop;
				if (itVertV->first > (itVertM->first +use_gop))
				{
					itV->second = 'v';
					itV->first = itVertV->first;
				}else
				{
					itV->second = 'm';
					itV->first = itVertM->first+use_gop;
				}
				itV->first += gep;

				//calculate deletion value
				use_gop = (i==(dim1-1))? 0 : gop;
				if (itPrevH->first > (itM->first +use_gop))
				{
					itH->second = 'h';
					itH->first = itPrevH->first;
				} else
				{
					itH->second = 'm';
					itH->first = itM->first+use_gop;
				}
				itH->first += gep;

				//calculate match value
				++itM;
				match_score+=itM->first;
				if (itV->first > itH->first)
				{
					itM->second = 'v';
					itM->first = itV->first;
				}
				else
				{
					itM->second = 'h';
					itM->first = itH->first;
				}

				if (match_score >= itM->first)
				{
					itM->second = 'm';
					itM->first = match_score;
				}
				++itH;
				++itPrevH;
				++itV;
			}
		}
	}

	/**
	 * \brief Runs a complete gotho algorithm.
	 *
	 * This function should behave like a combination of initGotoh and runGotoh.
	 * \tparam SeqType The sequence type (can be a DomainArrangement).
	 * \tparam DataType The matrix type used.
	 * \tparam SimilarityMatrix The Similarity matrix to use (should be DSM for domains).
	 * @param seq1 The first sequence/arrangement.
	 * @param seq2 The second sequence/arrangement.
	 * @param matrices The three matrices
	 * @param simMat The similarity matrix
	 * @param gop The gap opening costs
	 * @param gep The gap extension costs
	 */
	template<typename SeqType, typename DataType, typename SimilarityMatrix>
	void
	gotoh(const SeqType seq1, const SeqType &seq2, MatrixStack<3, std::pair<DataType, char> > &matrices, const SimilarityMatrix &simMat, DataType gop, DataType gep)
	{
		size_t dim1 = seq1.size() + 1;
		size_t dim2 = seq2.size() + 1;
		matrices.resize(dim1, dim2);
		Matrix<std::pair<DataType, char> > &matM = matrices[0];
		Matrix<std::pair<DataType, char> > &matH = matrices[1];
		Matrix<std::pair<DataType, char> > &matV = matrices[2];
		typedef typename std::vector<std::pair<DataType, char> >::iterator MatrixIterator;
		MatrixIterator itM=matM[0].begin();
		MatrixIterator itH=matH[0].begin();
		MatrixIterator itV=matV[0].begin();
		MatrixIterator itVertM, itPrevH, itVertV;

		DataType use_gop;
		DataType match_score;
		const DataType MINIMUM = (-1)*(std::numeric_limits<DataType>::max()-1000);
		// 0=match, 1=insert, 2=deletion

		itM->first = 0;
		itH->first = 0;
		itV->first = 0;

		for (size_t j=1; j<dim2; ++j)
		{
			(++itM)->first = (++itH)->first=(long int)j*gep;
			(++itV)->first = MINIMUM;
		}


		for (size_t i=1; i<dim1; ++i)
		{
			itM=matM[i].begin();
			itPrevH=itH=matH[i].begin();
			itV = matV[i].begin();
			itH->first = MINIMUM;
			itM->first = itV->first=matV[i-1][0].first+gep;
			itVertM=matM[i-1].begin();
			itVertV=matV[i-1].begin();
			++itH;
			++itV;
			for (size_t j=1; j<dim2; ++j)
			{
				//calculate insert value
				match_score=itVertM->first;
				++itVertM;
				++itVertV;
				use_gop = (j==(dim2-1))? 0 : gop;
				if (itVertV->first > (itVertM->first +use_gop))
				{
					itV->second = 'v';
					itV->first = itVertV->first;
				}else
				{
					itV->second = 'm';
					itV->first = itVertM->first+use_gop;
				}
				itV->first += gep;

				//calculate deletion value
				use_gop = (i==(dim1-1))? 0 : gop;
				if (itPrevH->first > (itM->first +use_gop))
				{
					itH->second = 'h';
					itH->first = itPrevH->first;
				} else
				{
					itH->second = 'm';
					itH->first = itM->first+use_gop;
				}
				itH->first += gep;

				//calculate match value
				++itM;
				match_score += simMat.val(seq1[i-1], seq2[j-1]); //itM->first;
				if (itV->first > itH->first)
				{
					itM->second = 'v';
					itM->first = itV->first;
				}
				else
				{
					itM->second = 'h';
					itM->first = itH->first;
				}

				if (match_score >= itM->first)
				{
					itM->second = 'm';
					itM->first = match_score;
				}
				++itH;
				++itPrevH;
				++itV;
			}
		}
	}

	/**
	 * \brief Performs the traceback of the Gotoh matrix.
	 * @param matrices The filled dynamic programming matrix.
	 * @param dim1 The size of the first dimension.
	 * @param dim2 The size of the second dimension.
	 * @param editString1 The gap string of sequence one.
	 * @param editString2 The gap string of sequence two.
	 */
	template<typename DataType>
	void
	gotohTraceback(MatrixStack<3, std::pair<DataType, char> > &matrices, size_t dim1, size_t dim2, std::string &editString1, std::string &editString2)
	{
		editString1.clear();
		editString2.clear();
		size_t i = dim1;
		size_t j = dim2;

		//char state='c';
		int mat = 0;
		while ((i!=0) && (j!=0))
		{
			char state = matrices[mat][i][j].second;
			if (mat==0)
			{
				if (state=='m')
				{
					--i;
					--j;
					editString1.push_back('m');
					editString2.push_back('m');
				}
				else
				{
					if (state=='v')
						mat=2;
					else
						mat=1;
				}
			}
			else
			{
				if (mat==2)
				{
					--i;
					editString1.push_back('m');
					editString2.push_back('-');
				}
				else
				{
					--j;
					editString1.push_back('-');
					editString2.push_back('m');
				}
				if (state=='m')
					mat = 0;
			}
		}

		while (j>0)
		{
			--j;
			editString1.push_back('-');
			editString2.push_back('m');
		}
		while (i>0)
		{
			--i;
			editString1.push_back('m');
			editString2.push_back('-');
		}
	}

	/**
	 * \brief Calculate identity and similarity.
	 * @param seq1 The first sequence.
	 * @param edit1 The edit sequence for the first sequence.
	 * @param seq2 The second sequence.
	 * @param edit2 The edit sequence for the second sequence.
	 * @param simMat The similarity matrix.
	 * @return Tuple containing (#identical matches, #similar matches, # of all matches)
	 */
	template<typename SeqType, typename DataType>
	std::tuple<size_t, size_t, size_t>
	idSim(const SeqType &seq1, const std::string &edit1, const SeqType &seq2, const std::string &edit2, const SimilarityMatrix<DataType> &simMat)
	{
		size_t id = 0;
		size_t sim = 0;
		size_t length = 0;
		size_t pos1 = edit1.size();
		size_t pos2 = edit2.size();
		size_t s1 = 0;
		size_t s2 = 0;
		while ((pos1>0) && (pos2>0))
		{
			--pos1;
			--pos2;
			if ((edit1[pos1] == 'm') && (edit2[pos2] == 'm'))
			{
				if (seq1[s1] == seq2[s2])
				{
					++id;
					++sim;
				}
				else
				{
					if (simMat.val(seq1[s1], seq2[s2])>0)
						++sim;
				}
				++length;
				++s1;
				++s2;
			}
			else
			{
				if (edit1[pos1] == 'm')
					++s1;
				else
					++s2;
			}
		}

		return std::tuple<size_t, size_t, size_t>(id, sim, length);
	}


	template<typename SeqType>
	std::tuple<size_t, size_t>
	id(const SeqType &seq1, const std::string &edit1, const SeqType &seq2, const std::string &edit2)
	{
		size_t id = 0;
		size_t length = 0;
		size_t pos1 = edit1.size();
		size_t pos2 = edit2.size();
		size_t s1 = 0;
		size_t s2 = 0;
		while ((pos1>0) && (pos2>0))
		{
			--pos1;
			--pos2;
			if ((edit1[pos1] == 'm') && (edit2[pos2] == 'm'))
			{
				if (seq1[s1] == seq2[s2])
					++id;
				++length;
				++s1;
				++s2;
			}
			else
			{
				if (edit1[pos1] == 'm')
					++s1;
				else
					++s2;
			}
		}

		return std::tuple<size_t, size_t>(id,length);
	}

	/**
	 * \brief Calculates the identity of two sequences.
	 * @param seq1 The first sequence.
	 * @param seq2 The second sequence.
	 * @param edit The edit string.
	 * @return  (#identical matches, # of all matches)
	 */
	template<typename SeqType>
	std::tuple<size_t, size_t>
	id(const SeqType &seq1, const SeqType &seq2, const EditString &edit)
	{
		size_t id = 0;
		size_t length = 0;
		const std::string &edit1 = edit.s1;
		const std::string &edit2 = edit.s2;
		size_t alnLength = edit1.size();

		size_t pos1 = edit.start1;
		size_t pos2 = edit.start2;
		for (size_t i=0; i<alnLength; ++i)
		{
			if ((edit1[i] == 'm') && (edit2[i] == 'm'))
			{
				if (seq1[pos1] == seq2[pos2])
					++id;
				++length;
				++pos1;
				++pos2;
			}
			else
			{
				if (edit1[i] == 'm')
					++pos1;
				else
					++pos2;
			}
		}

		return std::tuple<size_t, size_t>(id,length);
	}

	/**
	 * \brief Calculates the identity and similarity of two sequences.
	 * @param seq1 The first sequence.
	 * @param seq2 The second sequence.
	 * @param edit The edit string.
	 * @param simMat The similarity matrix.
	 * @return Tuple containing (#identical matches, #similar matches, # of all matches)
	 */
	template<typename SeqType, typename DataType>
	std::tuple<size_t, size_t, size_t>
	idSim(const SeqType &seq1, const SeqType &seq2, const EditString &edit, const SimilarityMatrix<DataType> &simMat)
	{
		size_t id = 0;
		size_t sim = 0;
		size_t length = 0;
		const std::string &edit1 = edit.s1;
		const std::string &edit2 = edit.s2;
		size_t alnLength = edit1.size();

		size_t pos1 = edit.start1;
		size_t pos2 = edit.start2;
		for (size_t i=0; i<alnLength; ++i)
		{
			if ((edit1[i] == 'm') && (edit2[i] == 'm'))
			{
				if (seq1[pos1] == seq2[pos2])
				{
					++id;
					++sim;
				}
				else
				{
					if (simMat.val(seq1[pos1], seq2[pos2]) > 0)
						++sim;
				}

				++length;
				++pos1;
				++pos2;
			}
			else
			{
				if (edit1[i] == 'm')
					++pos1;
				else
					++pos2;
			}
		}
		return std::tuple<size_t, size_t, size_t>(id, sim, length);
	}
}








#endif /* NW_GOTOH_HPP_ */
