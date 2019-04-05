/*
 * LineMatrix.hpp
 *
 *  Created on: 23 Nov 2013
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
 * \file LineMatrix.hpp
 * \brief File containing the LineMatrix.
 */
#ifndef LINEMATRIX_HPP_
#define LINEMATRIX_HPP_

#include<vector>

namespace BioSeqDataLib
{

/**
 * \brief A two-dimensional matrix in a single line.
 */
template<typename DataType>
class LineMatrix
{

private:
	size_t dim1_;
	size_t dim2_;
	std::vector<DataType> matrix_;

public:

	/**
	 * \brief Constructor for a matrix.
	 * @param i The size of the first matrix.
	 * @param j The size of the second matrix.
	 */
	LineMatrix(size_t i, size_t j) : dim1_(i), dim2_(j), matrix_(i*j)
	{}

	/**
	 * \brief Constructor for a matrix with a certain initialization.
	 * @param i The size of the first matrix.
	 * @param j The size of the second matrix.
	 * @param value Value to initialize the matrix with.
	 */
	LineMatrix(size_t i, size_t j, const DataType &value) : dim1_(i), dim2_(j), matrix_(i*j, value)
	{}

	/**
	 * \brief Access to a value in the matrix.
	 * @param i The index in the first dimension.
	 * @param j The index in the second dimension.
	 * @return Reference to the value.
	 */
	DataType&
	val(size_t i, size_t j)
	{
		return matrix_[i*dim2_+j];
	}

	/**
	 * \brief Access to a value in the matrix.
	 * @param i The index in the first dimension.
	 * @param j The index in the second dimension.
	 * @return Const reference to the value.
	 */
	const DataType&
	val(size_t i, size_t j) const
	{
		return matrix_[i*dim2_+j];
	}

	/**
	 * Sets a new value.
	 * @param i The index in the first dimension.
	 * @param j The index in the second dimension.
	 * @param val The new value.
	 */
	void
	val(size_t i, size_t j, DataType &&val)
	{
		matrix_[i*dim2_+j]=std::forward<DataType>(val);
	}


	/**
	 * Resizes the matrix two a new dimension.
	 * @param i New size of the first dimension.
	 * @param j New size of the second dimension.
	 * \warning The indices change, so the old values will be lost.
	 */
	void
	resize(size_t i, size_t j)
	{
		dim1_ = i;
		dim2_ = j;
		matrix_.resize(i*j);
	}

	/**
	 * Returns the size of the first dimension.
	 * @return The size of the first dimension.
	 */
	size_t dim1() const
	{
		return dim1;
	}

	/**
	 * Returns the size of the second dimension.
	 * @return The size of the second dimension.
	 */
	size_t dim2() const
	{
		return dim2;
	}

	/**
	 * \brief Returns an iterator to the beginning of the matrix.
	 * @return Iterator to the beginning of the matrix.
	 */
	typename std::vector<DataType>::iterator
	begin()
	{
		return matrix_.begin();
	}

	typename std::vector<DataType>::iterator
	begin(size_t position)
	{
		return (matrix_.begin()+position*dim2_);
	}


	/**
	 * \brief Iterator to a specific position.
	 * @param position The position in the matrix.
	 * @return Iterator to a specific position.
	 */
	typename std::vector<DataType>::iterator
	it(size_t position)
	{
		return (matrix_.begin()+position);
	}



	/**
	 * \brief Returns an iterator to the end of the matrix.
	 * @return Iterator to the end of the matrix.
	 */
	typename std::vector<DataType>::iterator
	end()
	{
		return matrix_.end();
	}

	typename std::vector<DataType>::iterator
	end(size_t position)
	{
		return (matrix_.begin()+position*(dim2_+1));
	}
};

}

#endif /* LINEMATRIX_HPP_ */
