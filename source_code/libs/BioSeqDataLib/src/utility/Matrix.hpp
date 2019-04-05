/*
 * Matrix.hpp
 *
 *  Created on: 25 Oct 2013
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
 * \file Matrix.hpp
 * \brief File containing the Matrix class.
 */
#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <cstdlib>
#include <vector>


namespace BioSeqDataLib
{


/**
 * \brief The matrix class.
 * \details Simple two dimensional matrix.
 * \tparam The data to be stored.
 */
template<typename DataType>
class Matrix
{
private:
	std::vector<std::vector<DataType> > matrix_;

public:
	/**
	 * \brief Standard constructor
	 */
	Matrix();

	/**
	 * \brief Constructor initialising to a certain size.
	 * @param dim1 The size of the first dimension.
	 * @param dim2 The size of the second dimension.
	 */
	Matrix(size_t dim1, size_t dim2);

	/**
	 * \brief Constructor initialising the matrix with a certain value.
	 * @param dim1 The size of the first dimension.
	 * @param dim2 The size of the second dimension.
	 * @param init The value to initialize the matrix with.
	 */
	Matrix(size_t dim1, size_t dim2, const DataType &init);

	/**
	 * \brief Standard destructor
	 */
	virtual ~Matrix();

	/**
	 * \brief Access operator
	 * @param index The index to acess.
	 * @return Reference to the field.
	 */
	std::vector<DataType> &operator[](unsigned int index)
	{
		return matrix_[index];
	}

	/**
	 * \brief Access operator
	 * @param index The index to acess.
	 * @return Reference to the field.
	 */
	const std::vector<DataType> &operator[](unsigned int index) const
	{
		return matrix_[index];
	}

	/**
	 *  \brief Returns the size of the first dimension.
	 * @return The size of the first dimension.
	 */
	size_t
	dim1() const
	{
		return matrix_.size();
	}

	/**
	 *  \brief Returns the size of the second dimension.
	 * @return The size of the second dimension.
	 */
	size_t
	dim2() const
	{
		if (!matrix_.empty())
			return matrix_[0].size();
		return 0;
	}

	/**
	 * \brief Resizes the matrix.
	 * @param dim1 New size of the first dimension.
	 * @param dim2 New size of the second dimension.
	 */
	void
	resize(size_t dim1, size_t dim2)
	{
		size_t i;
		matrix_.resize(dim1);
		for (i=0; i<dim1; ++i)
			matrix_[i].resize(dim2);
	}

	typename std::vector<DataType>::iterator
	begin()
	{
		return matrix_.begin();
	}

	typename std::vector<DataType>::iterator
	it(size_t diff)
	{
		return (matrix_.begin()+diff);
	}


	typename std::vector<DataType>::iterator
	end()
	{
		return matrix_.end();
	}

	/**
	 * Fills the whole matrix with a single value.
	 * @param value The value to fill the matrix with.
	 */
	void
	fill(const DataType &value);
};


template<typename DataType>
Matrix<DataType>::Matrix():matrix_()
{}

template<typename DataType>
Matrix<DataType>::Matrix(size_t dim1, size_t dim2):matrix_(std::vector<std::vector<DataType> >(dim1, std::vector<DataType>(dim2)))
{}

template<typename DataType>
Matrix<DataType>::Matrix(size_t dim1, size_t dim2, const DataType &init):matrix_(std::vector<std::vector<DataType> >(dim1, std::vector<DataType>(dim2, init)))
{}

template<typename DataType>
Matrix<DataType>::~Matrix()
{}


template<typename DataType>
void
Matrix<DataType>::fill(const DataType &value)
{
	size_t j;
	size_t dim2, dim1=matrix_.size();
	if (dim1)
		dim2 = matrix_[0].size();
	for (size_t i=0; i<dim1; ++i)
	{
		for (j=0; j<dim2; ++j)
			matrix_[i][j] = value;
	}
}

} // namespace BioSeqDataLib

#endif /* MATRIX_HPP_ */
