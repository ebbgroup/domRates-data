/*
 * MatrixStack.hpp
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
 * \file MatrixStack.hpp
 * \brief Class to contain several matrices.
 */
#ifndef MATRIXSTACK_HPP_
#define MATRIXSTACK_HPP_

#include <vector>

#include "Matrix.hpp"

namespace BioSeqDataLib
{

/**
 * \brief Class to store multiple matrices.
 * \tparam n The number of matrices.
 * \tparam DataType The type to store.
 */
template<int n, typename DataType>
class MatrixStack
{

private:
	std::vector<Matrix<DataType> > _stack;

public:

	/**
	 * \brief Standard constructor
	 */
	MatrixStack();

	/**
	 * \brief Constructor setting the sizes.
	 * @param dim1 The size of the first dimension.
	 * @param dim2 The size of the second dimension.
	 */
	MatrixStack(size_t dim1, size_t dim2);

	/**
	 * \brief Constructor filling the matrices with a certain value.
	 * @param dim1 The size of the first dimension.
	 * @param dim2 The size of the second dimension.
	 * @param init The value to fill the matrices.
	 */
	MatrixStack(size_t dim1, size_t dim2, const DataType &init);


	/**
	 * \brief Resizes all matrices in the stack to the new dimensions.
	 * @param dim1 New size of the first dimension.
	 * @param dim2 New size of the second dimension.
	 */
	void
	resize(size_t dim1, size_t dim2)
	{
		for (size_t i=0; i<n; ++i)
			_stack[i].resize(dim1, dim2);
	}

	/**
	 * \brief Access operator
	 * @param index The index to acess.
	 * @return Reference to the field.
	 */
	Matrix<DataType> &operator[](unsigned int index)
	{
		return _stack[index];
	}

	/**
	 * \brief Access operator
	 * @param index The index to acess.
	 * @return Reference to the field.
	 */
	const Matrix<DataType> &operator[](unsigned int index) const
	{
		return _stack[index];
	}

	/**
	 * \brief Returns the size of the first dimension.
	 * @return The size of the first dimension.
	 */
	size_t dim1() const
	{
		return _stack[0].dim1();
	}

	/**
	 * \brief Returns the size of the second dimension.
	 * @return The size of the second dimension.
	 */
	size_t dim2() const
	{
		return _stack[0].dim2();
	}

};


template<int n, typename DataType>
MatrixStack<n, DataType>::MatrixStack():_stack(3, Matrix<DataType>())
{}

template<int n, typename DataType>
MatrixStack<n, DataType>::MatrixStack(size_t dim1, size_t dim2):_stack(3, Matrix<DataType>(dim1, dim2))
{}

template<int n, typename DataType>
MatrixStack<n, DataType>::MatrixStack(size_t dim1, size_t dim2, const DataType &init) :_stack(3, Matrix<DataType>(dim1, dim2, init))
{}

}


#endif /* MATRIXSTACK_HPP_ */
