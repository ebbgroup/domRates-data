/*
 * Matrix_Test.hpp
 *
 *  Created on: 26 Oct 2013
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
#ifndef MATRIX_TEST_HPP_
#define MATRIX_TEST_HPP_


#include <boost/test/unit_test.hpp>

#include "../../src/utility/Matrix.hpp"


BOOST_AUTO_TEST_SUITE(Matrix_Test)


BOOST_AUTO_TEST_CASE( Matrix_Test)
{

	BioSeqDataLib::Matrix<int> mat;
	BOOST_CHECK_EQUAL(mat.dim1(), 0);
	BOOST_CHECK_EQUAL(mat.dim2(), 0);

	BioSeqDataLib::Matrix<int> mat1(4, 5);
	BOOST_CHECK_EQUAL(mat1.dim1(), 4);
	BOOST_CHECK_EQUAL(mat1.dim2(), 5);

	BioSeqDataLib::Matrix<int> mat2(4, 5, 7);
	BOOST_CHECK_EQUAL(mat2[0][0], 7);
	BOOST_CHECK_EQUAL(mat2[3][4], 7);

	mat2.resize(9,8);
	BOOST_CHECK_EQUAL(mat2.dim1(), 9);
	BOOST_CHECK_EQUAL(mat2.dim2(), 8);

	mat2.fill(9);
	BOOST_CHECK_EQUAL(mat2[0][0], 9);
	BOOST_CHECK_EQUAL(mat2[8][7], 9);


}

BOOST_AUTO_TEST_SUITE_END()

#endif /* MATRIX_TEST_HPP_ */
