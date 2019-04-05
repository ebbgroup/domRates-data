/*
 * MatrixStack_Test.hpp
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
#ifndef MATRIXSTACK_TEST_HPP_
#define MATRIXSTACK_TEST_HPP_


#include <boost/test/unit_test.hpp>

#include "../../src/utility/MatrixStack.hpp"


BOOST_AUTO_TEST_SUITE(MatrixStack_Test)


BOOST_AUTO_TEST_CASE( SimilarityMatrix_Test)
{
	BioSeqDataLib::MatrixStack<3, int> matrices;
	matrices.resize(3, 4);
	matrices[0][0][0] = 4;
	matrices[2][2][3] = 7;
	BOOST_CHECK_EQUAL(matrices[0][0][0], 4);
	BOOST_CHECK_EQUAL(matrices[2][2][3], 7);
	BOOST_CHECK_EQUAL(matrices.dim1(), 3);
	BOOST_CHECK_EQUAL(matrices.dim2(), 4);

	BioSeqDataLib::MatrixStack<3, int> matrices2(4, 4);
	BOOST_CHECK_EQUAL(matrices2.dim1(), 4);
	BOOST_CHECK_EQUAL(matrices2.dim2(), 4);

	BioSeqDataLib::MatrixStack<3, int> matrices3(7, 7, 8);
	BOOST_CHECK_EQUAL(matrices3.dim1(), 7);
	BOOST_CHECK_EQUAL(matrices3.dim2(), 7);
	BOOST_CHECK_EQUAL(matrices3[0][0][0], 8);
	BOOST_CHECK_EQUAL(matrices3[2][6][6], 8);
}

BOOST_AUTO_TEST_SUITE_END()





#endif /* MATRIXSTACK_TEST_HPP_ */
