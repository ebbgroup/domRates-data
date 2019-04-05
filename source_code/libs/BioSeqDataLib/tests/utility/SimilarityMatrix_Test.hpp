/*
 * SimilarityMatrix.hpp
 *
 *  Created on: 27 Oct 2013
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
#ifndef SIMILARITYMATRIX_TEST_HPP_
#define SIMILARITYMATRIX_TEST_HPP_



#include <boost/test/unit_test.hpp>

#include "../../src/utility/SimilarityMatrix.hpp"


BOOST_AUTO_TEST_SUITE(SimilarityMatrix_Test)


BOOST_AUTO_TEST_CASE( SimilarityMatrix_Test)
{

	BioSeqDataLib::SimilarityMatrix<int> mat;
	mat.read("../tests/utility/data/BLOSUM62.txt");
	BOOST_CHECK_EQUAL(mat.val('A','a'), 4);
	BOOST_CHECK_EQUAL(mat.val('R','R'), 5);
	BOOST_CHECK_EQUAL(mat.val('A','r'), -1);
	BOOST_CHECK_EQUAL(mat.val('x','x'), -1);
	BOOST_CHECK_EQUAL(mat.val('*','*'), 1);
	BOOST_CHECK_EQUAL(mat.val('J','v'), -4);

}

BOOST_AUTO_TEST_CASE( SimilarityMatrixDefault_Test)
{

	BioSeqDataLib::SimilarityMatrix<int> mat;
	mat.setBlosum62();
	BOOST_CHECK_EQUAL(mat.val('A','a'), 4);
	BOOST_CHECK_EQUAL(mat.val('R','R'), 5);
	BOOST_CHECK_EQUAL(mat.val('A','r'), -1);
	BOOST_CHECK_EQUAL(mat.val('x','x'), -1);
	BOOST_CHECK_EQUAL(mat.val('*','*'), 1);
	BOOST_CHECK_EQUAL(mat.val('J','v'), -4);

}

BOOST_AUTO_TEST_SUITE_END()



#endif /* SIMILARITYMATRIX_HPP_ */
