/*
 * SwTest.hpp
 *
 *  Created on: Aug 15, 2016
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
#ifndef TESTS_ALIGN_SWTEST_HPP_
#define TESTS_ALIGN_SWTEST_HPP_




#include <boost/test/unit_test.hpp>
#include <iostream>
#include <utility>


#include "../../src/sequence/Sequence.hpp"
#include "../../src/sequence/SequenceSet.hpp"
#include "../../src/utility/Matrix.hpp"
#include "../../src/utility/SimilarityMatrix.hpp"
#include "../../src/utility/utility.hpp"
#include "../../src/align/nw_gotoh.hpp"
#include "../../src/align/sw.hpp"

BOOST_AUTO_TEST_SUITE(SW_Test)

// most frequently you implement test cases as a free functions with automatic registration




BOOST_AUTO_TEST_CASE( SW_align_Test )
{
	BioSeqDataLib::Sequence<> seq1("seq1", "WWWWWWWWWWQGPYELSDETLQAPVLNDEWGTEAVFELLSNAVWWWWWWWWWWWWW", "", "test sequence");
	BioSeqDataLib::Sequence<> seq2("seq2", "PPPPPPPPPPPQGPYELSDDTNQAPVLNDEGTEAVFELLSNAVPPPPPPPPPPPP", "", "test sequence");
	BioSeqDataLib::SimilarityMatrix<float> simMat("../tests/align/data/BLOSUM62.txt");

	BioSeqDataLib::Matrix< std::pair<float,char> > matrix(1,2);
	BioSeqDataLib::initSWMatrix(seq1, seq2, simMat, matrix);
	float gapPenalty = -3;
	BioSeqDataLib::sw(matrix, seq1.size(), seq2.size(), gapPenalty);
	BioSeqDataLib::EditString editString;
	BioSeqDataLib::swTraceback(matrix, seq1.size(), seq2.size(), editString);
	BOOST_CHECK_EQUAL(editString.start1, 10);
	BOOST_CHECK_EQUAL(editString.start2, 11);
	BOOST_CHECK_EQUAL(editString.end1, 42);
	BOOST_CHECK_EQUAL(editString.end2, 42);
	BOOST_CHECK_EQUAL(editString.s1,"mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm");
	BOOST_CHECK_EQUAL(editString.s2,"mmmmmmmmmmmmmmmmmmm-mmmmmmmmmmmmm");

	auto pair = BioSeqDataLib::id(seq1, seq2, editString);
	BOOST_CHECK_EQUAL(std::get<0>(pair), 30);
	BOOST_CHECK_EQUAL(std::get<1>(pair), 32);

	auto triple = BioSeqDataLib::idSim(seq1, seq2, editString, simMat);
	BOOST_CHECK_EQUAL(std::get<0>(triple), 30);
	BOOST_CHECK_EQUAL(std::get<1>(triple), 31);
	BOOST_CHECK_EQUAL(std::get<2>(triple), 32);
}


BOOST_AUTO_TEST_SUITE_END()


#endif /* TESTS_ALIGN_SWTEST_HPP_ */
