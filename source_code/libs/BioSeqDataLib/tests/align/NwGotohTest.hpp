/*
 * Nw_Gotoh_Test.hpp
 *
 *  Created on: 20 Nov 2013
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


#ifndef NW_GOTOH_TEST_HPP_
#define NW_GOTOH_TEST_HPP_


#include <boost/test/unit_test.hpp>
#include <iostream>
#include <utility>


#include "../../src/sequence/Sequence.hpp"
#include "../../src/sequence/SequenceSet.hpp"
#include "../../src/utility/Settings.hpp"
#include "../../src/utility/Matrix.hpp"
#include "../../src/utility/LineMatrix.hpp"
#include "../../src/utility/SimilarityMatrix.hpp"
#include "../../src/utility/DSM.hpp"
#include "../../src/utility/utility.hpp"
#include "../../src/DomainModule.hpp"
#include "../../src/align/nw_gotoh.hpp"


BOOST_AUTO_TEST_SUITE(Nw_Gotoh_Test)

// most frequently you implement test cases as a free functions with automatic registration


BOOST_AUTO_TEST_CASE( Nw_align_Test )
{
	BioSeqDataLib::Sequence<> seq1("seq1", "LMLDSGSEPKLIAEPLXPQGPYELSDETLQAPVLNDEGTEAVFELLSNAVEV", "", "test sequence");
	BioSeqDataLib::Sequence<> seq2("seq2", "LMLDSKLIAEPLPPQGPYELSDETLQAPVLNDEGTEAVFELLSNAVEVTGKEPLP", "", "test sequence");
	BioSeqDataLib::SimilarityMatrix<float> simMat("../tests/align/data/BLOSUM62.txt");

	BioSeqDataLib::Matrix< std::pair<float,char> > matrix(1,2);
	BioSeqDataLib::initNwMatrix(seq1, seq2, simMat, matrix);
	float gapPenalty = -3;
	BioSeqDataLib::nw(matrix, seq1.size(), seq2.size(), gapPenalty);
	std::string editString1, editString2;
	BioSeqDataLib::nwTraceback(matrix, seq1.size(), seq2.size(), editString1, editString2);
	BOOST_CHECK_EQUAL(editString1,"-------mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm");
	BOOST_CHECK_EQUAL(editString2,"mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm----mmmmm");
}

BOOST_AUTO_TEST_CASE( Nw_alignDomain_Test )
{

	BioSeqDataLib::DomainExt dom;
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::DomainExt> seq1, seq2;

	dom.accession("PF02458");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00001");
	seq2.push_back(dom);
	dom.accession("PF00545");
	seq1.push_back(dom);
	seq2.push_back(dom);

	boost::filesystem::path path = BioSeqDataLib::getEnv("DOMAINWORLD_DATA");
	BioSeqDataLib::DSM simMat(path / "dsm" / "pfam-31.dsm");
	simMat.useNegative(true);
	BioSeqDataLib::Matrix< std::pair<float,char> > matrix(1,2);

	BioSeqDataLib::initNwMatrix(seq1, seq2, simMat, matrix);
	float gapPenalty = -100;
	std::string editString1, editString2;
	BioSeqDataLib::nw(matrix, seq1.size(), seq2.size(), gapPenalty);
	BioSeqDataLib::nwTraceback(matrix, seq1.size(), seq2.size(), editString1, editString2);

	BOOST_CHECK_EQUAL(editString1,"m-m");
	BOOST_CHECK_EQUAL(editString2,"mmm");

}



BOOST_AUTO_TEST_CASE( Nw_LineMatrix_align_Test )
{
	BioSeqDataLib::Sequence<> seq1("seq1", "LMLDSGSEPKLIAEPLXPQGPYELSDETLQAPVLNDEGTEAVFELLSNAVEV", "", "test sequence");
	BioSeqDataLib::Sequence<> seq2("seq2", "LMLDSKLIAEPLPPQGPYELSDETLQAPVLNDEGTEAVFELLSNAVEVTGKEPLP", "", "test sequence");
	BioSeqDataLib::SimilarityMatrix<float> simMat("../tests/align/data/BLOSUM62.txt");

	BioSeqDataLib::LineMatrix< std::pair<float,char> > matrix(1,2);
	BioSeqDataLib::initNwMatrix(seq1, seq2, simMat, matrix);
	float gapPenalty = -3;
	BioSeqDataLib::nw(matrix, seq1.size(), seq2.size(), gapPenalty);
	std::string editString1, editString2;
	BioSeqDataLib::nwTraceback(matrix, seq1.size(), seq2.size(), editString1, editString2);

	BOOST_CHECK_EQUAL(editString1,"-------mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm");
	BOOST_CHECK_EQUAL(editString2,"mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm----mmmmm");
}


BOOST_AUTO_TEST_CASE( split_gotoh_align_Test )
{
	BioSeqDataLib::Sequence<> seq1("seq1", "LMLDSGSEPKLIAEPLXPQGPYELSDETLQAPVLNDEGTEAVFELLSNAVEV", "", "test sequence");
	BioSeqDataLib::Sequence<> seq2("seq2", "LLDSKLIAEPLPPQGPYELSDETLQAPVLNDEGTEAVFELLSNAVEVTGKEPLP", "", "test sequence");
	BioSeqDataLib::SimilarityMatrix<float> simMat("../tests/align/data/BLOSUM62.txt");

	BioSeqDataLib::MatrixStack<3, std::pair<float,char> > matrix(1,2);
	BioSeqDataLib::initGotohMatrix(seq1, seq2, simMat, matrix);
	float gop=-11, gep=-1;
	BioSeqDataLib::runGotoh(matrix, seq1.size(), seq2.size(), gop, gep);
	std::string editString1, editString2;
	BioSeqDataLib::gotohTraceback(matrix, seq1.size(), seq2.size(), editString1, editString2);


	BOOST_CHECK_EQUAL(editString1,"-------mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm");
	BOOST_CHECK_EQUAL(editString2,"mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm----mmmm-");
}

BOOST_AUTO_TEST_CASE( gotoh_align_Test )
{
	BioSeqDataLib::Sequence<> seq1("seq1", "LMLDSGSEPKLIAEPLXPQGPYELSDETLQAPVLNDEGTEAVFELLSNAVEV", "", "test sequence");
	BioSeqDataLib::Sequence<> seq2("seq2", "LLDSKLIAEPLPPQGPYELSDETLQAPVLNDEGTEAVFELLSNAVEVTGKEPLP", "", "test sequence");
	BioSeqDataLib::SimilarityMatrix<float> simMat("../tests/align/data/BLOSUM62.txt");

	BioSeqDataLib::MatrixStack<3, std::pair<float,char> > matrix(1,2);
	float gop=-11, gep=-1;
	BioSeqDataLib::gotoh(seq1, seq2, matrix, simMat, gop, gep);
	std::string editString1, editString2;
	BioSeqDataLib::gotohTraceback(matrix, seq1.size(), seq2.size(), editString1, editString2);


	BOOST_CHECK_EQUAL(editString1,"-------mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm");
	BOOST_CHECK_EQUAL(editString2,"mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm----mmmm-");
}


BOOST_AUTO_TEST_CASE( Gotoh_alignDomain_Test )
{
	BioSeqDataLib::DomainExt dom;
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::DomainExt> seq1, seq2;

	dom.accession("PF02458");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00001");
	seq2.push_back(dom);
	dom.accession("PF00545");
	seq1.push_back(dom);
	seq2.push_back(dom);

	BioSeqDataLib::Settings settings;
	settings.readSettings();
	fs::path matrixName =  settings["dsm"] / "pfam-31.dsm";
	BioSeqDataLib::DSM simMat(matrixName);
	simMat.useNegative(true);
	BioSeqDataLib::MatrixStack<3, std::pair<float,char> > matrix(1,2);
	BioSeqDataLib::initGotohMatrix(seq1, seq2, simMat, matrix);
	BioSeqDataLib::runGotoh(matrix, seq1.size(), seq2.size(), -100, -10);
	std::string editString1, editString2;
	BioSeqDataLib::gotohTraceback(matrix, seq1.size(), seq2.size(), editString1, editString2);

	BOOST_CHECK_EQUAL(editString1,"m-m");
	BOOST_CHECK_EQUAL(editString2,"mmm");
}


BOOST_AUTO_TEST_CASE( Gotoh_alignDomain_Test2 )
{
	BioSeqDataLib::Settings settings;
	settings.readSettings();
	fs::path matrixName =  settings["dsm"] / "pfam-31.dsm";
	BioSeqDataLib::DSM simMat(matrixName);
	simMat.useNegative(true);

	BioSeqDataLib::DomainExt dom;
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> seq1, seq2;


	dom.accession("PF00244");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00244");
	seq1.push_back(dom);
	seq1.push_back(dom);
	seq2.push_back(dom);
	std::string editString1, editString2;
	BioSeqDataLib::MatrixStack<3, std::pair<int,char> > matrix(1,2);
	BioSeqDataLib::gotoh(seq1, seq2, matrix, simMat, -50, -10);
	BioSeqDataLib::gotohTraceback(matrix, seq1.size(), seq2.size(), editString1, editString2);

	BOOST_CHECK_EQUAL(editString1,"mmm");
	BOOST_CHECK_EQUAL(editString2,"mm-");
	BOOST_CHECK_EQUAL(matrix[0][3][2].first, 190.0);
	editString1.clear();
	editString2.clear();

	seq1.clear();
	seq2.clear();
	dom.accession("PF00244");
	seq1.push_back(dom);
	dom.accession("PF12881");
	seq1.push_back(dom);
	seq2.push_back(dom);
	BioSeqDataLib::gotoh(seq1, seq2, matrix, simMat, -50, -10);
	BioSeqDataLib::gotohTraceback(matrix, seq1.size(), seq2.size(), editString1, editString2);

	BOOST_CHECK_EQUAL(editString1,"mm");
	BOOST_CHECK_EQUAL(editString2,"m-");
	BOOST_CHECK_EQUAL(matrix[0][seq1.size()][seq2.size()].first, 90.0);


	seq1.clear();
	seq2.clear();
	dom.accession("PF00244");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF12881");
	seq1.push_back(dom);
	BioSeqDataLib::gotoh(seq1, seq2, matrix, simMat, -50, -10);
	BioSeqDataLib::gotohTraceback(matrix, seq1.size(), seq2.size(), editString1, editString2);

	BOOST_CHECK_EQUAL(editString1,"mm");
	BOOST_CHECK_EQUAL(editString2,"-m");
	BOOST_CHECK_EQUAL(matrix[0][seq1.size()][seq2.size()].first, 90.0);


	seq1.clear();
	seq2.clear();
	dom.accession("PF02458");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00001");
	seq2.push_back(dom);
	dom.accession("PF00545");
	seq1.push_back(dom);
	seq2.push_back(dom);
	BioSeqDataLib::gotoh(seq1, seq2, matrix, simMat, -50, -5);
	BioSeqDataLib::gotohTraceback(matrix, seq1.size(), seq2.size(), editString1, editString2);
	BOOST_CHECK_EQUAL(editString1,"m-m");
	BOOST_CHECK_EQUAL(editString2,"mmm");
	BOOST_CHECK_EQUAL(matrix[0][seq1.size()][seq2.size()].first, 145.0);


	BioSeqDataLib::gotoh(seq2, seq1, matrix, simMat, -50, -5);
	BioSeqDataLib::gotohTraceback(matrix, seq2.size(), seq1.size(), editString1, editString2);
	BOOST_CHECK_EQUAL(editString1,"mmm");
	BOOST_CHECK_EQUAL(editString2,"m-m");
	BOOST_CHECK_EQUAL(matrix[0][seq2.size()][seq1.size()].first, 145.0);

}

BOOST_AUTO_TEST_CASE( idTest )
{
	BioSeqDataLib::Sequence<> seq1("seq1", "AIIEERI", "", "test sequence");
	BioSeqDataLib::Sequence<> seq2("seq2", "AIIQ--E", "", "test sequence");
	std::string edit1 = "mmmmmmm";
	std::string edit2 = "m--mmmm";
	BioSeqDataLib::SimilarityMatrix<float> simMat;
	simMat.setBlosum62();

	auto val = idSim(seq1, edit1, seq2, edit2, simMat);
	BOOST_CHECK_EQUAL(std::get<0>(val), 3);
	BOOST_CHECK_EQUAL(std::get<1>(val), 4);
	BOOST_CHECK_EQUAL(std::get<2>(val), 5);

	auto val2 = id(seq1, edit1, seq2, edit2);
	BOOST_CHECK_EQUAL(std::get<0>(val2), 3);
	BOOST_CHECK_EQUAL(std::get<1>(val2), 5);

}

BOOST_AUTO_TEST_SUITE_END()


#endif /* NW_GOTOH_TEST_HPP_ */
