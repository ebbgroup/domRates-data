/*
 * AlignmentMatrix_Test.hpp
 *
 *  Created on: 19 Feb 2018
 *      Author: Carsten Kemena
 *		 Email: c.kemena[@]uni-muenster.de
 *	 Copyright: 2018
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


#ifndef ALIGNMENTMATRIX_TEST_HPP_
#define ALIGNMENTMATRIX_TEST_HPP_


#include <boost/test/unit_test.hpp>
#include <iostream>
#include <utility>

#include "../../src/sequence/Sequence.hpp"
#include "../../src/sequence/SequenceSet.hpp"
#include "../../src/utility/Settings.hpp"
#include "../../src/utility/SimilarityMatrix.hpp"
#include "../../src/utility/utility.hpp"
#include "../../src/align/AlignmentMatrix.hpp"
#include "../../src/DomainModule.hpp"


BOOST_AUTO_TEST_SUITE(AlignementMatrix_Test)

BOOST_AUTO_TEST_CASE( general_tests)
{
	BioSeqDataLib::AlignmentMatrix<float, BioSeqDataLib::SimilarityMatrix<float> > mat;
	mat.gop(-0.5);
	BOOST_CHECK_EQUAL(mat.gop(), -0.5);

	mat.gep(-6);
	BOOST_CHECK_EQUAL(mat.gep(), -6);
}



BOOST_AUTO_TEST_CASE( nw_sequence_test)
{
	BioSeqDataLib::Sequence<> seq1("seq1", "LMLDSGSEPKLIAEPLXPQGPYELSDETLQAPVLNDEGTEAVFELLSNAVEV", "", "test sequence");
	BioSeqDataLib::Sequence<> seq2("seq2", "LMLDSKLIAEPLPPQGPYELSDETLQAPVLNDEGTEAVFELLSNAVEVTGKEPLP", "", "test sequence");

	BioSeqDataLib::SimilarityMatrix<float> simMat("../tests/align/data/BLOSUM62.txt");
	BioSeqDataLib::AlignmentMatrix<float, BioSeqDataLib::SimilarityMatrix<float> > mat(-10, simMat);
 	mat.nw(seq1, seq2);
	std::string editString1, editString2;

	std::vector<long int> expected1 {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
		33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,-1,-1,-1,-1,-1,-1,-1};
	std::vector<long int> expected2 {0,1,2,3,4,-1,-1,-1,-1,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
		29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54};
	auto result = mat.result();
	BOOST_CHECK_EQUAL(result.size(), 59);
	BOOST_CHECK_EQUAL(result.start1, 0);
	BOOST_CHECK_EQUAL(result.start2, 0);
	BOOST_CHECK_EQUAL(result.end1, 51);
	BOOST_CHECK_EQUAL(result.end2, 54);
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());
}

BOOST_AUTO_TEST_CASE( nw_domain_test)
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

	BioSeqDataLib::AlignmentMatrix<int, BioSeqDataLib::DSM> mat;
	mat.gep(-100);
	mat.scoring(simMat);
	mat.nw(seq1, seq2);

	std::vector<long int> expected1 { 0, -1, 1};
	std::vector<long int> expected2 { 0, 1, 2};
	auto result = mat.result();
	BOOST_CHECK_EQUAL(mat.score(), 100);
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());
	
	auto alnStrings = result.alnStrings(seq1, seq2);
	std::string alnString1 = "PF02458 ******* PF00545";
	std::string alnString2 = "PF02458 PF00001 PF00545";
	BOOST_CHECK_EQUAL(alnStrings.first, alnString1);
	BOOST_CHECK_EQUAL(alnStrings.second, alnString2);

	BioSeqDataLib::AlignmentMatrix<int, BioSeqDataLib::DSM> mat2(-100, simMat);
	mat2.nw(seq1, seq2);

	result = mat.result();
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());
}

BOOST_AUTO_TEST_CASE( gotoh_align_Test )
{
	BioSeqDataLib::Sequence<> seq1("seq1", "LMLDSGSEPKLIAEPLXPQGPYELSDETLQAPVLNDEGTEAVFELLSNAVEV", "", "test sequence");
	BioSeqDataLib::Sequence<> seq2("seq2", "LLDSKLIAEPLPPQGPYELSDETLQAPVLNDEGTEAVFELLSNAVEVTGKEPLP", "", "test sequence");
	BioSeqDataLib::SimilarityMatrix<float> simMat("../tests/align/data/BLOSUM62.txt");

	BioSeqDataLib::AlignmentMatrix<float, BioSeqDataLib::SimilarityMatrix<float>> mat(-11, -1, simMat);
	mat.gotoh(seq1, seq2);

	std::vector<long int> expected1 {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
		31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,-1,-1,-1,-1,-1,-1,-1};
	std::vector<long int> expected2 {-1, 0,1,2,3,-1,-1,-1,-1,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
		31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53};

	auto result = mat.result();
	BOOST_CHECK_EQUAL(result.size(), 59);
	BOOST_CHECK_EQUAL(result.start1, 0);
	BOOST_CHECK_EQUAL(result.start2, 0);
	BOOST_CHECK_EQUAL(result.end1, 51);
	BOOST_CHECK_EQUAL(result.end2, 53);
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());
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

	BioSeqDataLib::AlignmentMatrix<int, BioSeqDataLib::DSM> mat(-50, -10, simMat);
 	mat.gotoh(seq1, seq2);

	std::vector<long int> expected1 { 0, 1, 2};
	std::vector<long int> expected2 { -1, 0, 1};
	auto result = mat.result();
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());

	BOOST_CHECK_EQUAL(mat.score(), 190.0);

	seq1.clear();
	seq2.clear();
	dom.accession("PF00244");
	seq1.push_back(dom);
	dom.accession("PF12881");
	seq1.push_back(dom);
	seq2.push_back(dom);
	mat.gotoh(seq1, seq2);

	expected1 = {0, 1};
	expected2 = {-1, 0};
	result = mat.result();
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());
	BOOST_CHECK_EQUAL(mat.score(), 90.0);


	seq1.clear();
	seq2.clear();
	dom.accession("PF00244");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF12881");
	seq1.push_back(dom);
	mat.gotoh(seq1, seq2);

	expected1 = { 0, 1};
	expected2 = { 0, -1};
	result = mat.result();
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());
	BOOST_CHECK_EQUAL(mat.score(), 90.0);

	mat.gep(-5);
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
	mat.gotoh(seq1, seq2);
	
	expected1 = { 0, -1, 1};
	expected2 = { 0, 1, 2};
	result = mat.result();
	BOOST_CHECK_EQUAL(result.size(), 3);
	BOOST_CHECK_EQUAL(result.start1, 0);
	BOOST_CHECK_EQUAL(result.start2, 0);
	BOOST_CHECK_EQUAL(result.end1, 1);
	BOOST_CHECK_EQUAL(result.end2, 2);
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());
	BOOST_CHECK_EQUAL(mat.score(), 145.0);

	mat.gotoh(seq2, seq1);
	result = mat.result();
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected2.begin(), expected2.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL(mat.score(), 145.0);
}

BOOST_AUTO_TEST_CASE(sw_align_test )
{
	BioSeqDataLib::Sequence<> seq1("seq1", "WWWWWWWWWWQGPYELSDTLQAPVLNDEWGTEAVFELLSNAVWWWWWWWWWWWWW", "", "test sequence");
	BioSeqDataLib::Sequence<> seq2("seq2", "PPPPPPPPPPPQGPYELSDDTNQAPVLNDEGTEAVFELLSNAVPPPPPPPPPPPP", "", "test sequence");
	BioSeqDataLib::SimilarityMatrix<float> simMat("../tests/align/data/BLOSUM62.txt");

	BioSeqDataLib::Matrix< std::pair<float,char> > matrix(1,2);
	BioSeqDataLib::AlignmentMatrix<float, BioSeqDataLib::SimilarityMatrix<float> > mat(-3, simMat);
	mat.sw(seq1, seq2);
	std::vector<long int> eS1, eS2;


	std::vector<long int> expected1 { 10, 11, 12, 13, 14, 15, 16, 17, -1, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41};
	std::vector<long int> expected2 { 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42};
	auto result = mat.result();
	BOOST_CHECK_EQUAL(result.size(), 33);
	BOOST_CHECK_EQUAL(result.start1, 10);
	BOOST_CHECK_EQUAL(result.start2, 11);
	BOOST_CHECK_EQUAL(result.end1, 41);
	BOOST_CHECK_EQUAL(result.end2, 42);
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());
}

/*
BOOST_AUTO_TEST_CASE(raspodom_align_domain_test )
{

	BioSeqDataLib::Settings settings;
	settings.readSettings();
	fs::path matrixName =  settings["dsm"] / "pfam-31.dsm";
	BioSeqDataLib::DSM simMat(matrixName);
	simMat.useNegative(true);

	BioSeqDataLib::DomainExt dom;
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> seq1, seq2;
	dom.accession("PF00001");
	seq1.push_back(dom);
	dom.accession("PF00002");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00003");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00001");
	seq2.push_back(dom);

	BioSeqDataLib::AlignmentMatrix<int, BioSeqDataLib::DSM> mat(-10, simMat);
 	mat.raspodom_nw(seq1, seq2);
	std::vector<long int> expected1 { 0, 1, 2};
	std::vector<long int> expected2 { 2, 0, 1};
	auto result = mat.result();
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());
	BOOST_CHECK_EQUAL(mat.isCP(), true);
	BOOST_CHECK_EQUAL(mat.score(), 500);
	BOOST_CHECK_EQUAL(result.size(), 3);

	seq1.clear();
	seq2.clear();
	dom.accession("PF00001");
	seq1.push_back(dom);
	dom.accession("PF00002");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00003");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00004");
	seq1.push_back(dom);
	dom.accession("PF00001");
	seq2.push_back(dom);

 	mat.raspodom_nw(seq1, seq2);
	expected1 = { 3, 0, 1, 2};
	expected2 = { -1, 2, 0, 1};
	result = mat.result();
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());
	BOOST_CHECK_EQUAL(mat.isCP(), true);
	BOOST_CHECK_EQUAL(mat.score(), 490);


	seq1.clear();
	seq2.clear();
	dom.accession("PF00001");
	seq1.push_back(dom);
	dom.accession("PF00002");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00003");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00004");
	seq2.push_back(dom);
	dom.accession("PF00001");
	seq2.push_back(dom);

 	mat.raspodom_nw(seq1, seq2);
	expected1 = { -1, 0, 1, 2};
	expected2 = { 2, 3, 0, 1};
	result = mat.result();
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());
	BOOST_CHECK_EQUAL(mat.isCP(), true);



	seq1.clear();
	seq2.clear();
	dom.accession("PF00001");
	seq1.push_back(dom);
	dom.accession("PF00002");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00003");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00001");
	seq2.push_back(dom);
	dom.accession("PF00002");
	seq2.push_back(dom);
	dom.accession("PF00003");
	seq2.push_back(dom);

 	mat.raspodom_nw(seq1, seq2);
	expected1 = { -1, -1, 0, 1, 2};
	expected2 = { 0, 1, 2,  3, 4};
	result = mat.result();
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS1.begin(), result.eS1.end(), expected1.begin(), expected1.end());
	BOOST_CHECK_EQUAL_COLLECTIONS(result.eS2.begin(), result.eS2.end(), expected2.begin(), expected2.end());
	BOOST_CHECK_EQUAL(mat.isCP(), false);

}*/

BOOST_AUTO_TEST_SUITE_END()

#endif /* ALIGNMENTMATRIX_TEST_HPP_ */
