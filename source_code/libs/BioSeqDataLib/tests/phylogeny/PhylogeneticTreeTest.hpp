/*
 * PhylogeneticTree.hpp
 *
 *  Created on: 10 Aug 2014
 *      Author: Carsten Kemena
 *	 Copyright: 2014
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
#ifndef PhylogeneticTreeTEST_HPP_
#define PhylogeneticTreeTEST_HPP_


#include "../../src/utility/Matrix.hpp"
#include "../../src/phylogeny/PhylogeneticTree.hpp"


BOOST_AUTO_TEST_SUITE(PhylogeneticTree_Test)

// most frequently you implement test cases as a free functions with automatic registration
BOOST_AUTO_TEST_CASE( TreeNJ_Test )
{

	BioSeqDataLib::Matrix<float> distMat;
	distMat.resize(4,4);
	distMat.fill(0);
	std::vector<std::string> names = {"A", "B", "C", "D"};

	distMat[0][1] = distMat[1][0] = 3;
	distMat[0][2] = distMat[2][0] = 14;
	distMat[0][3] = distMat[3][0] = 12;
	distMat[1][2] = distMat[2][1] = 13;
	distMat[1][3] = distMat[3][1] = 11;
	distMat[2][3] = distMat[3][2] = 4;
	distMat[0][0] = 0;
	distMat[1][1] = 0;
	distMat[2][2] = 0;
	distMat[3][3] = 0;
	BioSeqDataLib::PhylogeneticTree<int> tree;
	tree.nj(distMat, names);
	BOOST_CHECK_EQUAL(tree.str() ,"((A:2.000000,B:1.000000):9.000000,C:3.000000,D:1.000000);");
}

BOOST_AUTO_TEST_CASE( stringTest_Test )
{
	BioSeqDataLib::PhylogeneticTree<int> tree;
	tree.str2tree("((A:2.000000,B:1.000000)F:9.000000,C:3.000000,D:1.000000)X;");
	BOOST_CHECK_EQUAL(tree.str() ,"((A:2.000000,B:1.000000)F:9.000000,C:3.000000,D:1.000000)X;");

}


BOOST_AUTO_TEST_CASE( TreeUPGMA_Test )
{

	BioSeqDataLib::Matrix<float> distMat;
	distMat.resize(6,6);
	distMat.fill(0);
	std::vector<std::string> names = {"A", "B", "C", "D", "E", "F"};

	distMat[0][1] = distMat[1][0] = 2;
	distMat[0][2] = distMat[2][0] = 4;
	distMat[0][3] = distMat[3][0] = 6;
	distMat[0][4] = distMat[4][0] = 6;
	distMat[0][5] = distMat[5][0] = 8;
	distMat[1][2] = distMat[2][1] = 4;
	distMat[1][3] = distMat[3][1] = 6;
	distMat[1][4] = distMat[4][1] = 6;
	distMat[1][5] = distMat[5][1] = 8;
	distMat[2][3] = distMat[3][2] = 6;
	distMat[2][4] = distMat[4][2] = 6;
	distMat[2][5] = distMat[5][2] = 8;
	distMat[3][4] = distMat[4][3] = 4;
	distMat[3][5] = distMat[5][3] = 8;
	distMat[4][5] = distMat[5][4] = 8;
	distMat[0][0] = 0;
	distMat[1][1] = 0;
	distMat[2][2] = 0;
	distMat[3][3] = 0;
	distMat[4][4] = 0;
	distMat[5][5] = 0;
	BioSeqDataLib::PhylogeneticTree<int> tree;
	tree.upgma(distMat, names);
	BOOST_CHECK_EQUAL(tree.str() ,"(F:4.000000,(((A:1.000000,B:1.000000):1.000000,C:2.000000):1.000000,(D:2.000000,E:2.000000):1.000000):1.000000);");
}

BOOST_AUTO_TEST_CASE( TreeIterator_Test )
{
	BioSeqDataLib::Matrix<float> distMat;
	distMat.resize(6,6);
	distMat.fill(0);
	std::vector<std::string> names = {"A", "B", "C", "D", "E", "F"};
	distMat[0][1] = distMat[1][0] = 2;
	distMat[0][2] = distMat[2][0] = 4;
	distMat[0][3] = distMat[3][0] = 6;
	distMat[0][4] = distMat[4][0] = 6;
	distMat[0][5] = distMat[5][0] = 8;
	distMat[1][2] = distMat[2][1] = 4;
	distMat[1][3] = distMat[3][1] = 6;
	distMat[1][4] = distMat[4][1] = 6;
	distMat[1][5] = distMat[5][1] = 8;
	distMat[2][3] = distMat[3][2] = 6;
	distMat[2][4] = distMat[4][2] = 6;
	distMat[2][5] = distMat[5][2] = 8;
	distMat[3][4] = distMat[4][3] = 4;
	distMat[3][5] = distMat[5][3] = 8;
	distMat[4][5] = distMat[5][4] = 8;
	distMat[0][0] = 0;
	distMat[1][1] = 0;
	distMat[2][2] = 0;
	distMat[3][3] = 0;
	distMat[4][4] = 0;
	distMat[5][5] = 0;
	BioSeqDataLib::PhylogeneticTree<int> tree;
	tree.upgma(distMat, names);
	auto it = tree.preorderBegin();
	auto itEnd = tree.preorderEnd();
	std::vector<std::string> nodeNames = {"", "F", "", "", "", "A", "B", "C", "", "D", "E"};
	int pos = -1;
	while (it != itEnd)
		BOOST_CHECK_EQUAL((it++)->name, nodeNames[++pos]);
	for (size_t i=0; i<nodeNames.size(); ++i)
		BOOST_CHECK_EQUAL((--it)->name, nodeNames[pos--]);
}

BOOST_AUTO_TEST_CASE( TreeRead_Test )
{
	std::string treeLine = "((Aha,bee)Arag,(cc,ff,d))whatever;";
	BioSeqDataLib::PhylogeneticTree<int> tree;
	tree.str2tree(treeLine);
	std::vector<std::string> nodeNames = {"whatever", "Arag", "Aha", "bee", "", "cc", "ff", "d"};
	int pos = -1;
	auto it = tree.preorderBegin();
	auto itEnd = tree.preorderEnd();
	while (it != itEnd)
		BOOST_CHECK_EQUAL((it++)->name, nodeNames[++pos]);

	pos = -1;
	std::vector<std::string> nodeNames2 = {"Aha", "bee", "Arag", "cc","ff", "d", "", "whatever"};
	auto postIt = tree.postorderBegin();
	auto postItEnd = tree.postorderEnd();
	while (postIt != postItEnd)
		BOOST_CHECK_EQUAL((postIt++)->name, nodeNames2[++pos]);
	for (size_t i=0; i<nodeNames.size(); ++i)
	{
		BOOST_CHECK_EQUAL((--postIt)->name, nodeNames2[pos--]);
	}


	std::string test = "((A:1.5,B:4.0)F:3,C):1.0002;";
	BioSeqDataLib::PhylogeneticTree<int> tree3;

	std::vector<std::string> nodeNames3 = {"A", "B", "F", "C",""};
	std::vector<double> edgeLen3 = {1.5, 4, 3, 1, 1.0002};
	tree3.str2tree(test);
	postIt = tree3.postorderBegin();
	postItEnd = tree3.postorderEnd();
	pos =-1;
	while (postIt != postItEnd)
	{
		BOOST_CHECK_EQUAL((postIt)->name, nodeNames3[++pos]);
		BOOST_CHECK_CLOSE((postIt++)->edgeLength, edgeLen3[pos], 0.00001);
	}

	BioSeqDataLib::PhylogeneticTree<int> tree2;
	tree2.read("../tests/phylogeny/data/tree.nwk");
	postIt = tree2.postorderBegin();
	postItEnd = tree2.postorderEnd();
	pos =-1;
	while (postIt != postItEnd)
	{
		BOOST_CHECK_EQUAL((postIt)->name, nodeNames3[++pos]);
		BOOST_CHECK_CLOSE((postIt++)->edgeLength, edgeLen3[pos], 0.00001);
	}

}



BOOST_AUTO_TEST_SUITE_END()

#endif /* PhylogeneticTreeTEST_HPP_ */
