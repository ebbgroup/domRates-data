/*
 * fitch.hpp
 *
 *  Created on: 09.08.16
 *      Author: Elias Dohmen
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
#ifndef FitchTEST_HPP_
#define FitchTEST_HPP_

#include <vector>

#include "../../src/phylogeny/PhylogeneticTree.hpp"
#include "../../src/phylogeny/fitch.hpp"

BOOST_AUTO_TEST_SUITE(Fitch_Test)

// most frequently you implement test cases as a free functions with automatic registration
BOOST_AUTO_TEST_CASE( Fitch_Test )
{
        BioSeqDataLib::PhylogeneticTree<std::vector<int>> nTree;
        nTree.str2tree("((A:2.000000,(B:1.000000,C:3.000000)X:9.000000)Y:8.000000,D:1.000000)R;");

        for (auto aNode=nTree.preorderBegin();aNode!=nTree.preorderEnd();++aNode) {
            if(aNode->name == "A")
            {
                aNode->data = {1,1,-1,1};
            }
            else if (aNode->name == "B")
            {
                aNode->data = {1,1,1,-1};
            }
            else if (aNode->name == "C")
            {
                aNode->data = {-1,1,-1,-1};
            }
            else if (aNode->name == "D")
            {
                aNode->data = {-1,1,-1,1};
            }
         }

        fitch(nTree);

        BOOST_CHECK_EQUAL(nTree.str() ,"((A:2.000000,(B:1.000000,C:3.000000)X:9.000000)Y:8.000000,D:1.000000)R;");

        for (auto bNode=nTree.preorderBegin();bNode!=nTree.preorderEnd();++bNode) {
            if(bNode->name == "R"){
                BOOST_CHECK_EQUAL(bNode->data[0],1);
                BOOST_CHECK_EQUAL(bNode->data[1],1);
                BOOST_CHECK_EQUAL(bNode->data[2],-1);
                BOOST_CHECK_EQUAL(bNode->data[3],1);
            }
            else if(bNode->name == "Y"){
                BOOST_CHECK_EQUAL(bNode->data[0],1);
                BOOST_CHECK_EQUAL(bNode->data[1],1);
                BOOST_CHECK_EQUAL(bNode->data[2],-1);
                BOOST_CHECK_EQUAL(bNode->data[3],1);
            }
            else if(bNode->name == "X"){
                BOOST_CHECK_EQUAL(bNode->data[0],1);
                BOOST_CHECK_EQUAL(bNode->data[1],1);
                BOOST_CHECK_EQUAL(bNode->data[2],-1);
                BOOST_CHECK_EQUAL(bNode->data[3],-1);
            }
            else if(bNode->name == "A"){
                BOOST_CHECK_EQUAL(bNode->data[0],1);
                BOOST_CHECK_EQUAL(bNode->data[1],1);
                BOOST_CHECK_EQUAL(bNode->data[2],-1);
                BOOST_CHECK_EQUAL(bNode->data[3],1);
            }
            else if(bNode->name == "B"){
                BOOST_CHECK_EQUAL(bNode->data[0],1);
                BOOST_CHECK_EQUAL(bNode->data[1],1);
                BOOST_CHECK_EQUAL(bNode->data[2],1);
                BOOST_CHECK_EQUAL(bNode->data[3],-1);
            }
            else if(bNode->name == "C"){
                BOOST_CHECK_EQUAL(bNode->data[0],-1);
                BOOST_CHECK_EQUAL(bNode->data[1],1);
                BOOST_CHECK_EQUAL(bNode->data[2],-1);
                BOOST_CHECK_EQUAL(bNode->data[3],-1);
            }
            else if(bNode->name == "D"){
                BOOST_CHECK_EQUAL(bNode->data[0],-1);
                BOOST_CHECK_EQUAL(bNode->data[1],1);
                BOOST_CHECK_EQUAL(bNode->data[2],-1);
                BOOST_CHECK_EQUAL(bNode->data[3],1);
            }

        }
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* FitchTEST_HPP_ */
