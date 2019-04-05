/*
 * dollo.hpp
 *
 *  Created on: 31.08.16
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
#ifndef DOLLOTEST_HPP_
#define DOLLOTEST_HPP_

#include <vector>

#include "../../src/phylogeny/PhylogeneticTree.hpp"
#include "../../src/phylogeny/dollo.hpp"

BOOST_AUTO_TEST_SUITE(Dollo_Test)

// most frequently you implement test cases as a free functions with automatic registration
BOOST_AUTO_TEST_CASE( Dollo_Test )
{
        BioSeqDataLib::PhylogeneticTree<std::vector<int>> nTree;
        nTree.str2tree("(((((A:1,B:1)AB:1,(C:1,D:1)CD:1)AC:1,E:1)AE:1,F:1)AF:1,G:1)R;");

        for (auto aNode=nTree.preorderBegin();aNode!=nTree.preorderEnd();++aNode) {
            if(aNode->name == "A")
            {
                aNode->data = {1,-1,-1};
            }
            else if (aNode->name == "B")
            {
                aNode->data = {-1,-1,-1};
            }
            else if (aNode->name == "C")
            {
                aNode->data = {1,-1,1};
            }
            else if (aNode->name == "D")
            {
                aNode->data = {-1,-1,-1};
            }
            else if (aNode->name == "E")
            {
                aNode->data = {-1,-1,-1};
            }
            else if (aNode->name == "F")
            {
                aNode->data = {-1,1,1};
            }
            else if (aNode->name == "G")
            {
                aNode->data = {1,-1,-1};
            }
         }

        dollo(nTree);

        BOOST_CHECK_EQUAL(nTree.str() ,"(((((A:1.000000,B:1.000000)AB:1.000000,(C:1.000000,D:1.000000)CD:1.000000)AC:1.000000,E:1.000000)AE:1.000000,F:1.000000)AF:1.000000,G:1.000000)R;");

        for (auto bNode=nTree.preorderBegin();bNode!=nTree.preorderEnd();++bNode) {
            if(bNode->name == "R"){
                BOOST_CHECK_EQUAL(bNode->data[0],1);
                BOOST_CHECK_EQUAL(bNode->data[1],-1);
                BOOST_CHECK_EQUAL(bNode->data[2],-1);
            }
            else if(bNode->name == "AF"){
                BOOST_CHECK_EQUAL(bNode->data[0],1);
                BOOST_CHECK_EQUAL(bNode->data[1],-1);
                BOOST_CHECK_EQUAL(bNode->data[2],1);
            }
            else if(bNode->name == "AE"){
                BOOST_CHECK_EQUAL(bNode->data[0],1);
                BOOST_CHECK_EQUAL(bNode->data[1],-1);
                BOOST_CHECK_EQUAL(bNode->data[2],1);
            }
            else if(bNode->name == "AC"){
                BOOST_CHECK_EQUAL(bNode->data[0],1);
                BOOST_CHECK_EQUAL(bNode->data[1],-1);
                BOOST_CHECK_EQUAL(bNode->data[2],1);
            }
            else if(bNode->name == "AB"){
                BOOST_CHECK_EQUAL(bNode->data[0],1);
                BOOST_CHECK_EQUAL(bNode->data[1],-1);
                BOOST_CHECK_EQUAL(bNode->data[2],-1);
            }
            else if(bNode->name == "CD"){
                BOOST_CHECK_EQUAL(bNode->data[0],1);
                BOOST_CHECK_EQUAL(bNode->data[1],-1);
                BOOST_CHECK_EQUAL(bNode->data[2],1);
            }
        }
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* DolloTEST_HPP_ */
