/*
 * TreeTest.hpp
 *
 *  Created on: 24 Mar 2014
 *      Author: Carsten Kemena
 *		 Email: c.kemena[@]uni-muenster.de
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


/**
 * \file TreeTest.hpp
 * \brief .
 */

#ifndef TREETEST_HPP_
#define TREETEST_HPP_

#include "../../src/utility/Matrix.hpp"
#include "../../src/phylogeny/Tree.hpp"
#include "../../src/phylogeny/PhylogeneticTree.hpp"

BOOST_AUTO_TEST_SUITE(Tree_Test)

BOOST_AUTO_TEST_CASE( Bifurcation_Test )
{
    // bifurcating tree
    BioSeqDataLib::PhylogeneticTree<int> bifurTree;
    bifurTree.str2tree("(((C:2.000000,(A:1.000000,B:3.000000)),(D:1.000000,E:1.500000)),OG);"); // bifurcation

    bool bifurTreeBifur = BioSeqDataLib::isBifurcatingTree(bifurTree);
    BOOST_CHECK_EQUAL(bifurTreeBifur, true);

    // monofurcating tree
    BioSeqDataLib::PhylogeneticTree<int> monofurTree;
    monofurTree.str2tree("((((C),(A,B)),(D,E)),OG);"); // monofurcation

    bool monoTreeBifur = BioSeqDataLib::isBifurcatingTree(monofurTree);
    BOOST_CHECK_EQUAL(monoTreeBifur, false);

    BOOST_CHECK_THROW(BioSeqDataLib::isBifurcatingTree(monofurTree, true), std::runtime_error);

    // multifurcating tree
    BioSeqDataLib::PhylogeneticTree<int> multifurTree;
    multifurTree.str2tree("(((C,B,A),(D,E)),OG);"); // multifurcation

    bool multiTreeBifur = BioSeqDataLib::isBifurcatingTree(multifurTree);
    BOOST_CHECK_EQUAL(multiTreeBifur, false);

    BOOST_CHECK_THROW(BioSeqDataLib::isBifurcatingTree(multifurTree, true), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* TREETEST_HPP_ */
