/*
 * fitch.hpp
 *
 *  Created on: Aug 2016
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

/**
 * \file fitch.hpp
 * \brief Fitch parsimony implementation.
 */
#ifndef FITCH_HPP
#define FITCH_HPP


#include <string>
#include <vector>

#include "PhylogeneticTree.hpp"


namespace BioSeqDataLib {

    /** \addtogroup PhyloGroup
    *  @{
    */

    /**
    * \brief Infers all states of inner nodes of a phylogenetic tree according to a fitch parsimony.
    * @param phyTree The input tree. Has to be a BioSeqDataLib::PhylogeneticTree with a vector<int> in the data field.
    * \relates PhylogeneticTree
    */
    void fitch(PhylogeneticTree<std::vector < int> > &phyTree);

}

#endif //FITCH_HPP
