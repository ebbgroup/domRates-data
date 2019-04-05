/*
 * properties.hpp
 *
 *  Created on: 26 Apr 2014
 *      Author: Carsten Kemena
 *		 Email: c.kemena[@]uni-muenster.de
 *	 Copyright: 20142013
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
 * \file properties.hpp
 * \brief .
 */

#ifndef PROPERTIES_HPP_
#define PROPERTIES_HPP_


#include <vector>


namespace BioSeqDataLib
{
/**
 * Kyte, J. and Doolittle, R. 1982. A simple method for displaying the hydropathic character of a protein. J. Mol. Biol. 157: 105-132.
 */


std::vector<float>
hydropathyScores();

}


#endif /* PROPERTIES_HPP_ */
