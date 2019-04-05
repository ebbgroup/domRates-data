/*
 * Alphabet.hpp
 *
 *  Created on: 16 Oct 2013
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

#ifndef ALPHABET_HPP_
#define ALPHABET_HPP_

namespace BioSeqDataLib
{

enum class DNA {A, C, G, T, R, Y, S, W, K, M, B, D, H, V, N, GAP, GAP2};

enum class RNA {A, C, G, U, R, Y, S, W, K, M, B, D, H, V, N, GAP, GAP2};

enum class AminoAcid {A, B, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, X, Y, Z, GAP, GAP2};

} // BioSeqDataLib

#endif /* ALPHABET_HPP_ */
