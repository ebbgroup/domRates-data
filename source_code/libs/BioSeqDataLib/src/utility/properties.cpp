/*
 * properties.cpp
 *
 *  Created on: 3 Aug 2014
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

// C header
#include "properties.hpp"


namespace BioSeqDataLib
{
/**
 * Kyte, J. and Doolittle, R. 1982. A simple method for displaying the hydropathic character of a protein. J. Mol. Biol. 157: 105-132.
 */


std::vector<float>
hydropathyScores()
{
	std::vector<float>	scores(256, 0);
	scores['i'] = scores['I'] = 4.5;
	scores['v'] = scores['V'] = 4.2;
	scores['l'] = scores['L'] = 3.8;
	scores['f'] = scores['F'] = 2.8;
	scores['c'] = scores['C'] = 2.5;
	scores['m'] = scores['M'] = 1.9;
	scores['a'] = scores['A'] = 1.8;
	scores['g'] = scores['G'] = -0.4;
	scores['t'] = scores['T'] = -0.7;
	scores['w'] = scores['W'] = -0.9;
	scores['s'] = scores['S'] = -0.8;
	scores['y'] = scores['Y'] = -1.3;
	scores['p'] = scores['P'] = -1.6;
	scores['h'] = scores['H'] = -3.2;
	scores['e'] = scores['E'] = -3.5;
	scores['q'] = scores['Q'] = -3.5;
	scores['d'] = scores['D'] = -3.5;
	scores['n'] = scores['N'] = -3.5;
	scores['k'] = scores['K'] = -3.9;
	scores['r'] = scores['R'] = -4.5;
	return scores;
}

}


