/*
 * TwoValues.cpp
 *
 *  Created on: 24 Jul 2014
 *      Author: ckemena
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
#include "TwoValues.hpp"

namespace BioSeqDataLib
{

TwoValues::TwoValues():first_(0),second_(0)
{
	// TODO Auto-generated constructor stub

}

TwoValues::TwoValues(size_t first, size_t second):first_(first), second_(second)
{

}

TwoValues::~TwoValues()
{
	// TODO Auto-generated destructor stub
}

} /* namespace BioSeqDataLib */
