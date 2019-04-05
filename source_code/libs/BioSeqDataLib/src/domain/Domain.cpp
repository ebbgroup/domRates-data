/*
 * Domain.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: ckeme_01
 * 	 Copyright: 2016
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

#include "Domain.hpp"


namespace BioSeqDataLib
{

Domain::Domain() : start_(0), end_(0), evalue_(0.0), domDB_(DomainDB::unknown)
{
	// TODO Auto-generated constructor stub

}


Domain::Domain(const std::string &accession, unsigned long start, unsigned long end, double evalue, DomainDB db) : accession_(accession), start_(start), end_(end), evalue_(evalue), domDB_(db)
{}

Domain::Domain(std::string &&accession, unsigned long start, unsigned long end, double evalue, DomainDB db) : accession_(std::forward<std::string>(accession)), start_(start), end_(end), evalue_(evalue), domDB_(db)
{}


Domain::~Domain()
{
	// TODO Auto-generated destructor stub
}

} /* namespace BioSeqDataLib */
