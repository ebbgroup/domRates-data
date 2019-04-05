/*
 * DomainTest.hpp
 *
 *  Created on: 21 Nov 2013
 *      Author: Carsten Kemena
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


#ifndef BLASTHITTEST_HPP_
#define BLASTHITTEST_HPP_


#include <boost/test/unit_test.hpp>
#include <fstream>

#include "../../src/annotation/BlastHit.hpp"


BOOST_AUTO_TEST_SUITE(BlastHit_Test)

// most frequently you implement test cases as a free functions with automatic registration
BOOST_AUTO_TEST_CASE( BlastHit_Test )
{
	std::ifstream inS("../tests/annotation/data/BlastHits.txt");
	std::string line;
	for (size_t i=0; i<6; ++i)
		std::getline(inS, line);

	BioSeqDataLib::BlastHit hit(line);
	BOOST_CHECK_EQUAL(hit.queryID, "FBpp0082095");
	BOOST_CHECK_EQUAL(hit.subjectID, "scaffold979");
	BOOST_CHECK_CLOSE(hit.percentIdentity, 65.09, 0.0001);
	BOOST_CHECK_EQUAL(hit.alignmentLength, 106);
	BOOST_CHECK_EQUAL(hit.nMismatches, 37);
	BOOST_CHECK_EQUAL(hit.nGapOpenings, 0);
	BOOST_CHECK_EQUAL(hit.queryStart, 277);
	BOOST_CHECK_EQUAL(hit.queryEnd, 382);
	BOOST_CHECK_EQUAL(hit.queryStrand, '+');
	BOOST_CHECK_EQUAL(hit.subjectStart, 1240912);
	BOOST_CHECK_EQUAL(hit.subjectEnd, 1241229);
	BOOST_CHECK_EQUAL(hit.subjectStrand, '-');
	BOOST_CHECK_CLOSE(hit.eValue, 2e-41, 0.0001);
	BOOST_CHECK_EQUAL(hit.HSPBitScore, 163);
}

BOOST_AUTO_TEST_CASE(BlastRead_Test)
{
	std::map<std::string, std::vector<BioSeqDataLib::BlastHit> > hits;
	readBlast("../tests/annotation/data/BlastResult.txt", hits);
	BOOST_CHECK_EQUAL(hits.size(), 2);
	auto it = hits.begin();
	BOOST_CHECK_EQUAL(it->second.size(), 1);
	++it;
	BOOST_CHECK_EQUAL(it->second.size(), 2);
	BOOST_CHECK_EQUAL(it->first, "gnl|em|testQuery");
	BOOST_CHECK_EQUAL(it->second.begin()->queryID, "gnl|em|testQuery");
	auto it2 = it->second.begin();
	BOOST_CHECK_EQUAL(it2->subjectID, "gnl|em|test1");
	BOOST_CHECK_EQUAL(it2->eValue, 0);
}





BOOST_AUTO_TEST_SUITE_END()




#endif /* BLASTHITTEST_HPP_ */
