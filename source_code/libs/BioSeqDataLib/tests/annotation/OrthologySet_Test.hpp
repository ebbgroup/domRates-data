/*
 * OrthologySet_Test.hpp
 *
 *  Created on: 18 Aug 2014
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
#ifndef ORTHOLOGYSET_TEST_HPP_
#define ORTHOLOGYSET_TEST_HPP_


#include <boost/test/unit_test.hpp>
#include <iostream>
#include "../../src/annotation/OrthologySet.hpp"


BOOST_AUTO_TEST_SUITE(OrthologySet_Test)


BOOST_AUTO_TEST_CASE( OrthologySet_TestOrthoMCL)
{
	BioSeqDataLib::OrthologySet set;
	set.read("../tests/annotation/data/orthoMCL.txt");
	auto value = set["2ants1001"];
	BOOST_CHECK_EQUAL(value.size(), 2);
	BOOST_CHECK_EQUAL(value["sinv"].size(), 50 );
	BOOST_CHECK_EQUAL(value["pbar"].size(), 1);
	BOOST_CHECK_EQUAL(value["sinv"][0],"SINV10854-PA");

	const BioSeqDataLib::OrthologySet &set2 = set;
	value = set2["2ants1001"];
	BOOST_CHECK_EQUAL(value["sinv"].size(), 50 );

	value = set.group("pbar", "PB26385-RA");
	BOOST_CHECK_EQUAL(value.size(), 2);
	BOOST_CHECK_EQUAL(value["sinv"][1], "SINV16968-PA");
	//2ants1035: pbar|PB12860-RA pbar|PB12863-RA pbar|PB26385-RA sinv|SINV15062-PA sinv|SINV16968-PA sinv|SINV16974-PA sinv|SINV16975-PA sinv|SINV16978-PA sinv|SINV16981-PA sinv|SINV23715-PA pbar|PB12165-RA pbar|PB26388-RA

	const BioSeqDataLib::OrthoGroup *group = set.find("sinv", "SINV16968-PA");
	BOOST_CHECK_EQUAL(group->size(), 2);
	BOOST_CHECK_EQUAL(group->at("pbar")[0], "PB12860-RA");
	BOOST_CHECK_EQUAL(group->at("pbar").size(), 5);
	BOOST_CHECK_EQUAL(group->at("sinv").size(), 7);

	auto specSet = set.speciesSet();
	BOOST_CHECK_EQUAL(specSet.size(), 2);
	BOOST_CHECK(specSet.find("sinv") != specSet.end());
	BOOST_CHECK(specSet.find("pbar") != specSet.end());
}

BOOST_AUTO_TEST_CASE( OrthologySet_TestProteinortho)
{
	BioSeqDataLib::OrthologySet set;
	set.read("../tests/annotation/data/proteinortho.txt");
	auto value = set["1"];
	BOOST_CHECK_EQUAL(value.size(), 2);
	BOOST_CHECK_EQUAL(value["Trichoplax_adhaerens"].size(), 1 );
	BOOST_CHECK_EQUAL(value["Trichoplax_adhaerens"][0],"XP_002107592.1");

	value = set["9"];
	BOOST_CHECK_EQUAL(value.size(), 13);
	BOOST_CHECK_EQUAL(value["Homo_sapiens"][0],"ENSP00000485772");
	BOOST_CHECK_EQUAL(value["Homo_sapiens"][1],"ENSP00000431800");

	auto specSet = set.speciesSet();
	BOOST_CHECK_EQUAL(specSet.size(), 14);
	BOOST_CHECK(specSet.find("Trichoplax_adhaerens") != specSet.end());
	BOOST_CHECK(specSet.find("Xenopus_tropicalis") != specSet.end());
	//Trichoplax_adhaerens	Schmidtea_mediterranea	Ciona_savignyi	Petromyzon_marinus	Monosiga_brevicollis	Nematostella_vectensis	Mnemiopsis_leidyi	Drosophila_melanogaster	Amphimedon_queenslandica	Hydra_vulgaris	Caenorhabditis_elegans	Saccoglossus_kowalevskii	Homo_sapiens	Xenopus_tropicalis
}


BOOST_AUTO_TEST_SUITE_END()




#endif /* ORTHOLOGYSET_TEST_HPP_ */
