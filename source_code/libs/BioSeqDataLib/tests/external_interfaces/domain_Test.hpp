/*
 * Sequence_Test.cpp
 *
 *  Created on: 15 Oct 2013
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


#include <boost/test/unit_test.hpp>
#include <iostream>


#include "../../src/external_interfaces/domainProgs.hpp"
#include "../../src/utility/Settings.hpp"
#include "../../src/domain/DomainArrangementSet.hpp"
#include "../../src/sequence/SequenceSet.hpp"

BOOST_AUTO_TEST_SUITE(DBSearch_Test)



BOOST_AUTO_TEST_CASE( Domain_Tests )
{
	// with standard evalue
	BioSeqDataLib::Settings settings;
	settings.readSettings();
	BOOST_CHECK_EQUAL(settings["pfam_db"].empty(),false);
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> domSet;
	BioSeqDataLib::runPfamScan("../tests/external_interfaces/data/hem.fasta", "", settings["pfam_db"], domSet, 1, "log.txt");
	BOOST_CHECK_EQUAL(domSet.size(), 4);
	BOOST_CHECK_EQUAL(domSet["sp|P68871|HBB_HUMAN"].size(), 1);



	//check with changed evalu
	BOOST_CHECK_EQUAL(settings["pfam_db"].empty(),false);
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> domSet2;
	BioSeqDataLib::runPfamScan("../tests/external_interfaces/data/hem.fasta", "", settings["pfam_db"], 10, domSet2, 1, "log.txt");
	BOOST_CHECK_EQUAL(domSet2.size(), 4);
	BOOST_CHECK_EQUAL(domSet2["sp|P68871|HBB_HUMAN"].size(), 2);
}



BOOST_AUTO_TEST_CASE( DomainInfo_Test )
{
	std::map<std::string, BioSeqDataLib::DomainInfo> info;
	fs::path home = BioSeqDataLib::getEnv("DOMAINWORLD_DATA");
	if (home.empty())
	{
		home = BioSeqDataLib::getEnv("HOME");
		home /= ".domainWorld";
	}
	home /= "external/pfam_info_31.txt";
	readDomainInfo(home, info);
	BOOST_CHECK_EQUAL(info.size(), 16712);
	auto it = info.find("PF00093");
	BOOST_CHECK_EQUAL(it->second.name, "VWC");
	BOOST_CHECK_EQUAL(it->second.type, "Family");
	BOOST_CHECK_EQUAL(it->second.clan, "CL0451");
}





BOOST_AUTO_TEST_SUITE_END()
