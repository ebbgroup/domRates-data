/*
 * ArrangementTest.hpp
 *
 *  Created on: Sep 7, 2016
 *      Author: ckeme_01
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
#ifndef TESTS_DOMAIN_DOMAINARRANGEMENTTEST_HPP_
#define TESTS_DOMAIN_DOMAINARRANGEMENTTEST_HPP_


#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../../src/DomainModule.hpp"

BOOST_AUTO_TEST_SUITE(DomainArrangement_Test)



BOOST_AUTO_TEST_CASE(DA_Insert_Test)
{
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> da;

	// check overlapping of different domains
	da.emplace_back("A", 10, 100, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("B", 80, 200, 0.4, BioSeqDataLib::DomainDB::superfamily);
	da.emplace_back("C", 95, 200, 0.4, BioSeqDataLib::DomainDB::gene3d);
	BOOST_CHECK_EQUAL(da.size(), 3);

	auto it = da.begin();

	++it;
	BioSeqDataLib::Domain x("X", 10, 100, 0.4, BioSeqDataLib::DomainDB::pfam);
	it = da.insert(it, x);

	BOOST_CHECK_EQUAL(da.size(), 4);
	BOOST_CHECK_EQUAL(it-da.begin(), 1);
	BOOST_CHECK_EQUAL(it->accession(), "X");


	++it;
	it = da.emplace(it, "Y", 10, 100, 0.4, BioSeqDataLib::DomainDB::pfam);

	BOOST_CHECK_EQUAL(da.size(), 5);
	BOOST_CHECK_EQUAL(it-da.begin(), 2);
	BOOST_CHECK_EQUAL(it->accession(), "Y");

}

BOOST_AUTO_TEST_CASE( cos_test )
{
	boost::filesystem::path path = BioSeqDataLib::getEnv("DOMAINWORLD_DATA");
	BioSeqDataLib::DSM dsm(path / "dsm" / "pfam-31.dsm");
	BioSeqDataLib::DomainExt dom;
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::DomainExt> seq1, seq2;
	dom.accession("PF02458");
	seq1.push_back(dom);
	seq2.push_back(dom);
	dom.accession("PF00001");
	seq2.push_back(dom);
	dom.accession("PF00545");
	seq1.push_back(dom);
	seq2.push_back(dom);
	BOOST_CHECK_CLOSE(cos(seq1, seq2, dsm), 0.81649, 0.001);
}

BOOST_AUTO_TEST_CASE( DB_Overlap_Test )
{
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> da;

	// check overlapping of different domains
	da.emplace_back("A", 10, 100, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("B", 80, 200, 0.4, BioSeqDataLib::DomainDB::superfamily);
	da.emplace_back("C", 95, 200, 0.4, BioSeqDataLib::DomainDB::gene3d);
	da.emplace_back("D", 210, 250, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("F", 245, 255, 0.4, BioSeqDataLib::DomainDB::superfamily);
	da.emplace_back("E", 251, 290, 0.4, BioSeqDataLib::DomainDB::pfam);

	// check overlap removal based on e-value
	da.emplace_back("G", 300, 350, 0.1, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("H", 330, 400, 0.9, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("I", 400, 450, 0.9, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("J", 430, 500, 0.1, BioSeqDataLib::DomainDB::pfam);

	// check removal of domanis from db not in db list
	da.emplace_back("X", 900, 901, 0.05, BioSeqDataLib::DomainDB::prodom);
	BOOST_CHECK_EQUAL(da.size(), 11);

	std::vector<BioSeqDataLib::DomainDB> dbImportance = {BioSeqDataLib::DomainDB::pfam, BioSeqDataLib::DomainDB::superfamily, BioSeqDataLib::DomainDB::gene3d};
	da.solveDbOverlaps(dbImportance, 10, 0.10);

	BOOST_CHECK_EQUAL(da.size(), 6);
	BOOST_CHECK_EQUAL(da[0].accession(), "A");
	BOOST_CHECK_EQUAL(da[1].accession(), "C");
	BOOST_CHECK_EQUAL(da[2].accession(), "D");
	BOOST_CHECK_EQUAL(da[3].accession(), "E");
	BOOST_CHECK_EQUAL(da[4].accession(), "G");
	BOOST_CHECK_EQUAL(da[5].accession(), "J");

}

/*
 * Check that the begin and end functions of are working properly.
 * Iterator do not need to be tested themselves as they are simply hidden STL containers.
 */
BOOST_AUTO_TEST_CASE(DA_Itertator_Test)
{
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> da;
	da.emplace_back("A", 200, 400, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("B", 10, 100, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("C", 150, 170, 0.4, BioSeqDataLib::DomainDB::pfam);

	// iterator and reverse_iterator
	auto it = da.begin();
	BOOST_CHECK_EQUAL(it->accession(), "A");
	it = --da.end();
	BOOST_CHECK_EQUAL(it->accession(), "C");

	auto rit = da.rbegin();
	BOOST_CHECK_EQUAL(rit->accession(), "C");
	rit = --da.rend();
	BOOST_CHECK_EQUAL(rit->accession(), "A");

	// const_iterator and const_reverse_iterator
	auto const &const_da = da;
	auto c_it = const_da.begin();
	BOOST_CHECK_EQUAL(c_it->accession(), "A");
	c_it = --const_da.end();
	BOOST_CHECK_EQUAL(it->accession(), "C");

	auto c_rit = const_da.rbegin();
	BOOST_CHECK_EQUAL(c_rit->accession(), "C");
	c_rit = --const_da.rend();
	BOOST_CHECK_EQUAL(c_rit->accession(), "A");
}

/*
 * Check the sorting function of DomainArrangement
 */
BOOST_AUTO_TEST_CASE(Sort_check)
{
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> da;
	da.emplace_back("A", 200, 400, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("B", 10, 100, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("C", 150, 170, 0.4, BioSeqDataLib::DomainDB::pfam);
	BOOST_CHECK_EQUAL(da[0].accession(), "A");

	sort(da.begin(), da.end());
	BOOST_CHECK_EQUAL(da[0].accession(), "B");
	BOOST_CHECK_EQUAL(da[1].accession(), "C");
	BOOST_CHECK_EQUAL(da[2].accession(), "A");
}



BOOST_AUTO_TEST_CASE(Collapse_check)
{
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> da;
	da.emplace_back("A", 1, 40, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("A", 100, 101, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("C", 1000, 3000, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.collapse(false);
	BOOST_CHECK_EQUAL(da.size(), 2);
	BOOST_CHECK_EQUAL(da[0].accession(), "A");
	BOOST_CHECK_EQUAL(da[1].accession(), "C");
	
	BOOST_CHECK_EQUAL(da.str(), "A-C");
	
	
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> da2;
	da2.emplace_back("PF00007", 1, 40, 0.4, BioSeqDataLib::DomainDB::pfam);
	da2.emplace_back("PF00008", 100, 101, 0.4, BioSeqDataLib::DomainDB::pfam);
	da2.emplace_back("PF00009", 1000, 3000, 0.4, BioSeqDataLib::DomainDB::pfam);
	da2.emplace_back("PF00010", 10000, 30000, 0.4, BioSeqDataLib::DomainDB::pfam);
	da2.collapse(false);
	BOOST_CHECK_EQUAL(da2[1].accession(), "PF00008");
	BOOST_CHECK_EQUAL(da2.size(), 4);
	BOOST_CHECK_EQUAL(da2.str(), "PF00007-PF00008-PF00009-PF00010");
	
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> da3;
	da3.emplace_back("PF00007", 1, 40, 0.4, BioSeqDataLib::DomainDB::pfam);
	da3.emplace_back("PF00008", 100, 101, 0.4, BioSeqDataLib::DomainDB::pfam);
	da3.emplace_back("PF00008", 102, 104, 0.4, BioSeqDataLib::DomainDB::pfam);
	da3.emplace_back("PF00009", 1000, 3000, 0.4, BioSeqDataLib::DomainDB::pfam);
	da3.emplace_back("PF00010", 10000, 30000, 0.4, BioSeqDataLib::DomainDB::pfam);
	da3.emplace_back("PF00010", 300000, 3000000, 0.4, BioSeqDataLib::DomainDB::pfam);
	da3.collapse(false);
	BOOST_CHECK_EQUAL(da3[1].accession(), "PF00008");
	BOOST_CHECK_EQUAL(da3.size(), 4);
	BOOST_CHECK_EQUAL(da3.str(), "PF00007-PF00008-PF00009-PF00010");
	
	
	da.clear();
	da.emplace_back("A", 1, 40, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("A", 100, 101, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.emplace_back("C", 1000, 3000, 0.4, BioSeqDataLib::DomainDB::pfam);
	da.collapse(true);
	BOOST_CHECK_EQUAL(da.size(), 2);
	BOOST_CHECK_EQUAL(da[0].accession(), "A");
	BOOST_CHECK_EQUAL(da[1].accession(), "C");
	
	BOOST_CHECK_EQUAL(da.str(), "A-C");
	
	
	da2.clear();
	da2.emplace_back("PF00007", 1, 40, 0.4, BioSeqDataLib::DomainDB::pfam);
	da2.emplace_back("PF00008", 100, 101, 0.4, BioSeqDataLib::DomainDB::pfam);
	da2.emplace_back("PF00009", 1000, 3000, 0.4, BioSeqDataLib::DomainDB::pfam);
	da2.emplace_back("PF00010", 10000, 30000, 0.4, BioSeqDataLib::DomainDB::pfam);
	da2.collapse(true);
	BOOST_CHECK_EQUAL(da2[1].accession(), "PF00008");
	BOOST_CHECK_EQUAL(da2.size(), 4);
	BOOST_CHECK_EQUAL(da2.str(), "PF00007-PF00008-PF00009-PF00010");
	
	da3.clear();
	da3.emplace_back("PF00007", 1, 40, 0.4, BioSeqDataLib::DomainDB::pfam);
	da3.emplace_back("PF00008", 100, 101, 0.4, BioSeqDataLib::DomainDB::pfam);
	da3.emplace_back("PF00008", 102, 104, 0.4, BioSeqDataLib::DomainDB::pfam);
	da3.emplace_back("PF00009", 1000, 3000, 0.4, BioSeqDataLib::DomainDB::pfam);
	da3.emplace_back("PF00010", 10000, 30000, 0.4, BioSeqDataLib::DomainDB::pfam);
	da3.emplace_back("PF00010", 300000, 3000000, 0.4, BioSeqDataLib::DomainDB::pfam);
	BOOST_CHECK_EQUAL(da3.size(), 6);
	da3.collapse(true);
	BOOST_CHECK_EQUAL(da3[1].accession(), "PF00008");
	BOOST_CHECK_EQUAL(da3.size(), 4);
	BOOST_CHECK_EQUAL(da3.str(), "PF00007-PF00008-PF00009-PF00010");
}


BOOST_AUTO_TEST_SUITE_END()



#endif /* TESTS_DOMAIN_DOMAINARRANGEMENTTEST_HPP_ */
