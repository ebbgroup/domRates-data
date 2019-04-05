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


#ifndef DOMAINTEST_HPP_
#define DOMAINTEST_HPP_


#include <boost/test/unit_test.hpp>
#include <iostream>


#include "../../src/sequence/Sequence.hpp"
#include "../../src/sequence/SequenceSet.hpp"
#include "../../src/sequence/SeqFunctions.hpp"
#include "../../src/DomainModule.hpp"
#include "../../src/domain/DomainArrangementSet.hpp"


BOOST_AUTO_TEST_SUITE(Domain_Test)



BOOST_AUTO_TEST_CASE( Domain_Test )
{
	// constructor tests
	std::string acc = "Acc";
	BioSeqDataLib::Domain dom(acc, 2, 30, 0.5);
	BOOST_CHECK_EQUAL(dom.accession(), "Acc");
	BOOST_CHECK_EQUAL(dom.start(), 2);
	BOOST_CHECK_EQUAL(dom.end(), 30);
	BOOST_CHECK_CLOSE(dom.evalue(), 0.5, 0.0001);
	BOOST_CHECK_EQUAL(acc, "Acc");

	BioSeqDataLib::Domain dom2(std::move(acc), 2, 30, 0.5);
	BOOST_CHECK_EQUAL(dom2.accession(), "Acc");
	BOOST_CHECK_EQUAL(dom2.start(), 2);
	BOOST_CHECK_EQUAL(dom2.end(), 30);
	BOOST_CHECK_CLOSE(dom2.evalue(), 0.5, 0.0001);
	BOOST_CHECK_EQUAL(acc, "");

	// basic tests
	dom.accession("PF02");
	dom.start(4);
	dom.end(5);
	dom.evalue(0.6);
	BOOST_CHECK_EQUAL(dom.accession(), "PF02");
	BOOST_CHECK_EQUAL(dom.start(), 4);
	BOOST_CHECK_EQUAL(dom.end(), 5);
	BOOST_CHECK_CLOSE(dom.evalue(), 0.6, 0.0001);

	BOOST_CHECK(dom2 < dom);
	BOOST_CHECK(dom > dom2);
}





BOOST_AUTO_TEST_CASE( DomainExt_Test )
{
	BioSeqDataLib::DomainExt dom("Acc","Name", 12, 20, 10, 25, 1, 20, 3, 1.1, 1.2);
	BOOST_CHECK_EQUAL(dom.accession(), "Acc");
	BOOST_CHECK_EQUAL(dom.name(), "Name");
	BOOST_CHECK_EQUAL(dom.start(), 12);
	BOOST_CHECK_EQUAL(dom.end(), 20);
	BOOST_CHECK_EQUAL(dom.envStart(), 10);
	BOOST_CHECK_EQUAL(dom.envEnd(), 25);
	BOOST_CHECK_EQUAL(dom.hmmStart(), 1);
	BOOST_CHECK_EQUAL(dom.hmmEnd(), 20);
	BOOST_CHECK_EQUAL(dom.hmmLength(), 3);
	BOOST_CHECK_CLOSE(dom.bitscore(), 1.1, 0.0001);
	BOOST_CHECK_CLOSE(dom.evalue(), 1.2, 0.0001);

	BioSeqDataLib::DomainExt dom2("Acc","Name", 90, 100, 90, 100, 1, 20, 3, 1.1, 1.2);
	BOOST_CHECK(dom < dom2);
	BOOST_CHECK(dom2 > dom);
}




/*
BOOST_AUTO_TEST_CASE( DomainArrangement_comparison_Test )
{
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::DomainExt> a1, a2;
	a1.push_back(BioSeqDataLib::PfamDomain("A"));
	a2.push_back(BioSeqDataLib::PfamDomain("A"));
	BOOST_CHECK_EQUAL(a1==a2, true);
	BOOST_CHECK_EQUAL(a1!=a2, false);
	BOOST_CHECK_EQUAL(a1<a2, false);
	BOOST_CHECK_EQUAL(a1>a2, false);
	BOOST_CHECK_EQUAL(a1<=a2, true);
	BOOST_CHECK_EQUAL(a1>=a2, true);
	a2.push_back(BioSeqDataLib::DomainExt("B"));
	BOOST_CHECK_EQUAL(a1==a2, false);
	BOOST_CHECK_EQUAL(a1!=a2, true);
	BOOST_CHECK_EQUAL(a1<a2, true);
	BOOST_CHECK_EQUAL(a1>a2, false);
	BOOST_CHECK_EQUAL(a1<=a2, true);
	BOOST_CHECK_EQUAL(a1>=a2, false);
	a1.push_back(BioSeqDataLib::DomainExt("B"));
	BOOST_CHECK_EQUAL(a1==a2, true);
	BOOST_CHECK_EQUAL(a1!=a2, false);
	BOOST_CHECK_EQUAL(a1<a2, false);
	BOOST_CHECK_EQUAL(a1>a2, false);
	BOOST_CHECK_EQUAL(a1<=a2, true);
	BOOST_CHECK_EQUAL(a1>=a2, true);
	a1.push_back(BioSeqDataLib::DomainExt("C"));
	a2.push_back(BioSeqDataLib::DomainExt("D"));
	BOOST_CHECK_EQUAL(a1==a2, false);
	BOOST_CHECK_EQUAL(a1!=a2, true);
	BOOST_CHECK_EQUAL(a1<a2, true);
	BOOST_CHECK_EQUAL(a1>a2, false);
	BOOST_CHECK_EQUAL(a1<=a2, true);
	BOOST_CHECK_EQUAL(a1>=a2, false);

	BioSeqDataLib::DomainArrangement<BioSeqDataLib::PfamDomain> a3, a4;
	a3.push_back(BioSeqDataLib::PfamDomain("A"));
	a4.push_back(BioSeqDataLib::PfamDomain("A"));
	a4.push_back(BioSeqDataLib::PfamDomain("B"));
	a3.push_back(BioSeqDataLib::PfamDomain("B"));
	a4.push_back(BioSeqDataLib::PfamDomain("C"));
	a3.push_back(BioSeqDataLib::PfamDomain("D"));
	BOOST_CHECK_EQUAL(a3==a4, false);
	BOOST_CHECK_EQUAL(a3!=a4, true);
	BOOST_CHECK_EQUAL(a3<a4, false);
	BOOST_CHECK_EQUAL(a3>a4, true);
	BOOST_CHECK_EQUAL(a3<=a4, false);
	BOOST_CHECK_EQUAL(a3>=a4, true);
}
*/


BOOST_AUTO_TEST_CASE( DomainArrangement_other_Test )
{
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> a;
	a.push_back(BioSeqDataLib::DomainExt("A", 0, 1, 0.0));
	a.push_back(BioSeqDataLib::DomainExt("A", 0, 1, 0.0));
	a.push_back(BioSeqDataLib::DomainExt("B", 0, 1, 0.0));
	a.push_back(BioSeqDataLib::DomainExt("C", 0, 1, 0.0));
	a.push_back(BioSeqDataLib::DomainExt("C", 0, 1, 0.0));
	a.push_back(BioSeqDataLib::DomainExt("D", 0, 1, 0.0));



	BOOST_CHECK_EQUAL("A-A-B-C-C-D", a.str());
	a.collapse();
	BOOST_CHECK_EQUAL("A-B-C-D", a.str());
	BOOST_CHECK_EQUAL(4, a.size());
	a.clear();
	BOOST_CHECK_EQUAL(0, a.size());

	a.push_back(BioSeqDataLib::DomainExt("A", 0, 1, 0.0));
	a.push_back(BioSeqDataLib::DomainExt("A", 0, 1, 0.0));
	a.push_back(BioSeqDataLib::DomainExt("B", 0, 1, 0.0));
	a.push_back(BioSeqDataLib::DomainExt("C", 0, 1, 0.0));
	a.push_back(BioSeqDataLib::DomainExt("C", 0, 1, 0.0));
	a.push_back(BioSeqDataLib::DomainExt("D", 0, 1, 0.0));
	a.push_back(BioSeqDataLib::DomainExt("D", 0, 1, 0.0));

	BOOST_CHECK_EQUAL("A-A-B-C-C-D-D", a.str());
	a.collapse(true);
	BOOST_CHECK_EQUAL("A-B-C-D", a.str());
	BOOST_CHECK_EQUAL(4, a.size());
	a.reconstruct();
	BOOST_CHECK_EQUAL("A-A-B-C-C-D-D", a.str());



	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
	arrangementSet.read("../tests/domain/data/BB20012.pfamScan");
	auto set = types(arrangementSet);
	BOOST_CHECK_EQUAL(set.size(), 12);

	BOOST_CHECK_EQUAL(set.find("PF00009")!=set.end(), true);
	BOOST_CHECK_EQUAL(set.find("PF01926")!=set.end(), true);
	BOOST_CHECK_EQUAL(set.find("PF03143")!=set.end(), true);
	BOOST_CHECK_EQUAL(set.find("PF03144")!=set.end(), true);
	BOOST_CHECK_EQUAL(set.find("PF07650")!=set.end(), true);
	BOOST_CHECK_EQUAL(set.find("PF09105")!=set.end(), true);
	BOOST_CHECK_EQUAL(set.find("PF09106")!=set.end(), true);
	BOOST_CHECK_EQUAL(set.find("PF09107")!=set.end(), true);
	BOOST_CHECK_EQUAL(set.find("PF09173")!=set.end(), true);
	BOOST_CHECK_EQUAL(set.find("PF11987")!=set.end(), true);
	BOOST_CHECK_EQUAL(set.find("PF14578")!=set.end(), true);
	BOOST_CHECK_EQUAL(set.find("PF14714")!=set.end(), true);


}

/*
BOOST_AUTO_TEST_CASE( DomainArrangement_Clean_Test )
{
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> da;
	da.emplace_back("SF1", 20,50);
	da.emplace_back("PF1", 20,50);
	da.emplace_back("PS1", 20,50);
	da.emplace_back("PS2", 2,70);
	da.emplace_back("PS3", 80,550);
	da.emplace_back("PF2", 100,200);
	da.emplace_back("PS4", 500,580);
	BOOST_CHECK_EQUAL(da.size(), 5);
	std::set<std::string> keep = {"SF1"};
	std::map<std::string, int> importance = {{"PF",1}, {"SF",2}, {"SM",3}, {"PS", 4}, {"G3", 5}};
	da.clean(keep, 10, importance);
	BOOST_CHECK_EQUAL(da.size(), 2);
	BOOST_CHECK_EQUAL(da[0].accession(), "SF1");
	BOOST_CHECK_EQUAL(da[1].accession(), "PF2");
}*/

BOOST_AUTO_TEST_CASE( SplitDom_Test )
{
	//splitdomaRec tests
	std::string testprotID;
	std::string testdomaID;
	unsigned int x=10;
	unsigned int y=150;
	unsigned int k;
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
	arrangementSet.read("../tests/domain/data/splitDomData.txt");
	BOOST_CHECK_EQUAL(arrangementSet["gnl|Amel_4.5|GB44679-PA"].size(), 4);
	BioSeqDataLib::splitDomRec(arrangementSet, x, y);
	std::map <std::string, BioSeqDataLib::DomainArrangement<BioSeqDataLib::PfamDomain>>::iterator ittestsplitdoma;
	for (ittestsplitdoma = arrangementSet.begin(); ittestsplitdoma != arrangementSet.end(); ++ittestsplitdoma)
	{
		testprotID=ittestsplitdoma->first;
		auto& testsplitdoma = ittestsplitdoma->second;
		for (k=0;k<testsplitdoma.size(); ++k)
		{
			unsigned long identifier1=testsplitdoma[k].end();
			size_t identifier2=testsplitdoma[k].hmmEnd();

			if ((identifier1==236) && (identifier2==29))
			{
				BOOST_CHECK_EQUAL(testsplitdoma[k].name(), "BAH");
				BOOST_CHECK_EQUAL(testsplitdoma[k].start(), 32);
				BOOST_CHECK_EQUAL(testsplitdoma[k].end(), 236);
				BOOST_CHECK_EQUAL(testsplitdoma[k].envStart(), 32);
				BOOST_CHECK_EQUAL(testsplitdoma[k].envEnd(), 236);
				BOOST_CHECK_EQUAL(testsplitdoma[k].hmmStart(), 0);
				BOOST_CHECK_EQUAL(testsplitdoma[k].hmmEnd(), 29);
				BOOST_CHECK_EQUAL(testsplitdoma[k].hmmLength(), 121);
			}
			if ((identifier1==699) && (identifier2==29))
			{
				BOOST_CHECK_EQUAL(testsplitdoma[k].name(), "BAH");
				BOOST_CHECK_EQUAL(testsplitdoma[k].start(), 499);
				BOOST_CHECK_EQUAL(testsplitdoma[k].end(), 699);
				BOOST_CHECK_EQUAL(testsplitdoma[k].envStart(), 32);
				BOOST_CHECK_EQUAL(testsplitdoma[k].envEnd(), 236);
				BOOST_CHECK_EQUAL(testsplitdoma[k].hmmStart(), 0);
				BOOST_CHECK_EQUAL(testsplitdoma[k].hmmEnd(), 29);
				BOOST_CHECK_EQUAL(testsplitdoma[k].hmmLength(), 121);
			}
			if ((identifier1==236) && (identifier2==120))
			{
				BOOST_CHECK_EQUAL(testsplitdoma[k].name(), "BAH");
				BOOST_CHECK_EQUAL(testsplitdoma[k].start(), 32);
				BOOST_CHECK_EQUAL(testsplitdoma[k].end(), 236);
				BOOST_CHECK_EQUAL(testsplitdoma[k].envStart(), 32);
				BOOST_CHECK_EQUAL(testsplitdoma[k].envEnd(), 236);
				BOOST_CHECK_EQUAL(testsplitdoma[k].hmmStart(), 0);
				BOOST_CHECK_EQUAL(testsplitdoma[k].hmmEnd(), 120);
				BOOST_CHECK_EQUAL(testsplitdoma[k].hmmLength(), 121);
			}
			if ((identifier1==699) && (identifier2==120))
			{
				BOOST_CHECK_EQUAL(testsplitdoma[k].name(), "BAH");
				BOOST_CHECK_EQUAL(testsplitdoma[k].start(), 499);
				BOOST_CHECK_EQUAL(testsplitdoma[k].end(), 699);
				BOOST_CHECK_EQUAL(testsplitdoma[k].envStart(), 32);
				BOOST_CHECK_EQUAL(testsplitdoma[k].envEnd(), 236);
				BOOST_CHECK_EQUAL(testsplitdoma[k].hmmStart(), 0);
				BOOST_CHECK_EQUAL(testsplitdoma[k].hmmEnd(), 120);
				BOOST_CHECK_EQUAL(testsplitdoma[k].hmmLength(), 121);
			}
		}
	}
	BOOST_CHECK_EQUAL(arrangementSet["gnl|Amel_4.5|GB44679-PA"].size(), 2);
}





BOOST_AUTO_TEST_CASE( DomainArrangement_Method_Test )
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
	arrangementSet.read("../tests/domain/data/BB20012.hmmScan");
	auto arrangementSet2 = std::move(arrangementSet);
	BOOST_CHECK_EQUAL(arrangementSet2.size(), 27);
	BOOST_CHECK_EQUAL(arrangementSet.size(), 0);
}


BOOST_AUTO_TEST_CASE( Function_Test )
{
	BioSeqDataLib::Domain dom("A", 2, 30, 0.5);
	BioSeqDataLib::Domain dom2("B", 2, 30, 0.5);
	BOOST_CHECK_EQUAL(dom.distance_overlap(dom2), -29);

	BioSeqDataLib::Domain dom3("A", 10, 30, 0.5);
	BioSeqDataLib::Domain dom4("B", 35, 50, 0.5);
	BOOST_CHECK_EQUAL(dom3.distance_overlap(dom4), 4);
	BOOST_CHECK_EQUAL(dom4.distance_overlap(dom3), 4);

	BioSeqDataLib::Domain dom5("A", 10, 30, 0.5);
	BioSeqDataLib::Domain dom6("B", 25, 50, 0.5);
	BOOST_CHECK_EQUAL(dom5.distance_overlap(dom6), -6);
	BOOST_CHECK_EQUAL(dom6.distance_overlap(dom5), -6);

	BioSeqDataLib::Domain dom7("A", 10, 30, 0.5);
	BioSeqDataLib::Domain dom8("B", 31, 50, 0.5);
	BOOST_CHECK_EQUAL(dom7.distance_overlap(dom8), 0);
	BOOST_CHECK_EQUAL(dom8.distance_overlap(dom7), 0);
}




BOOST_AUTO_TEST_SUITE_END()




#endif /* DOMAINTEST_HPP_ */
