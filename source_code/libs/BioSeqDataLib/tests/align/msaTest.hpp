
#ifndef TESTS_ALIGN_MSATEST_HPP_
#define TESTS_ALIGN_MSATEST_HPP_




#include <boost/test/unit_test.hpp>
#include <iostream>
#include <utility>


#include "../../src/sequence/Sequence.hpp"
#include "../../src/sequence/SequenceSet.hpp"
#include "../../src/utility/Matrix.hpp"
#include "../../src/utility/SimilarityMatrix.hpp"
#include "../../src/utility/utility.hpp"
#include "../../src/align/nw_gotoh.hpp"
#include "../../src/align/sw.hpp"
#include "../../src/align/msa.hpp"

BOOST_AUTO_TEST_SUITE(MSA_Test)

BOOST_AUTO_TEST_CASE( insertGaps_Test )
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::Domain> set;
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> da;
	da.emplace_back("1", 1, 2, 0.1);
	da.emplace_back("2", 1, 2, 0.1);
	da.emplace_back("3", 1, 2, 0.1);
	da.emplace_back("4", 1, 2, 0.1);
	set.emplace("1", da);



	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> da2;
	da2.emplace_back("1", 1, 2, 0.1);
	da2.emplace_back("2", 1, 2, 0.1);
	da2.emplace_back("3", 1, 2, 0.1);
	da2.emplace_back("4", 1, 2, 0.1);
	set.emplace("2", da2);

	std::string editString = "-mm-m-m-";
	BioSeqDataLib::insertGaps(set, editString);
	auto it=set.begin()->second.begin();
	BOOST_CHECK_EQUAL(it->accession(), "");
	BOOST_CHECK_EQUAL((++it)->accession(), "1");
	BOOST_CHECK_EQUAL((++it)->accession(), "");
	BOOST_CHECK_EQUAL((++it)->accession(), "2");
	BOOST_CHECK_EQUAL((++it)->accession(), "");
	BOOST_CHECK_EQUAL((++it)->accession(), "3");
	BOOST_CHECK_EQUAL((++it)->accession(), "4");
	BOOST_CHECK_EQUAL((++it)->accession(), "");

	auto it2=++(set.begin());
	BOOST_CHECK_EQUAL(it2->second.size(), 8);
}


BOOST_AUTO_TEST_SUITE_END()

#endif
