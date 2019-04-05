/*
 * DSM_Test.hpp
 *
 *  Created on: Nov 27, 2015
 *      Author: ckeme_01
 */

#ifndef TESTS_UTILITY_DSM_TEST_HPP_
#define TESTS_UTILITY_DSM_TEST_HPP_


#include "../../src/utility/DSM.hpp"
#include "../../src/domain/PfamDomain.hpp"

#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(LogNumber_Test)


BOOST_AUTO_TEST_CASE( DSM_Test)
{
	BioSeqDataLib::DSM mat;
	mat.read("../tests/utility/data/DSMtest.mat");
	BOOST_CHECK_EQUAL(mat.val("PF12279","PF02965"), 11);
	BOOST_CHECK_EQUAL(mat.val("PF02458","PF02458"), 100);
	BOOST_CHECK_EQUAL(mat.val("PF02458","PF12279"), 51);
	BOOST_CHECK_EQUAL(mat.val("PF02458","PF02965"), 0);

	BOOST_CHECK_EQUAL(mat.threshold(), 10);

	mat.useNegative(true);
	BOOST_CHECK_EQUAL(mat.val("PF02458","PF02965"), -100);
	BOOST_CHECK_EQUAL(mat.val("PF02458","PF02458"), 100);
	BOOST_CHECK_EQUAL(mat.val("PF12279","PF02965"), -78);
	BOOST_CHECK_EQUAL(mat.val("PF02458","PF12279"), 2);

	BioSeqDataLib::DomainExt dom1, dom2;
	dom1.accession("PF02458");
	dom2.accession("PF02965");
	BOOST_CHECK_EQUAL(mat.val(dom1, dom2), -100);
	BioSeqDataLib::PfamDomain dom3, dom4;
	dom3.accession("PF02458");
	dom4.accession("PF02458");
	BOOST_CHECK_EQUAL(mat.val(dom3, dom4), 100);

	BOOST_CHECK_THROW(mat.val("xxx", "xxx"), std::out_of_range);
}


BOOST_AUTO_TEST_SUITE_END()





#endif /* TESTS_UTILITY_DSM_TEST_HPP_ */
