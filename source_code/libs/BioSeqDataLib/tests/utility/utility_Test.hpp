/*
 * utility_Test.hpp
 *
 *  Created on: Feb 4, 2016
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

#ifndef TESTS_UTILITY_UTILITY_TEST_HPP_
#define TESTS_UTILITY_UTILITY_TEST_HPP_




#include <boost/test/unit_test.hpp>

#include "../../src/utility/utility.hpp"


BOOST_AUTO_TEST_SUITE(Utility_Test)


BOOST_AUTO_TEST_CASE( getEnv_Test)
{
	fs::path home = BioSeqDataLib::getEnv("HOME");
	BOOST_CHECK(home.string().empty() == false);
	BOOST_CHECK_THROW (BioSeqDataLib::getEnv("xyz1234uiaebiautrinetuiasevlb"), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* TESTS_UTILITY_UTILITY_TEST_HPP_ */
