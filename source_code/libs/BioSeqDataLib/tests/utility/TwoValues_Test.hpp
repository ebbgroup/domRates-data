/*
 * TwoValues_Test.hpp
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


#ifndef TWOVALUES_TEST_HPP_
#define TWOVALUES_TEST_HPP_



#include <boost/test/unit_test.hpp>

#include "../../src/utility/TwoValues.hpp"


BOOST_AUTO_TEST_SUITE(TwoValues_Test)


BOOST_AUTO_TEST_CASE( TwoValues_Test)
{
	BioSeqDataLib::TwoValues val(3,4);
	BOOST_CHECK_EQUAL(val.first(), 3);
	BOOST_CHECK_EQUAL(val.second(), 4);

	val += BioSeqDataLib::TwoValues(2,6);
	BOOST_CHECK_EQUAL(val.first(), 5);
	BOOST_CHECK_EQUAL(val.second(), 10);
	BOOST_CHECK_EQUAL(val.value(), 0.5);

	BioSeqDataLib::TwoValues val2(val);
	BOOST_CHECK_EQUAL(val2.first(), 5);
	BOOST_CHECK_EQUAL(val2.second(), 10);

	BioSeqDataLib::TwoValues val3(std::move(val));
	BOOST_CHECK_EQUAL(val3.first(), 5);
	BOOST_CHECK_EQUAL(val3.second(), 10);

	val3.first(4);
	val2=val3;
	BOOST_CHECK_EQUAL(val2.first(), 4);
	BOOST_CHECK_EQUAL(val2.second(), 10);

	val=std::move(val3);
	BOOST_CHECK_EQUAL(val.first(), 4);
	BOOST_CHECK_EQUAL(val.second(), 10);

	BioSeqDataLib::TwoValues val4;
	BOOST_CHECK_EQUAL(val4.first(), 0);
	BOOST_CHECK_EQUAL(val4.second(), 0);
}

BOOST_AUTO_TEST_SUITE_END()





#endif /* TWOVALUES_TEST_HPP_ */
