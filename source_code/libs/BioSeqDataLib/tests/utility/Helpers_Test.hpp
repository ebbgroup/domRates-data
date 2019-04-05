/*
 * Helpers_Test.hpp
 *
 *  Created on: 25 Feb 2014
 *      Author: Carsten Kemena
 *		 Email: c.kemena[@]uni-muenster.de
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


/**
 * \file Helpers_Test.hpp
 * \brief .
 */


#include "../../src/utility/stringHelpers.hpp"


#ifndef HELPERS_TEST_HPP_
#define HELPERS_TEST_HPP_

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Helpers_Test)


BOOST_AUTO_TEST_CASE( Helpers_Test)
{
	std::string line= "This is  a Test-to check.";
	std::vector<std::string> tokens = BioSeqDataLib::split(line, " -");
	BOOST_CHECK_EQUAL(tokens.size(), 6);
	BOOST_CHECK_EQUAL(tokens[0], "This");
	BOOST_CHECK_EQUAL(tokens[2], "a");
	BOOST_CHECK_EQUAL(tokens[3], "Test");

	line= "This is  a Test-to check. ";
	tokens = BioSeqDataLib::split(line, " -");
	BOOST_CHECK_EQUAL(tokens.size(), 6);
	BOOST_CHECK_EQUAL(tokens[5], "check.");

	line= "This is  a Test-to check.";
	tokens.clear();
	BioSeqDataLib::split(line, " -", tokens);
	BOOST_CHECK_EQUAL(tokens.size(), 6);
	BOOST_CHECK_EQUAL(tokens[0], "This");
	BOOST_CHECK_EQUAL(tokens[2], "a");
	BOOST_CHECK_EQUAL(tokens[3], "Test");

	line= "This is  a Test-to check. ";
	tokens.clear();
	BioSeqDataLib::split(line, " -", tokens);
	BOOST_CHECK_EQUAL(tokens.size(), 6);


	line= " test \t ";
	BioSeqDataLib::trimRight(line);
	BOOST_CHECK_EQUAL(line, " test");

	line= " test";
	BioSeqDataLib::trimRight(line);
	BOOST_CHECK_EQUAL(line, " test");

	line= " test \t ";
	BioSeqDataLib::removeSpaces(line);
	BOOST_CHECK_EQUAL(line, "test");

	line= "test";
	BioSeqDataLib::removeSpaces(line);
	BOOST_CHECK_EQUAL(line, "test");

	line= " te st -";
	BioSeqDataLib::strip(line, 0);
	BOOST_CHECK_EQUAL(line, "test");
}


BOOST_AUTO_TEST_SUITE_END()

#endif /* HELPERS_TEST_HPP_ */
