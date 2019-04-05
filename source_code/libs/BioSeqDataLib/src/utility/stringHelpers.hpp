/*
 * stringHelpers.hpp
 *
 *  Created on: 27 Oct 2013
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

/**
 * \file stringHelpers.hpp
 * \brief File containing some helper functions.
 */
#ifndef HELPERS_HPP_
#define HELPERS_HPP_

#include <cctype>
#include <string>
#include <vector>
#include <algorithm>

namespace BioSeqDataLib
{

	/**
	 * \brief Splits a string according to the given delimiters.
	 * @param[in] s The string to split.
	 * @param[in] delimiters The delimiters to use.
	 * @return A vector of the resulting tokens.
	 */
	std::vector<std::string>
	split(const std::string &s, const std::string &delimiters, bool keepEmpty = false);

	/**
	 * \brief Splits a string according to the given delimiters.
	 * @param[in] s The string to split.
	 * @param[in] delimiters The delimiters to use.
	 * @param[out] tokens  A vector to store the resulting tokens.
	 */
	void
	split(const std::string &s, const std::string &delimiters, std::vector<std::string> &tokens, bool keepEmpty = false);

	/**
	 * \brief Removes all whitespaces from the end of a string.
	 * @param str The string to trim
	 */
	void
	trimRight(std::string &str);

	/**
	 * \brief Removes all spaces from a string.
	 * @param str The string.
	 * @param startPos The position to start from.
	 */
	void
	removeSpaces(std::string &str, int startPos=0);

	/**
	 * \brief Removes all non alpha from string.
	 * @param line The line to clean.
	 */
	void
	strip(std::string &line, size_t startPos=0);
}


#endif /* HELPERS_HPP_ */
