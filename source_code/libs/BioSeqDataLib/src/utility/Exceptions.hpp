/*
 * Exceptions.hpp
 *
 *  Created on: 22 Nov 2013
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

/**
 * \file Exceptions.hpp
 * \brief File containing several exception classes.
 */
#include <stdexcept>

#ifndef EXCEPTIONS_HPP_
#define EXCEPTIONS_HPP_


namespace BioSeqDataLib
{

/**
 * \brief Exception to throw when a problem with the format occurred.
 */
class FormatException : public std::runtime_error
{
public:
	/**
	 * \brief Message to return with what()
	 * @param message The message.
	 */

	explicit FormatException(const char *message) : std::runtime_error(message)
	{}

	explicit FormatException(const std::string &message) : std::runtime_error(message)
	{}


	FormatException(const FormatException&) = default;
	explicit FormatException(FormatException&&) = default;
	virtual ~FormatException() noexcept = default;
};

}


#endif /* EXCEPTIONS_HPP_ */
