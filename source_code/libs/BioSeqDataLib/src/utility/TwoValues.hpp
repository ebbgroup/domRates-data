/*
 * TwoValues.hpp
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

/**
 * \file TwoValues.hpp
 * \brief File containing the TwoValues class.
 */
#ifndef TWOVALUES_H_
#define TWOVALUES_H_

#include <cstdlib>

namespace BioSeqDataLib
{

/**
 * \brief Simple class to store to values.
 */
class TwoValues
{
private:
	size_t first_;
	size_t second_;

public:
	/**
	 * \brief Standard constructor
	 */
	TwoValues();

	/**
	 * \brief Constructor
	 * @param first The first value.
	 * @param second The second value.
	 */
	TwoValues(size_t first, size_t second);

	/**
	 * \brief Copy constructor
	 * @param Element to copy.
	 */
	TwoValues(const TwoValues &) = default;

	/**
	 * \brief Move constructor
	 * @param Element to move.
	 */
	TwoValues(TwoValues &&) = default;

	/**
	 * \brief Standard destructor
	 */
	virtual ~TwoValues();

	/**
	 * \brief Assignment operator
	 * @param The new value
	 * @return The object
	 */
	TwoValues & operator= ( const TwoValues & ) = default;

	/**
	 * \brief Move assignment operator
	 * @param The new value
	 * @return The object
	 */
	TwoValues & operator= ( TwoValues && ) = default;

	/**
	 * \brief Adds a value to the TwoValues
	 * \warning This not a fraction sum
	 * @param rhs The value to add
	 * @return The new value
	 */
	TwoValues& operator+=(const TwoValues& rhs)
	{
		first_ += rhs.first_;
		second_ += rhs.second_;
		return *this;
	}

	/**
	 * \brief Sums two values of type TwoValues.
	 * \warning This not a fraction sum
	 * @param lhs The first value.
	 * @param rhs The second value
	 * @return The new value.
	 */
	friend TwoValues operator+(TwoValues lhs, const TwoValues& rhs)
	{
		return lhs += rhs;
	}

	/**
	 * \brief Comparison operator
	 * @param lhs Left value
	 * @param rhs Right values
	 * @return true if the two values are the same, false else
	 */
	friend bool operator==(const TwoValues &lhs, const TwoValues& rhs)
	{
		return ((lhs.first_==rhs.first_) && (lhs.second_==rhs.second_));
	}


	/**
	 * \brief Returns the first value.
	 * @return The first value.
	 */
	size_t
	first() const
	{
		return first_;
	}

	/**
	 * \brief Sets a new value for the first element.
	 * @param newVal The new value.
	 */
	void
	first(size_t newVal)
	{
		first_ = newVal;
	}

	/**
	 * \brief Returns the second value.
	 * @return The second value.
	 */
	size_t
	second() const
	{
		return second_;
	}

	/**
	 * \brief Sets a new value for the second element.
	 * @param newVal The new value.
	 */
	void
	second(size_t newVal)
	{
		second_ = newVal;
	}

	/**
	 * \brief Divides the first value with the second.
	 * @return The result of dividing the first value by the second.
	 */
	double
	value()
	{
		return (first_*1.0)/second_;
	}


};

} /* namespace BioSeqDataLib */

#endif /* TWOVALUES_H_ */
