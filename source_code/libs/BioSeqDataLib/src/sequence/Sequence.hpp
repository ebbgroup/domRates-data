/*
 * Sequence.hpp
 *
 *  Created on: 13 Oct 2013
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
 * \file Sequence.hpp
 * \brief Contains the Sequence class.
 */

#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_


// C header
#include <cstdlib>
#include <iostream>

// C++ header
#include <string>
#include <type_traits>
#include <utility>

#include "Alphabet.hpp"

namespace BioSeqDataLib
{


/**
 * \brief Class to represent a biological sequence.
 * \details This is the most basic sequence class. It saves the name, a comment and the sequence itself. Additionally
 * it can contain an id.
 */
template<class... Args>
class Sequence : public Args...
{
private:
	std::string name_;
	std::string sequence_;
	std::string accession_;
	std::string comment_;
	size_t id_;

	using Sequence_t = std::string;
public:

	using iterator = Sequence_t::iterator;
	using const_iterator = Sequence_t::const_iterator;
	using reverse_iterator = Sequence_t::reverse_iterator;
	using const_reverse_iterator = Sequence_t::const_reverse_iterator;

	/**************************************************************************
	 *                   Constructors & Destructors                           *
	 **************************************************************************/


	/**
	 * \brief Standard constructor
	 */
	Sequence();


	/**
	 * \brief Constructor initialising all members.
	 * @param name The sequence name.
	 * @param seq The sequence itself.
	 * @param accession The accession number.
	 * @param comment The sequence comment.
	 * @param id The sequence id.
	 */
	template <typename T1, typename T2, typename T3, typename T4>
	Sequence(T1 &&name, T2 &&seq, T3 &&accession=std::string(""), T4 &&comment=std::string(""), size_t id=0) noexcept(std::is_nothrow_copy_constructible<T1>::value && std::is_nothrow_copy_constructible<T2>::value && std::is_nothrow_copy_constructible<T3>::value && std::is_nothrow_copy_constructible<T4>::value);

	/**
	 * \brief Constructor reserving space for the sequence.
	 * @param name The sequence name.
	 * @param length The size to reserve.
	 * @param accession The accession number.
	 * @param comment The sequence comment.
	 * @param id The sequence id.
	 */
	template <typename T1, typename T2, typename T3>
	Sequence(T1 &&name, size_t length, T2 &&accession="", T3 &&comment="", size_t id=0) noexcept(std::is_nothrow_copy_constructible<T1>::value && std::is_nothrow_copy_constructible<T2>::value && std::is_nothrow_copy_constructible<T3>::value);

	/**
	 * \brief Copy constructor.
	 * @param
	 */
	Sequence(const Sequence&) = default;

	/**
	 * \brief Copy Move constructor.
	 * @param  The Sequence to copy
	 */
	Sequence(Sequence &&) = default;

	/**
	 * \brief Standard destructor.
	 */
	virtual ~Sequence();


	/**************************************************************************
	 *                          Operators                                     *
	 **************************************************************************/

	/**
	 * \brief Assignment operator.
	 * @param Sequence to assign.
	 * @return The sequence.
	 */
	Sequence&
	operator=(const Sequence&) = default;

	/**
	 * \brief Move assignment operator.
	 * @param Sequence to assign.
	 * @return The sequence.
	 */
	Sequence&
	operator=(Sequence&&) = default;


	/**
	 * \brief Access operator.
	 * \param index The index of to access.
	 */
	char
	&operator[](size_t index)
	{
		return sequence_[index];
	}

	/**
	 * \overload
	 */
	const char
	&operator[](size_t index) const
	{
		return sequence_[index];
	}


	/**
	 * \brief Comparison operators.
	 * \param a The first sequence.
	 * \param b The second sequence.
	 * \returns true if the both sequences are the same. Name and comment are ignored.
	 */
	friend bool operator ==(const Sequence &a, const Sequence &b)
	{
		return(a.seq() == b.seq());
	}

	/**
	 * \brief Comparison operators.
	 * \param a The first sequence.
	 * \param b The second sequence.
	 * \returns true if the both sequences are the same. Name and comment are ignored.
	 */
	friend bool operator !=(const Sequence &a, const Sequence &b)
	{
		return(a.seq() != b.seq());
	}

	/**
	 * \brief Returns the lengths of the sequence.
	 * @return The length.
	 */
	size_t
	size() const noexcept
	{
		return sequence_.size();
	}


	/**
	 * \brief Returns the lengths of the sequence.
	 * @return The length.
	 */
	size_t
	length() const noexcept
	{
		return sequence_.size();
	}
	/**
	 * \brief Returns if a sequence has length 0.
	 * \param True, if sequence length is 0, else false.
	 */
	size_t
	empty() const noexcept
	{
		return sequence_.empty();
	}

	/**
	 * \brief Returns the name of the sequence.
	 * @return The name of the sequence.
	 */
	std::string
	name() const
	{
		return name_;
	}

	/**
	 * \brief Returns the name of the sequence.
	 * @return The name of the sequence.
	 */
	template<typename T>
	void
	name(T &&newname_)
	{
		name_=std::forward<T>(newname_);
	}

	/**
	 * \brief Returns the accession of the sequence.
	 * @return The accession.
	 */
	std::string accession() const
	{
		return accession_;
	}

	/**
	 * \brief Sets the accession number.
	 * @param newAccession_ The accession number.
	 */
	template<typename T>
	void
	accession(T &&newAccession_)
	{
		accession_=std::forward<T>(newAccession_);
	}

	/**
	 * \brief Returns the comment of the sequence.
	 * @return The comment.
	 */
	std::string comment() const
	{
		return comment_;
	}

	/**
	 * \brief Sets the comment.
	 * @param newComment_ The new comment.
	 */
	template<typename T>
	void
	comment(T &&newComment_)
	{
		comment_=std::forward<T>(newComment_);
	}

	/**
	 * \brief Returns the sequence itself.
	 * @return The sequence.
	 */
	std::string
	seq() const
	{
		return sequence_;
	}


	/**
	 * \brief Returns the sequence id.
	 * @return The sequence id.
	 */
	size_t
	id() const
	{
		return id_;
	}

	/**
	 * \brief Sets the id of the sequence.
	 * @param newid_ The new id.
	 */
	void
	id(size_t newid_)
	{
		id_=newid_;
	}


	/**
	 * \brief Appends a sequence to the existing one.
	 * @param segment The segment to attach.
	 */
	void
	append(const std::string &segment)
	{
		sequence_.append(segment);
	}

	/**
	 * \brief Appends a character to the existing one.
	 * @param c The character to attach.
	 */
	void
	append(char c)
	{
		sequence_.push_back(c);
	}

	/**
	 * \brief Resizes the sequence.
	 * @param newSize The new sequence size.
	 */
	void
	resize(size_t newSize)
	{
		sequence_.resize(newSize);
	}


	/**
	 * \brief Removes a character from the sequence.
	 * @param c The character that will be removed.
	 */
	void
	remove(const char c)
	{
		size_t pos = 0;
		size_t len = sequence_.size();
		for (size_t i=0; i<len; ++i)
		{
			if (sequence_[i] != c)
				sequence_[pos++] = sequence_[i];
		}
		sequence_.resize(pos);
	}

	/**************************************************************************
	 *                             Iterator                                   *
	 **************************************************************************/

 	/**
 	 * \brief Returns an iterator to the first element in the container.
 	 * @return An iterator to the first element in the container.
 	 */
	iterator begin() throw()
	{
		return sequence_.begin();
	}

	/**
	 * \brief Returns an iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	iterator end() throw()
	{
		return sequence_.end();
	}

 	/**
 	 * \brief Returns a const_iterator to the first element in the container.
 	 * @return An iterator to the first element in the container.
 	 */
	const_iterator begin() const throw()
	{
		return sequence_.begin();
	}


	/**
	 * \brief Returns a const_iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	const_iterator end() const throw()
	{
		return sequence_.end();
	}

	/**
	 * \brief Returns a reverse_iterator to the last element in the BinarySearchTree.
	 * @return An iterator to the last element.
	 */
	reverse_iterator rbegin() throw()
	{
		return sequence_.rbegin();
	}

	/**
	 * \brief Returns a const_reverse_iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	reverse_iterator rend() throw()
	{
		return sequence_.rend();
	}

 	/**
 	 * \brief Returns a const_iterator to the first element in the container.
 	 * @return An iterator to the first element in the container.
 	 */
	const_iterator cbegin() const throw()
	{
		return sequence_.cbegin();
	}


	/**
	 * \brief Returns a const_iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	const_iterator cend() const throw()
	{
		return sequence_.cend();
	}

 	/**
 	 * \brief Returns a const_iterator to the first element in the container.
 	 * @return An iterator to the first element in the container.
 	 */
	const_reverse_iterator crbegin() throw()
	{
		return sequence_.crbegin();
	}

	/**
	 * \brief Returns a const_iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	const_reverse_iterator crend() const throw()
	{
		return sequence_.crend();
	}

};

using BasicSeq = Sequence<>;

/******************************************************************************
 *                   Constructors & Destructors                               *
 ******************************************************************************/


template<class... Args>
Sequence<Args ...>::Sequence():name_(""), sequence_(""), accession_(""), comment_(""), id_(0)
{}

template<class... Args>
template <typename T1, typename T2, typename T3, typename T4>
Sequence<Args ...>::Sequence(T1 &&name, T2 &&seq, T3 &&accession, T4 &&comment, size_t id)  noexcept(std::is_nothrow_copy_constructible<T1>::value && std::is_nothrow_copy_constructible<T2>::value && std::is_nothrow_copy_constructible<T3>::value && std::is_nothrow_copy_constructible<T4>::value):name_(std::forward<T1>(name)), sequence_(std::forward<T2>(seq)), accession_(std::forward<T3>(accession)), comment_(std::forward<T4>(comment)), id_(id)
{}

template<class... Args>
template <typename T1, typename T2, typename T3>
Sequence<Args ...>::Sequence(T1 &&name, size_t length, T2 &&accession, T3 &&comment, size_t id)  noexcept(std::is_nothrow_copy_constructible<T1>::value && std::is_nothrow_copy_constructible<T2>::value && std::is_nothrow_copy_constructible<T3>::value):name_(std::forward<T1>(name)), sequence_(), accession_(std::forward<T2>(accession)), comment_(std::forward<T3>(comment)), id_(id)
{
	sequence_.reserve(length);
}

template<class... Args>
Sequence<Args ...>::~Sequence()
{}



/******************************************************************************
 *                          print function                                  *
 ******************************************************************************/
template<class... Args>
std::ostream& operator<< (std::ostream &out, const Sequence<Args ...> &seq)
{
	out << ">" << seq.name() << "\n";
	out << seq.seq();
	return out;
}

} //BioSeqDataLib

#endif /* SEQUENCE_HPP_ */
