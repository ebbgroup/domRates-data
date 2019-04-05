/*
 * DSM.hpp
 *
 *  Created on: Nov 23, 2015
 *      Author: ckeme_01
 */

/**
 * \file DSM.hpp
 * \brief File containing the domain similarity class.
 */
#ifndef SRC_UTILITY_DSM_HPP_
#define SRC_UTILITY_DSM_HPP_

// C++ header
#include <map>
#include <string>
#include <vector>
#include <fstream>

// Boost header
#include <boost/filesystem.hpp>

// BSDL header
#include "stringHelpers.hpp"
#include "../domain/Domain.hpp"


namespace fs = boost::filesystem;

namespace BioSeqDataLib
{

/**
 * \brief Class representing a domain similarity matrix (DSM).
 *
 * The domain matrix is stored in CRS format (Compressed Row Storage) to as most values in the file will be around zero. Since the matrix is symmetric only half of it is stored.
 */

class DSM
{
private:
	std::string description_;
	int nDomains_;
	std::map<std::string, int> acc2id_; // stores mapping from name (e.g. PF00039) to the row/column position (e.g. 5).
	std::vector<short> values_; // The values stored in the matrix.
	std::vector<int> colIDs_; // The column id.
	std::vector<int> rowIDs_; // The row id.
	short threshold_;
	bool useNegative_;

public:
	/**
	 * \brief Standard constructor
	 */
	DSM();

	/**
	 * \brief Constructor reading matrix from file.
	 * @param inFile The input file.
	 */
	explicit DSM(const fs::path&inFile) : useNegative_(false)
	{
		read(inFile);
	}

	/**
	 * \brief Standard destructor
	 */
	virtual ~DSM();

	short
	val(const std::string &d1, const std::string &d2) const;

	/**
	 * \brief Returns the matching value of two domains.
	 * @param d1 The first domain.
	 * @param d2 The second domain.
	 * @return The matching value if stored in matrix, else -1.
	 */

	short
	val(const Domain &d1, const Domain &d2) const
	{
		return val(d1.accession(), d2.accession());
	}



	/**
	 * \brief Reads matrix from file
	 * @param inFile The file containing the matrix.
	 */
	void
	read(const fs::path &inFile);

	/*
	 * \brief Sets the return on a given value to the given value.
	 * @param val The new return value.
	 *
	void
	setNonExistent(short val)
	{
		nonExistent_ = val;
	}*/

	/**
	 * \brief Can torn on scaling of values between -100 and 100.
	 * The nonExistent value is not affected.
	 * @param true turns on scaling, false turns it off.
	 */
	void
	useNegative(bool val)
	{
		useNegative_ = val;
	}

	/**
	 * \brief Returns threshold used when storing the matrix.
	 */
	short
	threshold() const
	{
		return threshold_;
	}

};






} /* namespace BioSeqDataLib */

#endif /* SRC_UTILITY_DSM_HPP_ */
