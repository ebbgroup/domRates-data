/*
 * DomainArrangementSet.hpp
 *
 *  Created on: 15 Nov 2013
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
 * \file DomainArrangementSet.hpp
 * \brief Header for the DomainArrangementSet
 */
#ifndef DOMAINARRANGEMENTSET_HPP_
#define DOMAINARRANGEMENTSET_HPP_

// C++ header
#include <algorithm>
#include <fstream>
#include <set>
#include <map>
#include <string>
#include <iostream>
#include <ios>
#include <iomanip>
#include <memory>
#include <type_traits>

// BioSeqDataLib header
#include "DomainArrangement.hpp"
#include "Domain.hpp"
#include "PfamDomain.hpp"
#include "SFDomain.hpp"

#include "../utility/stringHelpers.hpp"
#include "../utility/Exceptions.hpp"
#include "../utility/utility.hpp"
#include "../external/Input.hpp"


namespace AP=AlgorithmPack;

namespace BioSeqDataLib
{

typedef enum {unknown=-1, pfam=0, hmmscan_domtbl=1, xdom=2, ass=3, interpro_tsv=4, dama=5, radiant=6 } DomainFileFormat;


std::string
getFormatString(const DomainFileFormat &format);

/** \addtogroup DomainGroup
 *  @{
 */

/**
 * \brief Class to represent a set of domain arrangements.
 *
 * \tparam DomainType The domain type to use (e.g. PfamDomain)
 */
template<typename DomainType>
class DomainArrangementSet
{
private:
	std::map<std::string, DomainArrangement<DomainType> > arrangements_; // DomainArrangement storage by sequence name
	std::map<std::string, int> uniqueIds_; // Stores the unique ids of every domain type.
	using DomainArrangement_t = typename std::map<std::string, DomainArrangement<DomainType> >;

	template<typename DomainType2>
	static void
	pfamTokens2Domain_(std::string &&acc __attribute__((unused)), const std::vector<std::string> &tokens __attribute__((unused)), DomainArrangement<DomainType2> &da __attribute__((unused)))
	{
		throw std::runtime_error("This domain combination is currently not supported");
	}
	static void
	pfamTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<Domain> &da);
	static void
	pfamTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<DomainExt> &da);
	static void
	pfamTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<PfamDomain> &da);

	template<typename DomainType2>
	static void
	radiantTokens2Domain_(std::string &&acc __attribute__((unused)), const std::vector<std::string> &tokens __attribute__((unused)), DomainArrangement<DomainType2> &da __attribute__((unused)))
	{
		throw std::runtime_error("This domain combination is currently not supported");
	}
	static void
	radiantTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<Domain> &da);
	//static void
	//radiantTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<DomainExt> &da);
	//static void
	//radiantTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<PfamDomain> &da);


	template<typename DomainType2>
	static void
	hmmTokens2Domain_(std::string &&acc __attribute__((unused)), const std::vector<std::string> &tokens __attribute__((unused)), DomainArrangement<DomainType2> &da __attribute__((unused)))
	{
		throw std::runtime_error("This domain combination is currently not supported");
	}
	static void
	hmmTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<Domain> &da);
	static void
	hmmTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<DomainExt> &da);
	static void
	hmmTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<PfamDomain> &da);
	static void
	hmmTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<SFDomain> &da);



	template<typename DomainType2>
	static void
	assTokens2Domain_(const std::vector<std::string> &tokens __attribute__((unused)), DomainArrangement<DomainType2> &da __attribute__((unused)))
	{
		throw std::runtime_error("This domain combination is currently not supported");
	}
	static void
	assTokens2Domain_(const std::vector<std::string> &tokens, DomainArrangement<Domain> &da);
	static void
	assTokens2Domain_(const std::vector<std::string> &tokens, DomainArrangement<DomainExt> &da);
	static void
	assTokens2Domain_(const std::vector<std::string> &tokens, DomainArrangement<SFDomain> &da);



	template<typename DomainType2>
	static void
	interProTSVTokens2Domain_(const std::vector<std::string> &tokens __attribute__((unused)), DomainArrangement<DomainType2> &da __attribute__((unused)))
	{
		throw std::runtime_error("This domain combination is currently not supported");
	}
	static void
	interProTSVTokens2Domain_(const std::vector<std::string> &tokens, DomainArrangement<Domain> &da);

	template<typename DomainType2>
	static void
	damaTokens2Domain_(const std::vector<std::string> &tokens __attribute__((unused)), DomainArrangement<DomainType2> &da __attribute__((unused)))
	{
		throw std::runtime_error("This domain combination is currently not supported");
	}
	static void
	damaTokens2Domain_(const std::vector<std::string> &tokens, DomainArrangement<Domain> &da);

	static
	DomainFileFormat
	_identifyFormat(AP::Input &inFile);

	void
	_readPfamScanOutput(AP::Input &inFile);

	void
	_readRadiantOutput(AP::Input &inFile);

	void
	_readHMMScanDomtblOutput(AP::Input &inFile);

	void
	_readXDOMFormat(AP::Input &inFile);

	void
	_readInterProTSV(AP::Input &inFile);

	void
	_readDAMAOutput(AP::Input &inFile);

	// empty function, only exist to make specialization for cases it makes sense
	void
	_readASSFile(AP::Input &inFile);




	// empty function, only exist to make specialization for cases it makes sense
	void
	_writePfamScanOutput(std::ofstream &outFile __attribute__((unused)))
	{}


	template<typename SeqSet>
	void
	_writeXDomFormat(std::ofstream &outFile, const SeqSet &seqSet) const;

	void
	_writeXDomFormat(std::ofstream &outFile) const;

public:
	using iterator = typename DomainArrangement_t::iterator;
	using const_iterator = typename DomainArrangement_t::const_iterator;
	using reverse_iterator = typename DomainArrangement_t::reverse_iterator;
	using const_reverse_iterator = typename DomainArrangement_t::const_reverse_iterator;


	typedef DomainType value_type;

	/**
	 * \brief Standard constructor
	 */
	DomainArrangementSet();

	/**
	 * \brief Move constructor
	 * @param The domain arrangement set to move.
	 */
	DomainArrangementSet(DomainArrangementSet<DomainType> &&) = default;


	/**
	 * \brief Standard destructor
	 */
	virtual ~DomainArrangementSet();


	/**
	 * \brief The move assignment operator
	 * @param The DomainArrangementSet to move
	 * @return Reference to the DomainArrangementSet.
	 */
	DomainArrangementSet& operator=(DomainArrangementSet<DomainType>&&) = default;

	/**
	 * \brief Access the domain arrangement of a certain sequence.
	 * @param sourceId The sequence name.
	 * @return Reference to the domain arrangement.
	 */
	DomainArrangement<DomainType>
	&operator[](const std::string &sourceId)
	{
		return arrangements_[sourceId];
	}


	 /**
	 * \brief Access the domain arrangement of a certain sequence.
	 * @param sourceId The sequence name.
	 * @return Reference to the domain arrangement.
	 */
	const DomainArrangement<DomainType>
	&operator[](const std::string &sourceId) const
	{
		return arrangements_.find(sourceId)->second;
	}

	/**
	 * \brief Reads a domain file.
	 *
	 * In case of PFAM-B domains, the significance is set to -1, the clan will be empty just like in PFAM-A domains that do not belong to a clan.
	 * @param in_f The input file.
	 *
	 */
	void
	read(const fs::path &in_f);

	/**
	 * \brief Writes arrangement set to a file.
	 * @param outFileName The output file
	 * @param format The format to use.
	 */
	template<typename SeqSet>
	void
	write(const std::string &outFileName, const std::string &format, const SeqSet &seqSet);

	void
	write(const std::string &outFileName, const std::string &format);




	/**
	 * \brief Deletes elements from the DomainArrangementSet
	 * @param args The arguments.
	 */
	template<class... Args>
	void
	erase(Args&&... args)
	{
		arrangements_.erase(std::forward<Args>(args)...);
	}

	/**
	 * \brief Returns the number of arrangements.
	 * @return The number of arrangements.
	 */
	size_t
	size() const
	{
		return arrangements_.size();
	}

	template<typename... Args>
	std::pair<iterator, bool>
	emplace(Args&& ... args)
	{
		return arrangements_.emplace(std::forward<Args>(args)...);
	}

	std::pair<iterator, bool>
	insert(const std::pair<std::string, DomainArrangement<DomainType> > &da)
	{
		 return arrangements_.insert(da);
	}


	/**
	 * \brief Search for the arrangement of a certain sequence.
	 * @param key The sequence id.
	 * @return Iterator to an element.
	 */
	iterator
	find(const std::string &key)
	{
		return arrangements_.find(key);
	}


	/**
	 * \brief Search for the arrangement of a certain sequence.
	 * @param key The sequence id.
	 * @return Iterator to an element.
	 */
	const_iterator
	find(const std::string &key) const
	{
		return arrangements_.find(key);
	}



	/**
	 * \brief Moves all arrangement to the belonging sequence.
	 * \tparam The type of sequence set. The sequences stores need to have DomainArrangementWrapper as template.
	 * @param seqSet The sequence set.
	 */
	template<typename SeqSetType>
	void
	moveTo(SeqSetType &seqSet)
	{
		size_t nSeqs = seqSet.size();
		auto it = arrangements_.end();
		auto itEnd = it;
		for (size_t i=0; i<nSeqs; ++i)
		{
			if ((it=arrangements_.find(seqSet[i].name())) != itEnd)
				seqSet[i].domainArrangement(std::move(it->second));
		}
		arrangements_.clear();
	}

	/**
	 * \brief Returns a list of all arrangements containing the given order of domains.
	 * @param ids The domain accession ids.
	 * @return The list of arrangements.
	 */
	std::vector<std::string>
	search(const std::vector<std::string> &ids) const;


	/**
	 * \brief Returns a list of all arrangements containing the given domains.
	 * @param ids The domain accession ids.
	 * @return The list of arrangements.
	 */
	std::vector<std::string>
	contain(const std::vector<std::string> &ids) const;


	void
	solveDbOverlaps(std::vector<DomainDB> dbImportance, size_t maxAbsOverlap, float maxFracOverlap)
	{
		for (auto &da : arrangements_)
			da.second.solveDbOverlaps(dbImportance, maxAbsOverlap, maxFracOverlap);
	}


	/********************************************************************
	 *                ITERATOR FUNCTIONS                                *
	 ********************************************************************/

	/**
	 * \brief Gets an iterator to the begin of the set.
	 * @return An iterator to the begin of the set.
	 */
	iterator
	begin()
	{
		return arrangements_.begin();
	}


	/**
	 * \brief Gets an iterator to the begin of the set.
	 * @return An iterator to the begin of the set.
	 */
	const_iterator
	begin() const
	{
		return arrangements_.begin();
	}

	/**
	 * \brief Gets an iterator to the end of the set.
	 * @return An iterator to the end of the set.
	 */
	iterator
	end()
	{
		return arrangements_.end();
	}


	/**
	 * \brief Gets an iterator to the end of the set.
	 * @return An iterator to the end of the set.
	 */
	const_iterator
	end() const
	{
		return arrangements_.end();
	}

	/**
	 * \brief Returns an iterator to the reverse beginning.
	 * @return Iterator to the reverse beginning.
	 */
	reverse_iterator
	rbegin()
	{
		return arrangements_.rbegin();
	}

	/**
	 * \brief Gets an iterator to the reverse beginning of the set.
	 * @return An iterator to the reverse begin of the set.
	 */
	const_reverse_iterator
	rbegin() const
	{
		return arrangements_.rbegin();
	}

	/**
	 * \brief Gets an iterator to the reverse end of the set.
	 * @return An iterator to the reverse end of the set.
	 */
	reverse_iterator
	rend()
	{
		return arrangements_.rend();
	}


	/**
	 * \brief Gets an reverse iterator to the reverse end
	 * @return An iterator to the reverse end of the set.
	 */
	const_reverse_iterator
	rend() const
	{
		return arrangements_.rend();
	}

};


/* Template overloads for specific domain types */
/*
template<>
void
DomainArrangementSet<PfamDomain>::_writePfamScanOutput(std::ofstream &outFile);

template<>
void
DomainArrangementSet<PfamDomain>::_readPfamScanOutput(AP::Input &inFile);


template<>
void
DomainArrangementSet<SFDomain>::_readASSFile(AP::Input &inFile);

*/

template<typename DomainType>
DomainArrangementSet<DomainType>::DomainArrangementSet() : arrangements_(), uniqueIds_()
{}

template<typename DomainType>
DomainArrangementSet<DomainType>::~DomainArrangementSet()
{}


template<typename DomainType>
void
DomainArrangementSet<DomainType>::read(const fs::path &inFileName)
{
	DomainFileFormat format = unknown;
	AP::Input inFile;
	inFile.open(inFileName);
	try
	{
		format = _identifyFormat(inFile);
		inFile.seekg(0, std::ios_base::beg);
		switch ( format )
		{
			case pfam:
				_readPfamScanOutput(inFile);
				break;
			case hmmscan_domtbl:
				_readHMMScanDomtblOutput(inFile);
				break;
			case xdom:
				_readXDOMFormat(inFile);
				break;
			case ass:
				_readASSFile(inFile);
				break;
			case interpro_tsv:
				_readInterProTSV(inFile);
				break;
			case dama:
				_readDAMAOutput(inFile);
				break;
			case radiant:
				_readRadiantOutput(inFile);
				break;
			default:
				throw FormatException("Unknown Format in " + inFileName.string());
		}
		// make sure they are all ordered the correct way
		for (auto &arrangement : arrangements_)
			sort(arrangement.second.begin(), arrangement.second.end());

	}
	catch (const FormatException &formatException)
	{
		throw;
	}
	catch (const std::exception &e)
	{
		throw FormatException("Failed to read file: " + inFileName.string() + ". Format identified as: " + getFormatString(format));
	}

	for (auto &da : this->arrangements_)
	{
		std::sort(da.second.begin(), da.second.end(), [](const DomainType & a, const DomainType & b) -> bool
		{
		    return a.start() < b.start();
		});
	}
}



template<>
void
DomainArrangementSet<PfamDomain>::_writePfamScanOutput(std::ofstream &outFile);


template<typename DomainType>
template<typename SeqSet>
void
DomainArrangementSet<DomainType>::write(const std::string &outFileName, const std::string &format, const SeqSet &seqSet)
{
	try
	{
		std::ofstream outFile;
		outFile.exceptions (std::ofstream::failbit | std::ofstream::badbit );
		outFile.open(outFileName);
		outFile.exceptions ( std::ifstream::badbit);
		if (format == "pfam")
			_writePfamScanOutput(outFile);
		else
		{
			if (format == "xdom")
				_writeXDomFormat(outFile, seqSet);
			else
				throw FormatException("Unknown format " + format);
		}
	}
	catch (const std::exception &e)
	{
		throw FormatException("Failed to write file '" + outFileName + "' of format " + format);
	}
}

template<typename DomainType>
void
DomainArrangementSet<DomainType>::write(const std::string &outFileName, const std::string &format)
{
	try
	{
		std::ofstream outFile;
		outFile.exceptions (std::ofstream::failbit | std::ofstream::badbit );
		outFile.open(outFileName);
		outFile.exceptions ( std::ifstream::badbit);
		if (format == "pfam")
			_writePfamScanOutput(outFile);
		else
		{
			if (format == "xdom")
				_writeXDomFormat(outFile);
			else
				throw FormatException("Unknown format " + format);
		}
	}
	catch (const std::exception &e)
	{
		throw FormatException("Failed to write file '" + outFileName + "' of format " + format);
	}
}




template<typename DomainType>
template<typename SeqSet>
void
DomainArrangementSet<DomainType>::_writeXDomFormat(std::ofstream &outFile, const SeqSet &seqSet) const
{

	for (auto &pair : arrangements_)
	{
		outFile << ">" << pair.first << " " << seqSet[pair.first].length() << "\n";
		for (auto &dom : pair.second)
		{
			outFile << dom.start()+1 << " " << dom.end()+1 << " " << dom.accession() << " " << dom.evalue() << "\n";
		}
	}
}

template<typename DomainType>
void
DomainArrangementSet<DomainType>::_writeXDomFormat(std::ofstream &outFile) const
{

	for (auto &pair : arrangements_)
	{
		outFile << ">" << pair.first << "\n";
		for (auto &dom : pair.second)
		{
			outFile << dom.start()+1 << " " << dom.end()+1 << " " << dom.accession() << " " << dom.evalue() << "\n";
		}
	}
}



template<typename DomainType>
DomainFileFormat
DomainArrangementSet<DomainType>::_identifyFormat(AP::Input &inFile)
{
	std::string line;
	getline(inFile, line);
	// pfam format
	if (line.find("pfam_scan.pl") != std::string::npos)
		return DomainFileFormat::pfam;

	// hmmscan output
	if (line.find("--- full sequence --- -------------- this domain ------------") != std::string::npos)
		return DomainFileFormat::hmmscan_domtbl;

	// xdom format
	if (line[0] == '>')
		return DomainFileFormat::xdom;

	auto tokens = split(line, "\t", true);

	// superfamily output
	if (tokens.size() == 9)
		return DomainFileFormat::ass;

	// interpro scan TSV format
	if (tokens.size() >= 11)
		return DomainFileFormat::interpro_tsv;

	// DAMA output
	if (tokens.size() == 6)
		return DomainFileFormat::dama;

	// RADIANT format
	if (line.find("RADIANT") != std::string::npos)
		return DomainFileFormat::radiant;
	// unknown format
	return DomainFileFormat::unknown;
}


template<typename DomainType>
void
DomainArrangementSet<DomainType>::_readPfamScanOutput(AP::Input &inFile)
{
	/*
	 * pfam_scan.pl format
	 * # <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>
	 *1kjz_A          5    193      3    194 PF00009.22  GTP_EFTU          Domain     3   187   188    145.5   1.2e-42   1 CL0023
	 */
	std::string line;
	std::vector<std::string> domainInfo;
	std::string accessionNumber;
	std::string::size_type pos;
	auto it = uniqueIds_.end();
	int id;
	int line_counter_error = 0;
	while (getline(inFile, line))
	{
		++line_counter_error;
		if (line.empty() || (line[0] == '#'))
			continue;

		domainInfo.clear();
		split(line, " ", domainInfo);
		if (domainInfo.size() != 15)
			throw FormatException("Error! File not in proper PfamFormat. Error occured in line: " + std::to_string(line_counter_error) );
		if ((pos=domainInfo[5].find('.')) == std::string::npos)
			accessionNumber=domainInfo[5];
		else
			accessionNumber=domainInfo[5].substr(0, pos);
		if ((it=uniqueIds_.find(accessionNumber)) == uniqueIds_.end())
		{
			id = uniqueIds_.size();
			uniqueIds_[accessionNumber] = id;
		}
		else
			id= it->second;
		pfamTokens2Domain_(std::move(accessionNumber), domainInfo, arrangements_[domainInfo[0]]);//arrangements_[domainInfo[0]]);
	}
}


template<typename DomainType>
void
DomainArrangementSet<DomainType>::_readRadiantOutput(AP::Input &inFile)
{

/*
	# <seq id> <match start> <match end> <hmm acc> <hmm name> <type> <clan>

	ENSP00000362111.4      25     25 PF00335 Tetraspannin    Family      CL0347
*/

	std::string line;
	std::vector<std::string> domainInfo;
	std::string accessionNumber;
	auto it = uniqueIds_.end();
	int id;
	while (getline(inFile, line))
	{
		if (line.empty() || (line[0] == '#'))
			continue;

		domainInfo.clear();
		split(line, " ", domainInfo);
		accessionNumber=domainInfo[3];
		if ((it=uniqueIds_.find(accessionNumber)) == uniqueIds_.end())
		{
			id = uniqueIds_.size();
			uniqueIds_[accessionNumber] = id;
		}
		else
			id= it->second;
		radiantTokens2Domain_(std::move(accessionNumber), domainInfo, arrangements_[domainInfo[0]]);//arrangements_[domainInfo[0]]);
	}
}


template<typename DomainType>
void
DomainArrangementSet<DomainType>::_readDAMAOutput(AP::Input &inFile)
{
	/*
	 * 2.8e-26	404	508	sp|P20111|ACTN2_CHICK	PF00435	1
	 */
	std::string line;
	std::vector<std::string> domainInfo;
	auto it = uniqueIds_.end();
	int id;
	while (getline(inFile, line))
	{
		domainInfo.clear();
		split(line, "\t", domainInfo);
		if ((it=uniqueIds_.find(domainInfo[4])) == uniqueIds_.end())
		{
			id = uniqueIds_.size();
			uniqueIds_[domainInfo[4]] = id;
		}
		else
			id= it->second;
		damaTokens2Domain_(domainInfo, arrangements_[domainInfo[3]]);//arrangements_[domainInfo[0]]);
	}

	for (auto &arrangement : arrangements_)
		sort(arrangement.second.begin(), arrangement.second.end());
}


template<typename DomainType>
void
DomainArrangementSet<DomainType>::_readInterProTSV(AP::Input &inFile)
{
	//EFX70648	2c281c89dcf49d2e45ab9cf625a5ac16	1040	PRINTS	PR00419	Adrenodoxin reductase family signature	337	351	1.299999245E-14	T	24-02-2015

	std::string line;
	std::vector<std::string> domainInfo;
	std::string accessionNumber;
	auto it = uniqueIds_.end();
	int id;
	while (getline(inFile, line))
	{
		if (line.empty())
			continue;

		domainInfo.clear();
		split(line, "\t", domainInfo, true);
		accessionNumber=domainInfo[4];
		if ((it=uniqueIds_.find(accessionNumber)) == uniqueIds_.end())
		{
			id = uniqueIds_.size();
			uniqueIds_[accessionNumber] = id;
		}
		else
			id= it->second;
		interProTSVTokens2Domain_(domainInfo, arrangements_[domainInfo[0]]);
	}
}


template<typename DomainType>
void
DomainArrangementSet<DomainType>::_readHMMScanDomtblOutput(AP::Input &inFile)
{
	/*
	 * hmmscan tblout format
	 * #                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
	 * # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
	 * # ------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
	 * # GTP_EFTU             PF00009.22   188 1kjz_A               -            400   5.9e-43  146.4   0.1   1   1   1.8e-45   1.2e-42  145.5   0.1     3   187     5   193     3   194 0.95 Elongation factor Tu GTP binding domain
	 */

	std::string line;
	std::vector<std::string> domainInfo;
	std::string accessionNumber;
	std::string::size_type pos;
	auto it = uniqueIds_.end();
	int id;
	while (getline(inFile, line))
	{
		if (line.empty() || (line[0] == '#'))
			continue;
		domainInfo.clear();
		split(line, " ", domainInfo);
		if ((pos=domainInfo[1].find('.')) == std::string::npos)
			accessionNumber=domainInfo[1];
		else
			accessionNumber=domainInfo[1].substr(0, pos);
		if ((it=uniqueIds_.find(accessionNumber)) == uniqueIds_.end())
		{
			id = uniqueIds_.size();
			uniqueIds_[accessionNumber] = id;
		}
		else
			id= it->second;
		//arrangements_[domainInfo[3]].emplace_back(accessionNumber, domainInfo[0], std::stoi(domainInfo[17])-1, std::stoi(domainInfo[18])-1, std::stoi(domainInfo[19])-1, std::stoi(domainInfo[20])-1, std::stoi(domainInfo[15])-1, std::stoi(domainInfo[16])-1, std::stoi(domainInfo[2]), std::stod(domainInfo[7]), std::stod(domainInfo[6]), id);
		hmmTokens2Domain_(std::move(accessionNumber), domainInfo, arrangements_[domainInfo[3]]);
	}
}

template<typename DomainType>
void
DomainArrangementSet<DomainType>::_readASSFile(AP::Input &inFile)
{
	/*
	   Sequence ID
	   SUPERFAMILY model ID
	   Match region
	   Evalue score
	   Model match start position
	   Alignment to model
	   Family evalue
	   SCOP Family ID
	   SCOP domain ID of closest structure (px value)
	   gnl|Cobs_1.4|Cobs_00004-mRNA-1	0045816	390-451	5.72e-10	1	ICPYECPTNLVYQACGSsCKETCDTVNNSDnlKCTTSTPIEGCFCPEDYIFY-NGSCILKKNC	0.022	44897	57568
	*/

	std::string line;
	std::vector<std::string> domainInfo;
// 	std::string::size_type pos;
	auto it = uniqueIds_.end();
	int id;
	while (getline(inFile, line))
	{
		if (line.empty() || (line[0] == '#'))
			continue;
		domainInfo.clear();
		split(line, "\t", domainInfo);

		if ((it=uniqueIds_.find(domainInfo[8])) == uniqueIds_.end())
		{
			id = uniqueIds_.size();
			uniqueIds_[domainInfo[8]] = id;
		}
		else
			id= it->second;

		if ((it=uniqueIds_.find(domainInfo[8])) == uniqueIds_.end())
		{
			id = uniqueIds_.size();
			uniqueIds_[domainInfo[8]] = id;
		}
		else
			id= it->second;
		assTokens2Domain_( domainInfo, arrangements_[domainInfo[0]]);
	}
}

template<typename DomainType>
void
DomainArrangementSet<DomainType>::_readXDOMFormat(AP::Input &inFile)
{
	std::string line;
	auto it = arrangements_.end();
	while (getline(inFile, line))
	{
		if (line.empty())
			continue;
		if (line[0]=='>')
		{
			auto tokens = split(line, " >");
			auto pair = arrangements_.emplace(tokens[0], DomainArrangement<DomainType>());
			it = pair.first;
		}
		else
		{
			auto tokens = split(line, " ");
			it->second.emplace_back(tokens[2], stoul(tokens[0])-1, stoul(tokens[1])-1, stod(tokens[3]));
		}
	}
}



template<typename DomainType>
std::vector<std::string>
DomainArrangementSet<DomainType>::search(const std::vector<std::string> &ids) const
{
	std::vector<std::string> seqNames;
	size_t nDomains = ids.size();
	for (auto &pair : arrangements_)
	{
		if (pair.second.size() == nDomains)
		{
			size_t i;
			for (i=0; i<nDomains; ++i)
			{
				if (ids[i] != pair.second[i].accession())
					break;
			}
			if (i==nDomains)
				seqNames.push_back(pair.first);
		}
	}
	return seqNames;
}


template<typename DomainType>
std::vector<std::string>
DomainArrangementSet<DomainType>::contain(const std::vector<std::string> &ids) const
{
	std::vector<std::string> seqNames;
	size_t nDomains = ids.size();
	for (auto &pair : arrangements_)
	{
		std::map<std::string, int> counts = domainCounts(pair.second);
		auto itEnd = counts.end();
		size_t i;
		for (i=0; i<nDomains; ++i)
		{
			if (counts.find(ids[i]) == itEnd)
				break;
		}
		if (i==nDomains)
			seqNames.push_back(pair.first);
	}
	return seqNames;
}

/**
 * \brief Calculates the domain types occurring in the set.
 * \relates DomainArrangementSet
 * @param set The set to analyse.
 * @return A set containing the accession names of the different domains.
 */
template<typename DomainType>
std::set<std::string>
types(const DomainArrangementSet<DomainType> &set)
{
	std::set<std::string> dTypes;
	for (const auto &DA : set)
	{
		for (const auto &domain : DA.second)
			dTypes.insert(domain.accession());
	}
	return dTypes;
}


/**
 * \brief Counts the occurrence of each domain an a domain arrangement set.
 * @return The domain and the number of their occurrence.
 * \tparam DomainType The domain type
 * \relates DomainArrangementSet
 */
template<typename DomainType>
std::map<std::string, int>
domainCounts(const DomainArrangementSet<DomainType> &daSet)
{
	std::map<std::string, int> summary;
	for (const auto da : daSet)
	{
		for (const DomainType &domain : da.second)
			summary[domain.accession()] += 1;
	}
	return summary;
}


/**
 *
 * \brief Counts the occurrence of each clan in a domainArrangmentSet.
 * \param If useDomain is true, the domain accession number is used if the domain is not assigned to any clan.
 * @return The domain and the number of their occurrence.
 * \relates DomainArrangementSet
 */
std::map<std::string, int>
clanCounts(const DomainArrangementSet<PfamDomain> &daSet, bool useDomain=true);


/**
 *
 *\brief Function to reconstruct Split Domains into whole Domains
 *@param DomainSet The DomainArrangementSet
 *@param ot The overlap treshhold for domains
 *@param maxdist The maximal distance between two parts of a split domains.
 * \relates DomainArrangementSet
 */
template<typename DomainType>
void
splitDomRec(DomainArrangementSet<DomainType> &domainSet, unsigned int ot, unsigned int maxdist)
{
	typedef typename std::map <std::string,DomainArrangement<DomainType> >::iterator ItType;
	for (ItType itsplitdoma = domainSet.begin(); itsplitdoma != domainSet.end(); ++itsplitdoma)
		splitDomRec(itsplitdoma->second, ot, maxdist);
}
/** @} */ // Domain module


} /* namespace BioSeqDataLib */

#endif /* DOMAINARRANGEMENTSET_HPP_ */
