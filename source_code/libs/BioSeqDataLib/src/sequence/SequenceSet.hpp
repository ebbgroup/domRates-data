/*
 * SequenceSet.hpp
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
 * \file SequenceSet.hpp
 * \brief Contains the SequenceSet class
 */

// C header
#include <cctype>
#include <cstdlib>

// C++ header
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "boost/filesystem.hpp"

#include "Sequence.hpp"
#include "SeqFunctions.hpp"

#include "../utility/stringHelpers.hpp"
#include "../utility/utility.hpp"
#include "../external/Input.hpp"

#ifndef SEQUENCESET_HPP_
#define SEQUENCESET_HPP_


namespace fs=boost::filesystem;

namespace BioSeqDataLib
{

/**
 * \brief Class to handle a sequence set.
 * \tparam SequenceType The sequence type to use.
 *
 * \details This class handles different types of sequences. It is able to read and write sequences from/into different formats.
 */
template<typename SequenceType>
class SequenceSet
{
private:
	//typedef std::shared_ptr<SequenceType> SequencePtr;
	std::vector<SequenceType> sequences_;
	char seqType_;
	size_t id_;

	/**
	 * \brief Identifies the Alignment format.
	 * @param inF The opened stream pointing to the first line.
	 * @return string denoting the file format.
	 */
	std::string _identifyFileFormat(AlgorithmPack::Input &inF);

	/**
	 * \Functions for reading the different formats.
	 * @param inF The input string
	 * @param seqNames The sequences to include. If empty all sequences will be included.
	 * @param remove instead of beeing kept, the sequences in seqNames are removed.
	 */
	void _readFasta(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove);
	void _readGenbank(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove);
	void _readSwissprot(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove);
	void _readStockholm(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove);
	void _readMsf(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove);
	void _readPhylip(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove);

	// TODO:
	void _readAmps(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove);
	void _readClustalw(AlgorithmPack::Input &inF, const std::map<std::string, short> &seqNames, bool remove);
	void _readCodata(AlgorithmPack::Input &inF, const std::map<std::string, short> &seqNames, bool remove);


	/**
	 * \brief Functions for writing the different formats.
	 * @param outF The output file
	 */
	void _write(std::ostream &outS, const std::string &format, size_t linewidth=80) const;
	void _writeFasta(std::ostream &outF, size_t linewidth) const;

	// TODO:
	void _writeStockholm(std::ofstream &outF) const;
	void _writeGenbank(std::ofstream &outF) const;
	void _writePhylip(std::ofstream &outF) const;
	using SequenceSet_t = std::vector<SequenceType>;

public:
	using iterator = typename SequenceSet_t::iterator;
	using const_iterator = typename SequenceSet_t::const_iterator;
	using reverse_iterator = typename SequenceSet_t::reverse_iterator;
	using const_reverse_iterator = typename SequenceSet_t::const_reverse_iterator;


	/**
	 * \brief The SequenceType saved in the sequence set.
	 */
	typedef SequenceType value_type;

	/**
	 * \brief Standard constructor
	 */
	SequenceSet();

	/**
	 * \brief Constructor which reads a sequence set.
	 * @param input_f The input file.
	 */
	/**
	 * \brief Constructor which reads a sequence set.
	 * @param input_f The input file.
	 * @param seqNames The sequences to extract. If empty all sequences will be extracted.
	 * @param remove instead of beeing kept, the sequences in seqNames are removed.
	 * @param format The format of the input file. If auto the format will be automatically detected.
	 */
	SequenceSet(const std::string &input_f, const std::vector<std::string> &seqNames = std::vector<std::string>(), bool remove = false, const std::string &format = "auto");

	// No copy constructor
	/// @cond HIDDEN
	SequenceSet(const SequenceSet&) = delete;
	/// @endcond

	/**
	 * \brief Move constructor
	 * @param rvalue of a DomamainArrangment
	 */
	SequenceSet(SequenceSet &&) = default;

	/**
	 * \brief Destructor
	 */
	virtual ~SequenceSet();


	/**
	 * \brief Reference to last sequence in the set.
	 * \return Reference to last sequence in the set.
	 */
	SequenceType
	&back()
	{
		return sequences_.back();
	}


	/**
	 * \brief Access operator.
	 * \param index The index of to access.
	 */
	SequenceType
	&operator[](unsigned int index)
	{
		return (sequences_[index]);
	}

	/**
	 * \overload
	 */
	const SequenceType
	&operator[](unsigned int index) const
	{
		return (sequences_[index]);
	}

	/**
	 * \brief Access operator.
	 * \details linear in time.
	 * \param The name of the sequence.
	 */
	SequenceType
	&operator[](const std::string &name)
	{
		return const_cast<SequenceType&>((*static_cast<const SequenceSet<SequenceType> *>(this))[name]);
	}

	/**
	 * \overload
	 */
	const
	SequenceType &operator[](const std::string &name) const
	{
		size_t nSeqs = sequences_.size();
		size_t index;
		for (index=0; index<nSeqs; ++index)
		{
			if (name == sequences_[index].name())
				break;
		}
		if (index != nSeqs)
			return sequences_[index];
		else
			throw std::runtime_error(name + " not found");
	}



	/**
	 * \brief Standard move assignment constructor.
	 * @param The object to move.
	 * @return Reference to the object.
	 */
	SequenceSet<SequenceType > & operator= ( SequenceSet<SequenceType > && ) = default;
	SequenceSet<SequenceType > & operator= ( const SequenceSet<SequenceType > & ) = delete;

	/**
	 * \brief Returns the number of sequences in the set.
	 * @return The number of sequences.
	 */
	size_t
	size() const noexcept
	{
		return sequences_.size();
	}


	/**
	 * \brief Appends a sequence to the set.
	 * @param seq The sequence to add.
	 * \tparam SequenceType.
	 * \details Supports rvalue.
	 */
	template<typename T>
	void
	push_back(T &&seq)
	{
		sequences_.push_back(std::forward<T>(seq));
	}

	/**
	 * \brief Appends a sequence to the set.
	 * @param seq The sequence to add.
	 */
	template<typename ... Args>
	void
	emplace_back(Args&&... args)
	{
		sequences_.emplace_back(std::forward<Args>(args)...);
	}


	/**
	 * \brief Transfers all sequences from another sequence set to the current one.
	 * @param other The other sequence set
	 */
	template<typename SetType>
	void
	transfer(SetType &other)
	{
		size_t n_seqs=other.size();
		for (size_t i=0; i<n_seqs; ++i)
			sequences_.push_back(std::move(other.sequences_[i]));
		other.clear();
	}

	/**
	 * \brief Removes all sequences from the set.
	 */
	void
	clear()
	{
		sequences_.clear();
	}

	/**
	 * Checks if the set is empty.
	 * @return true if empty else false.
	 */
	bool
	empty() const
	{
		return sequences_.empty();
	}


	/**
	 * \brief Reads sequences from a file.
	 * @param input_f The file to read from.
	 * @param seqNames Sequences to extract.
	 * @param format Format of the file.
	 * \details Reads a sequences or a set of sequences from a file. The sequences are added to the existing set.
	 * If seqNames is empty all sequences will be read else only the ones defined in it. With format the sequence format
	 * can be given. If set to auto the format will be automatically described.
	 *
	 */
	virtual void
	read(const fs::path &input_f, const std::vector<std::string> &seqNames = std::vector<std::string>(), bool remove = false, const std::string &format = "auto");

	/**
	 * \brief Writes the sequence set to a file.
	 * @param output_f The output file.
	 * @param format The format to use.
	 */
	virtual void
	write(const fs::path &output_f, const std::string &format, size_t linewidth=80) const;


	virtual void
	write(const std::string &format, size_t linewidth=80) const
	{
		_write(std::cout, format, linewidth);
	}

	/**
	 * \brief find a sequence by name.
	 * @param name The sequence name to look for.
	 * @return Iterator to the found element else an iterator behind the last element.
	 */
	SequenceSet<SequenceType>::iterator
	find(const std::string &name)
	{
		//return const_cast<SequenceSet<SequenceType>::iterator>(static_cast<const SequenceSet<SequenceType> *>(this)->find(name));
		auto endIt =  this->end();
		for (auto it = this->begin(); it != endIt; ++it)
		{
			if (it->name() == name)
				return it;
		}
		return endIt;
	}

	SequenceSet<SequenceType>::const_iterator
	find(const std::string &name) const
	{
		auto endIt =  this->end();
		for (auto it = this->begin(); it != endIt; ++it)
		{
			if (it->name() == name)
				return it;
		}
		return endIt;
	}

	/**
	 * \brief Removes an element from the set.
	 * @param it The element to remove
	 */
	void
	erase(SequenceSet<SequenceType>::iterator it)
	{
		sequences_.erase(it);
	}

	/**
	 * \name Iterator functions
	 */
	/**@{*/

	/**
	 * \brief Returns an iterator to the first element in the container.
	 * @return An iterator to the first element in the container.
	 */
	iterator begin() throw()
	{
		return iterator(sequences_.begin());
	}

	/**
	 * \brief Returns an iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	iterator end() throw()
	{
		return iterator(sequences_.end());
	}

	/**
	 * \brief Returns a const_iterator to the first element in the container.
	 * @return An iterator to the first element in the container.
	 */
	const_iterator begin() const throw()
	{
		return const_iterator(sequences_.begin());
	}


	/**
	 * \brief Returns a const_iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	const_iterator end() const throw()
	{
		return const_iterator(sequences_.end());
	}

	/**
	 * \brief Returns a reverse_iterator to the last element in the container.
	 * @return An iterator to the last element.
	 */
	reverse_iterator rbegin() throw()
	{
		return reverse_iterator(sequences_.end());
	}

	/**
	 * \brief Returns a const_reverse_iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	reverse_iterator rend() throw()
	{
		return reverse_iterator(sequences_.begin());
	}

	/**
	 * \brief Returns a const_iterator to the first element in the container.
	 * @return An iterator to the first element in the container.
	 */
	const_iterator cbegin() throw()
	{
		return const_iterator(sequences_.cbegin());
	}

	/**
	 * \brief Returns a const_iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	const_iterator cend() const throw()
	{
		return const_iterator(sequences_.cend());
	}


	/**
	 * \brief Returns a const_iterator to the first element in the container.
	 * @return An iterator to the first element in the container.
	 */
	const_reverse_iterator crbegin() throw()
	{
		return const_reverse_iterator(sequences_.cend());
	}

	/**
	 * \brief Returns a const_iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	const_reverse_iterator crend() const throw()
	{
		return constreverse_iterator(sequences_.begin());
	}
	/**@}*/
};


using BasicSeqSet = SequenceSet<Sequence<> >;

/***********************************************************
 *                 PUBLIC FUNCTIONS                        *
 ***********************************************************/

//              Constructor & Destructor

template<typename SequenceType>
SequenceSet<SequenceType>::SequenceSet():sequences_(), seqType_('-'), id_(0)
{}

template<typename SequenceType>
SequenceSet<SequenceType>::SequenceSet(const std::string &input_f, const std::vector<std::string> &seqNames, bool remove, const std::string &format) : sequences_(), seqType_('-'), id_(0)
{
	read(input_f, seqNames, remove, format);
}

template<typename SequenceType>
SequenceSet<SequenceType>::~SequenceSet()
{}


//                  IO Functions

template<typename SequenceType>
void
SequenceSet<SequenceType>::read(const fs::path &inputF, const std::vector<std::string> &seqNames, bool remove, const std::string &inFormat)
{
	AlgorithmPack::Input inF(inputF);
	//openInFile(inputF, inF);
	//inF.exceptions ( AlgorithmPack::Input::failbit | AlgorithmPack::Input::badbit );
	//inF.open(input_f);
	//inF.exceptions ( AlgorithmPack::Input::badbit );
	std::map<std::string, short> extractNames;
	for (size_t i=0; i<seqNames.size(); ++i)
		extractNames[seqNames[i]]=0;
	std::string format = inFormat;
	if (format == "auto")
		format = _identifyFileFormat(inF);

	inF.seekg(0);
	if (format=="fasta")
		_readFasta(inF, extractNames, remove);
	else if (format=="genbank")
		_readGenbank(inF, extractNames, remove);
	else if (format=="swissprot")
		_readSwissprot(inF, extractNames, remove);
	else if (format=="stockholm")
		_readStockholm(inF, extractNames, remove);
	else if (format=="msf")
		_readMsf(inF, extractNames, remove);
	else if (format=="phylip")
		_readPhylip(inF, extractNames, remove);
	else if (format=="amps")
		_readAmps(inF, extractNames, remove);

	if (!remove)
	{
		for (auto pair : extractNames)
		{
			if (pair.second == 0)
				throw std::runtime_error("Sequence '" + pair.first + "' not found!");
		}
	}

	inF.close();
}


template<typename SequenceType>
void
SequenceSet<SequenceType>::_write(std::ostream &outF, const std::string &format, size_t linewidth) const
{
	if (format=="fasta")
		_writeFasta(outF, linewidth);
}

template<typename SequenceType>
void
SequenceSet<SequenceType>::write(const fs::path &output_f, const std::string &format, size_t linewidth) const
{
	std::ofstream outF;
	//outF.exceptions ( std::istream_iterator::failbit | std::istream_iterator::badbit );
	outF.open(output_f.string());
	//outF.exceptions ( std::istream_iterator::badbit );
	_write(outF, format, linewidth);
	outF.close();
}


/***********************************************************
 *                PRIVATE FUNCTIONS                        *
 ***********************************************************/


//              Read Functions

template<typename SequenceType>
std::string
SequenceSet<SequenceType>::_identifyFileFormat(AlgorithmPack::Input &inF)
{
	std::string line;
	while (getline(inF, line))
	{
		if (line[0]=='>')
		{
			getline(inF, line);
			if (line[0] != '>')
				return "fasta";
			else
				return "amps";
		}
		if (line.substr(0,5)=="LOCUS")
			return "genbank";
		if ((line[0]=='I') && (line[1]=='D'))
			return "swissprot";
		if (line == "# STOCKHOLM 1.0")
			return "stockholm";
		if (line.find("MSF")!= std::string::npos)
			return "msf";
		if (line[0] == ' ')
			return "phylip";
	}
	return "failed";
}


template<typename SequenceType>
void
SequenceSet<SequenceType>::_readFasta(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove)
{
	std::string line;
	SequenceType *seq = nullptr;
	std::string::size_type pos = 0;
	std::map<std::string, short>::const_iterator it_end=seqNames.end();
	std::map<std::string, short>::iterator it;
	bool extract = false;
	std::string name, comment;
	size_t seqId = sequences_.size()-1;
	while (getline(inF, line))
	{
		if (line.empty())
			continue;
		if (line[0] == '>')
		{
			pos = line.find(' ');
			if (pos!=std::string::npos) //sequence has a comment
			{
				name = std::string(line, 1, pos-1);
				comment = std::string(line, pos+1, line.size()-(pos+1));
			}
			else // sequence has no comment
			{
				name = std::string(line, 1);
				comment="";
			}
			it=seqNames.end();
			if (seqNames.empty() || ((!remove) && ((it=seqNames.find(name)) != it_end)) ||  ((remove) && (seqNames.find(name) == it_end)))
			{
				if (it != it_end)
					++it->second;
				//seq = new SequenceType(std::move(name), "", "", std::move(comment), ++seqId);
				sequences_.emplace_back(std::move(name), "", "", std::move(comment), ++seqId);
				seq = &sequences_.back();
				//sequences_.emplace_back(SequencePtr(seq));
				extract =true;
			}
			else
				extract=false;
		}
		else
		{
			if (extract)
				seq->append(line);
		}
	}
}





template<typename SequenceType>
void
SequenceSet<SequenceType>::_readGenbank(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove)
{
	std::string line;
	SequenceType *seq = nullptr;
	std::map<std::string, short>::const_iterator it_end=seqNames.end();
	std::map<std::string, short>::iterator it;
	bool extract = false;
	size_t seqId = sequences_.size()-1;
	size_t start =  sequences_.size();
	size_t seq_start;
	std::string locus;
	while (getline(inF, line))
	{
		if (line.empty())
			continue;
		if ((line.substr(0,5))=="LOCUS")
		{
			seq_start = 5;
			while (line[++seq_start] == ' ')
				;
			locus = line.substr(seq_start);
			it=seqNames.end();
			if (seqNames.empty() || ((!remove) && ((it=seqNames.find(locus)) != it_end)) ||  ((remove) && (seqNames.find(locus) == it_end)))
			{
				if (it != it_end)
					++it->second;
				//seq = new SequenceType(locus, "", "", "", ++seqId);
				//sequences_.push_back(SequencePtr(seq));
				sequences_.emplace_back(locus, "", "", "", ++seqId);
				seq = &sequences_.back();
				extract =true;
			}
			else
				extract=false;
		}
		else
		{
			if (!extract)
				continue;
			if ((line.substr(0,9))=="ACCESSION")
			{
				seq_start = 9;
				while (line[++seq_start] == ' ');
				seq->accession(line.substr(seq_start));
			}
			if ((line.substr(0,10))=="DEFINITION")
			{
				seq_start = 10;
				while (line[++seq_start] == ' ');
				seq->comment(line.substr(seq_start));
			}
			if (line.substr(0,6)=="ORIGIN")
			{
				while (getline(inF, line))
				{
					if (line[0] == '/')
						break;
					seq_start = 0;
					while (!isalpha(line[++seq_start]))
						;
					seq->append(line.substr(seq_start));
				}
			}


		}
	}
	size_t n_seqs = sequences_.size();
	for (size_t i=start; i<n_seqs; ++i)
		sequences_[i].remove(' ');
}


template<typename SequenceType>
void
SequenceSet<SequenceType>::_readSwissprot(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove)
{
	std::string line;
	SequenceType *seq = nullptr;
	std::map<std::string, short>::const_iterator it_end=seqNames.end();
	std::map<std::string, short>::iterator it;
	bool extract = false;
	size_t seqId = sequences_.size()-1;
	size_t start =  sequences_.size();
	size_t seq_start;
	size_t name_end;
	std::string locus;
	while (getline(inF, line))
	{
		if (line.empty())
			continue;
		if ((line.substr(0,2))=="ID")
		{
			seq_start = 3;
			while (line[++seq_start] == ' ');
			name_end=seq_start;
			while (line[++name_end] != ' ');
			locus = line.substr(seq_start, name_end-seq_start);
			it=seqNames.end();
			if (seqNames.empty() || ((!remove) && ((it=seqNames.find(locus)) != it_end)) ||  ((remove) && (seqNames.find(locus) == it_end)))
			{
				if (it != it_end)
					++it->second;
				//seq = new SequenceType(locus, "", "", "", ++seqId);
				//sequences_.push_back(SequencePtr(seq));
				sequences_.emplace_back(locus, "", "", "", ++seqId);
				seq = &sequences_.back();
				extract =true;
			}
			else
				extract=false;

		}
		else
		{
			if (!extract)
				continue;
			if ((line.substr(0,2))=="AC")
			{
				seq_start = 3;
				while (line[++seq_start] == ' ');
				seq->accession(line.substr(seq_start));
			}
			if ((line.substr(0,2))=="DE")
			{
				seq_start = 3;
				while (line[++seq_start] == ' ');
				seq->comment(line.substr(seq_start));
			}
			if (line.substr(0,2)=="SQ")
			{
				while (getline(inF, line))
				{
					if (line[0] == '/')
						break;
					seq->append(line);
				}
			}


		}
	}
	size_t n_seqs = sequences_.size();
	for (size_t i=start; i<n_seqs; ++i)
		sequences_[i].remove(' ');
}


template<typename SequenceType>
void
SequenceSet<SequenceType>::_readStockholm(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove)
{
	size_t seqId = sequences_.size()-1;
	const std::map<std::string, short>::const_iterator it_end=seqNames.end();
	std::map<std::string, short>::iterator it;
	size_t firstId = seqId+1;
	std::string line;
	bool append = false;
	std::vector<size_t> read_seq;
	//SequenceType *single_seq=NULL;
	size_t name_end, seq_start;
	size_t id=0;
	std::string name;
	while (getline(inF, line))
	{
		if (line.empty())
		{
			append = true;
			id = 0;
			continue;
		}
		if (line[0] == '#')
			continue;
		name_end =0;
		while (line[++name_end] != ' ');
		seq_start = name_end;
		while (line[++seq_start] == ' ');

		if (append)
		{
			if (read_seq[id]!=0)
				sequences_[read_seq[id++]-1].append(line.substr(seq_start));
		}
		else
		{
			name=line.substr(0,name_end);
			it=seqNames.end();
			if (seqNames.empty() || ((!remove) && ((it=seqNames.find(name)) != it_end)) ||  ((remove) && (seqNames.find(name) == it_end)))
			{
				if (it != it_end)
					++it->second;
				sequences_.emplace_back(std::move(name), line.substr(seq_start), "", "", ++seqId);
				read_seq.push_back(sequences_.size());
			}
			else
				read_seq.push_back(0);
		}
	}

	size_t n_seqs = sequences_.size();
	for (size_t i=firstId; i<n_seqs; ++i)
		replace(sequences_[i], '.', '-');

}

template<typename SequenceType>
void
SequenceSet<SequenceType>::_readMsf(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove)
{
	std::string line;
	size_t pos, pos2, aln_len=0;
	std::string name;
	std::vector<short> use;
	size_t seqId = sequences_.size()-1;
	size_t firstId = sequences_.size();
	const std::map<std::string, short>::const_iterator it_end=seqNames.end();
	std::map<std::string, short>::iterator it;

	// read header
	while (getline(inF, line))
	{
		if (line[0]=='/')
			break;
		if ((pos=line.find("MSF")) != std::string::npos)
		{
			pos+=3;
			while (isblank(line[++pos]));
			pos2=pos;
			while (!isblank(line[++pos2]));
			aln_len =std::stoul(line.substr(pos, pos2-pos));
		}
		else if (((pos=line.find("Name")) != std::string::npos) || ((pos=line.find("NAME")) != std::string::npos))
		{
			pos+=4;
			while (isblank(line[++pos]));
			pos2=pos;
			while (!isblank(line[++pos2]));
			name=line.substr(pos, pos2-pos);
			it=seqNames.end();
			if (seqNames.empty() || ((!remove) && ((it=seqNames.find(name)) != it_end)) ||  ((remove) && (seqNames.find(name) == it_end)))
			{
				if (it != it_end)
					++it->second;
				sequences_.emplace_back(name, aln_len, "", "",  ++seqId);
				use.push_back(seqId);
			}
			else
				use.push_back(-1);
		}

	}
	// read sequences
	short seq_number = -1;
	while (getline(inF, line))
	{
		if (line.empty())
			seq_number=-1;
		else
		{
			if (use[++seq_number]!=-1)
			{
				pos=0;
				while (!isblank(line[++pos]))
					;
				while (isblank(line[++pos]))
					;
				sequences_[use[seq_number]].append(line.substr(pos));
			}
		}
	}
	size_t n_seqs = sequences_.size();
	for (size_t i=firstId; i<n_seqs; ++i)
		sequences_[i].remove(' ');
	for (size_t i=firstId; i<n_seqs; ++i)
		replace(sequences_[i], '.', '-');
}


template<typename SequenceType>
void
SequenceSet<SequenceType>::_readPhylip(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove)
{
	std::string line;
	std::vector<std::string> tokens;
	getline(inF, line);
	split(line, "\t ", tokens);
	size_t nSeqs = std::stoul(tokens[0]);
	size_t seqLength = std::stoul(tokens[1]);
    int filePos = inF.tellg();

    // determine phylip format:
    // if spaces exist at first position => interleaved format
    bool isSequential = true;
    for (size_t i=0; i<nSeqs+1; ++i)
    	getline(inF, line);
    if (line[0] == ' ')
    	isSequential = false;
    else
    {
    	getline(inF, line);
    	if (line[0] == ' ')
    	    isSequential = false;
    }

    // read sequences
	const std::map<std::string, short>::const_iterator it_end=seqNames.end();
	std::map<std::string, short>::iterator it;
    inF.seekg (filePos, std::istream::beg);
	size_t seqId = sequences_.size()-1;
	size_t firstId = sequences_.size();
	std::string name="", seq;
	size_t currentLength=seqLength;
	SequenceType *currentSeq = nullptr;
    if (isSequential)
    {
		while (getline(inF, line))
		{
			if (currentLength == seqLength)
			{
				name = line.substr(0,10);
				trimRight(name);
				seq = line.substr(10);
				removeSpaces(seq);
				currentLength = seq.size();
				it=seqNames.end();
				if (seqNames.empty() || ((!remove) && ((it=seqNames.find(name)) != it_end)) ||  ((remove) && (seqNames.find(name) == it_end)))
				{
					if (it != it_end)
						++it->second;
					sequences_.emplace_back(name, seqLength, "", "",  ++seqId);
					currentSeq = &sequences_[sequences_.size()-1];
					currentSeq->append(seq);
				}
				else
					currentSeq = nullptr;
			}
			else
			{
				removeSpaces(line);
				currentLength += line.size();
				if (currentSeq)
					currentSeq->append(line);
			}
		}
    }
    else // not sequenctial
    {
    	std::vector<int> use(nSeqs,0);
    	// read first lines with names
    	for (size_t i=0; i<nSeqs; ++i)
    	{
    		getline(inF, line);
    		name = line.substr(0,10);
			trimRight(name);
			seq = line.substr(10);
			removeSpaces(seq);
			it=seqNames.end();
			if (seqNames.empty() || ((!remove) && ((it=seqNames.find(name)) != it_end)) ||  ((remove) && (seqNames.find(name) == it_end)))
			{
				if (it != it_end)
					++it->second;
				use[i] =1;
				sequences_.emplace_back(name, seqLength, "", "",  ++seqId);
				currentSeq = &sequences_[sequences_.size()-1];
				currentSeq->append(seq);
			}
			else
				currentSeq = NULL;
    	}

    	// read rest of file
    	size_t seqNo=0;
    	size_t id = 0;
		while (getline(inF, line))
		{
			removeSpaces(line);
			if (line.empty())
				continue;
			if (use[seqNo])
			{
				sequences_[firstId+id].append(line);
				++id;
			}
			if (++seqNo == nSeqs)
				seqNo = id = 0;
		}
    }

}

template<typename SequenceType>
void
SequenceSet<SequenceType>::_readAmps(AlgorithmPack::Input &inF, std::map<std::string, short> &seqNames, bool remove)
{
	const std::map<std::string, short>::const_iterator it_end=seqNames.end();
	std::map<std::string, short>::iterator it;
	std::string line;
	size_t seqId = sequences_.size()-1;
	size_t firstId = sequences_.size();
	std::string name;
	std::vector<int> use;
	while (getline(inF, line))
	{
		if (line[0] == '>')
		{
			name = line.substr(1);
			it=seqNames.end();
			if (seqNames.empty() || ((!remove) && ((it=seqNames.find(name)) != it_end)) ||  ((remove) && (seqNames.find(name) == it_end)))
			{
				if (it != it_end)
					++it->second;
				//sequences_.push_back(SequencePtr(new SequenceType(name, static_cast<size_t>(0), "", "",  ++seqId)));
				sequences_.emplace_back(name, static_cast<size_t>(0), "", "",  ++seqId);
				use.push_back(1);
			}
			else
				use.push_back(0);
		}
		if (line[0] == '*')
			break;
	}

	while (getline(inF, line))
	{
		if (line[0] == '*')
			break;
		size_t thisSeq = 0;
		for (size_t i=0; i<line.length(); ++i)
		{
			if (use[i])
				sequences_[firstId+(thisSeq++)].append(line[i]);
		}

	}
	size_t nSeqs = sequences_.size();
	for (size_t i=firstId; i<nSeqs; ++i)
		replace(sequences_[i], ' ', '-');
}

//              Write functions

template<typename SequenceType>
void
SequenceSet<SequenceType>::_writeFasta(std::ostream &outF, size_t linewidth) const
{
	size_t n_seqs=sequences_.size();
	for (size_t i=0; i<n_seqs; ++i)
	{
		size_t seqLength = sequences_[i].size();
		const SequenceType &seq = sequences_[i];
		if (sequences_[i].comment().empty())
			outF << ">" << sequences_[i].name()  << "\n";
		else
		 	outF << ">" << sequences_[i].name() << " " << sequences_[i].comment() << "\n";
		size_t pos = 0;
		while (pos < seqLength)
		{
			for (size_t endPos = std::min(pos+linewidth, seqLength); pos<endPos; ++pos )
				outF << seq[pos];
			outF << "\n";
		}
	}
}




} // BioSeqDataLib

#endif /* SEQUENCESET_HPP_ */
