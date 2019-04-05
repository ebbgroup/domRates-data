/*
 * Alignment.hpp
 *
 *  Created on: 22 Oct 2013
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
 * \file Alignment.hpp
 * \brief File containing the alignment class.
 */

#ifndef ALIGNMENT_HPP_
#define ALIGNMENT_HPP_

#include "SequenceSet.hpp"

namespace BioSeqDataLib
{

/**
 * \brief Class to represent an alignment.
 * \tparam SequenceType The sequence type.
 *
 * An alignment is very similar to a SequenceSet but has certain specifics. All sequences are of the same length, and no column consists of gap characters only.
 *
 */
template<typename SequenceType>
class Alignment : public SequenceSet<SequenceType>
{

private:
	/**
	 * \brief Writes the alignment in msf format.
	 * @param outF The output file.
	 */
	void _writeMSF(std::ofstream &outF) const;
	void _writePhylipS(std::ofstream &outF) const;
	/**
	 * \brief Checks if the alignment is a correct alignment.
	 * @return true, if its a correct alignment else false.
	 */
	bool _check() const;

public:
	/**
	 * \brief The SequenceType saved in the sequence set.
	 */
	typedef SequenceType value_type;

	/**
	 * \brief Reads an alignment from file.
	 * @param inputF The alignment file.
	 * @param seq_names The sequences to extract.
	 * @param remove instead of beeing kept, the sequences in seqNames are removed.
	 * @param format The format of the alignment. If "auto" it will be determined automatically.
	 * @param check Checks if the alignment is a proper alignment.
	 */
	Alignment(const std::string &input_f, const std::vector<std::string> &seqNames = std::vector<std::string>(), bool remove=false, const std::string &format = "auto", bool ceck=true);

	/**
	 * \brief Standard constructor.
	 */
	Alignment()
	{}

	/**
	 * \brief Standard destructor.
	 */
	virtual ~Alignment()
	{}

	/**
	 * \brief Reads an alignment from file.
	 * @param inputF The alignment file.
	 * @param seq_names The sequences to extract.
	 * @param remove instead of beeing kept, the sequences in seqNames are removed.
	 * @param format The format of the alignment. If "auto" it will be determined automatically.
	 * @param check Checks if the alignment is a proper alignment.
	 */
	void
	read(const std::string &inputF, const std::vector<std::string> &seq_names = std::vector<std::string>(), bool remove=false, const std::string &format = "auto", bool check=true)
	{
		SequenceSet<SequenceType>::read(inputF, seq_names, remove, format);
		if (check)
			_check();

	}

	/**
	 * \brief Writes an alignment.
	 * @param outputF The output file.
	 * @param format The format of the alignment to write.
	 * @param linewidth The linewidth to use.
	 */
	void
	write(const std::string &outputF, const std::string &format, size_t linewidth=80) const
	{
		if (format == "fasta")
		{
			SequenceSet<SequenceType>::write(outputF, format, linewidth);
		}
		else
		{
			std::ofstream outF;
			outF.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
			outF.open(outputF);
			outF.exceptions ( std::ifstream::badbit );
			if (format == "msf")
			{
				_writeMSF(outF);
			}
			else if (format == "phylips")
			{
				_writePhylipS(outF);
			}
			outF.close();
		}
	}

	/**
	 * \brief Returns the length of the alignment.
	 * @return The length of the alignment.
	 */
	size_t
	length() const
	{
		if (this->empty())
			return 0;
		else
			return (*this)[0].size();
	}


	/**
	 * \brief Transfers the sequences from one alignment to another.
	 * @param other The other alignment.
	 */
	void
	transfer(Alignment<SequenceType> &other)
	{

		if ((!this->empty()) && (length() != other.length()))
			throw std::invalid_argument("Sequences of different lengths");
		SequenceSet<SequenceType>::transfer(other);
	}

	/**
	 * \brief Return all sequences in ungapped format.
	 * @return The ungapped sequences.
	 */
	SequenceSet<SequenceType>
	seqSet()
	{
		SequenceSet<SequenceType> set;
		set.transfer(*this);
		size_t nSeqs = set.size();
		for (size_t i=0; i<nSeqs; ++i)
			set[i].remove('-');
		return set;
	}

};

template<typename SequenceType>
Alignment<SequenceType>::Alignment(const std::string &input_f, const std::vector<std::string> &seqNames, bool remove, const std::string &format, bool check)
{
	read(input_f, seqNames, remove, format, check);
}



template<typename SequenceType>
bool
Alignment<SequenceType>::_check() const
{
	size_t nSeqs = this->size();
	if (nSeqs == 0)
		return true;

	size_t len = (*this)[0].length();
	for (size_t i=0; i<nSeqs; ++i)
	{
		if ((*this)[i].length() != len)
			return false;
	}
	return true;
}

template<typename SequenceType>
void
Alignment<SequenceType>::_writeMSF(std::ofstream &outF) const
{
	outF << "PileUp\n\nMSF: " << (*this)[0].length() << " Type: P Check: 0 ..\n\n";
	size_t nSeqs=this->size();
	std::vector<std::string> names;
	names.reserve(nSeqs);
	size_t seqLength = this->length();
	for (size_t i=0; i<nSeqs; ++i)
	{
		names.push_back((*this)[i].name());
		names[i].resize(10, ' ');
	}
	for (size_t i=0; i<nSeqs; ++i)
		outF << "Name: " << names[i] << " Len: " << (*this)[0].length() << " Check: 0 Weight: 0\n";
	outF << "\n//\n\n";
	char c;
	size_t pos=0;
	size_t linewidth=50;
	while (pos < seqLength)
	{
		size_t endPos = std::min(pos+linewidth, seqLength);
		for (size_t i=0; i<nSeqs; ++i)
		{

			const SequenceType &seq = (*this)[i];
			outF << names[i] << " ";
			for (size_t tmpPos=pos;tmpPos<endPos; ++tmpPos )
			{
				if ((tmpPos % 10 == 0) && ((tmpPos % 50) != 0))
						outF << " ";
				if ((c=seq[tmpPos]) == '-')
					outF << ".";
				else
					outF << c;
			}
			outF << "\n";
		}
		outF << "\n";
		pos+=linewidth;
	}
}



/*
	PileUp

	 MSF:  383  Type: P    Check:  3696   ..

	 Name: ALEU_HORVU oo  Len:  383  Check:  9840  Weight:  34.4
	 Name: CATH_HUMAN oo  Len:  383  Check:  7134  Weight:  30.1
	 Name: CYS1_DICDI oo  Len:  383  Check:  6722  Weight:  35.4

	//


	ALEU_HORVU      MAHARVLLLA LAVLATAAVA VASSSSFADS NPIRPVTDRA ASTLESAVLG
	CATH_HUMAN      ......MWAT LPLLCAGAWL LG........ VPVCGAAELS VNSLEK....
	CYS1_DICDI      .....MKVIL LFVLAVFTVF VS........ .......SRG IPPEEQ....
*/

template<typename SequenceType>
void
Alignment<SequenceType>::_writePhylipS(std::ofstream &outF) const
{
	size_t nSeqs = this->size();
	size_t alnLength = this->length();
	outF << " " << nSeqs << " " << alnLength << "\n";
	size_t linewidth=80;
	for (size_t i=0; i<nSeqs; ++i)
	{
		std::string seq = (*this)[i].seq();
		std::string name = (*this)[i].name();
		name.resize(10, ' ');
		outF << name;
		size_t minVal = std::min(linewidth-10, alnLength);
		outF << seq.substr(0, minVal) << "\n";
		for (size_t j=linewidth-10; j<alnLength; j+=linewidth)
		{
			minVal = std::min(linewidth, alnLength-j);
			outF << seq.substr(j, minVal) << "\n";
		}
	}
}

} //BioSeqDataLib

#endif /* ALIGNMENT_HPP_ */
