/*
 * SeqFunctions.hpp
 *
 *  Created on: 19 Oct 2013
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
 * \file SeqFunctions.hpp
 * \brief Functions related to the Sequence class.
 *
 * \details This file contains several sequences to manipulate and analyse sequences.
 *
 */


#include<cctype>
#include<algorithm>
#include<string>
#include<set>
#include<map>
#include<set>

#include "../utility/properties.hpp"
#include "../utility/TwoValues.hpp"

#ifndef SEQFUNCTIONS_HPP_
#define SEQFUNCTIONS_HPP_

namespace BioSeqDataLib
{

/** @defgroup seqModifyFunction Sequence modifying functions
 *  Group with basic sequence modifying functions
 *  @{
 */



/**
 * \brief Replaces one character with another.
 * @param seq The sequence to modify.
 * @param c1 The character to replace.
 * @param c2 The new character.
 * \relates Sequence
 */
template<typename SequenceType>
void
replace(SequenceType &seq, char c1, char c2)
{
	size_t size=seq.size();
	for (size_t i=0; i<size; ++i)
	{
		if (seq[i]== c1)
			seq[i] = c2;
	}
}

/**
 * \brief Calculates the reverse complement of a DNA sequence.
 * @param seq The sequence to transform.
 *
 * \details The upper and lowercase structure of the sequence is maintained.
 * \relates Sequence
 */
template<typename SequenceType>
void
reverseComplement(SequenceType &seq)
{
	size_t length = seq.size();
	std::reverse(seq.begin(), seq.end());
	std::vector<char>transformation(128,'!');
	transformation['A']='T';
	transformation['a']='t';
	transformation['C']='G';
	transformation['c']='g';
	transformation['T']='A';
	transformation['t']='a';
	transformation['G']='C';
	transformation['g']='c';
	transformation['R']='Y';
	transformation['r']='y';
	transformation['Y']='R';
	transformation['y']='r';
	transformation['M']='K';
	transformation['m']='k';
	transformation['K']='M';
	transformation['k']='m';
	transformation['W']='W';
	transformation['w']='w';
	transformation['S']='S';
	transformation['s']='s';
	transformation['B']='V';
	transformation['b']='v';
	transformation['V']='B';
	transformation['v']='b';
	transformation['D']='H';
	transformation['d']='h';
	transformation['H']='D';
	transformation['h']='d';
	transformation['N']='N';
	transformation['n']='n';
	transformation['-']='-';

	for (size_t i=0; i<length; ++i)
		seq[i] = transformation[seq[i]];
}


/**
 * \brief Translation of a nucleotide sequence (DNA/RNA) into a protein sequence.
 * @param seq The sequence to translate.
 * @return The protein sequence.
 * \details Currently only the standard code is supported.
 * \relates Sequence
 */
template<typename SequenceType>
void
translate(SequenceType &seq)
{
	static const std::map<std::string, char> codon2aa
	{
		{"ATT",'I'}, {"ATC",'I'}, {"ATA",'I'}, {"CTT",'L'}, {"CTC",'L'}, {"CTA",'L'}, {"CTG",'L'}, {"TTA",'L'}, {"TTG",'L'}, {"GTT",'V'},
		{"GTC",'V'}, {"GTA",'V'}, {"GTG",'V'}, {"TTT",'F'}, {"TTC",'F'}, {"ATG",'M'}, {"TGT",'C'}, {"TGC",'C'}, {"GCT",'A'}, {"GCC",'A'},
		{"GCA",'A'}, {"GCG",'A'}, {"GGT",'G'}, {"GGC",'G'}, {"GGA",'G'}, {"GGG",'G'}, {"CCT",'P'}, {"CCC",'P'}, {"CCA",'P'}, {"CCG",'P'},
		{"ACT",'T'}, {"ACC",'T'}, {"ACA",'T'}, {"ACG",'T'}, {"TCT",'S'}, {"TCC",'S'}, {"TCA",'S'}, {"TCG",'S'}, {"AGT",'S'}, {"AGC",'S'},
		{"TAT",'Y'}, {"TAC",'Y'}, {"TGG",'W'}, {"CAA",'Q'}, {"CAG",'Q'}, {"AAT",'N'}, {"AAC",'N'}, {"CAT",'H'}, {"CAC",'H'}, {"GAA",'E'},
		{"GAG",'E'}, {"GAT",'D'}, {"GAC",'D'}, {"AAA",'K'}, {"AAG",'K'}, {"CGT",'R'}, {"CGC",'R'}, {"CGA",'R'}, {"CGG",'R'}, {"AGA",'R'},
		{"AGG",'R'}, {"TAA",'*'}, {"TAG",'*'}, {"TGA",'*'},
		{"AUU",'I'}, {"AUC",'I'}, {"AUA",'I'}, {"CUU",'L'}, {"CUC",'L'}, {"CUA",'L'}, {"CUG",'L'}, {"UUA",'L'}, {"UUG",'L'}, {"GUU",'V'},
		{"GUC",'V'}, {"GUA",'V'}, {"GUG",'V'}, {"UUU",'F'}, {"UUC",'F'}, {"AUG",'M'}, {"UGU",'C'}, {"UGC",'C'}, {"GCU",'A'}, {"GCC",'A'},
		{"GGU",'G'}, {"CCU",'P'}, {"ACU",'T'}, {"UCU",'S'}, {"UCC",'S'}, {"UCA",'S'}, {"UCG",'S'}, {"AGU",'S'},	{"UAU",'Y'}, {"UAC",'Y'},
		{"UGG",'W'}, {"AAU",'N'}, {"CAU",'H'}, {"GAU",'D'}, {"CGU",'R'}, {"UAA",'*'}, {"UAG",'*'}, {"UGA",'*'}
	};

	static const std::map<std::string, char> codonS2aa
	{
		{"GT",'V'}, {"GC",'A'}, {"GG",'G'}, {"TC",'S'}, {"CT",'L'}, {"CC",'P'}, {"CG",'R'}, {"AC",'T'}
	};

	std::map<std::string, char>::const_iterator it, itEnd=codon2aa.end(), itEndS=codonS2aa.end();
	size_t seqLength = seq.size();
	std::string codon="TTT";
	std::string codonS="TT";
	size_t pos = 0;
	for (size_t i=0; i<seqLength; i+=3)
	{
		codonS[0] = codon[0] = std::toupper(seq[i]);
		codonS[1] = codon[1] = std::toupper(seq[i+1]);
		codon[2]=std::toupper(seq[i+2]);

		if ((it=codonS2aa.find(codonS)) != itEndS)
			seq[pos++] = it->second;
		else
		{
			if ((it=codon2aa.find(codon)) != itEnd)
				seq[pos++] = it->second;
			else
				seq[pos++] = 'X';
		}
	}
	seq.resize(pos);
}

/**
 * \brief Changes all letters to upper case.
 * @param seq The sequence.
 * \relates Sequence
 */
template<typename SequenceType>
void
toUpper(SequenceType &seq)
{
	size_t length = seq.size();
	for (size_t i=0; i<length; ++i)
		seq[i]=std::toupper(seq[i]);
}

/**
 * \brief Changes all letters to lower case.
 * @param seq
 * \relates Sequence
 */
template<typename SequenceType>
void
toLower(SequenceType &seq)
{
	size_t length = seq.size();
	for (size_t i=0; i<length; ++i)
		seq[i]=std::tolower(seq[i]);
}

/**
 * \brief Produces a Sequence, that is a subsequence of an existing one.
 * @param seq The original sequence.
 * @param first The first position of the subsequence.
 * @param last The last position of the subsequence.
 * @return The new sequence.
 */
template<typename SequenceType>
SequenceType
subseq(const SequenceType &seq, size_t first, size_t last)
{
	if (last < seq.size())
		return SequenceType(seq.name(), seq.seq().substr(first, last-first+1), seq.accession(), "");//, std::to_string(first)+" "+std::to_string(last));
	else
		throw std::out_of_range("subseq");
}




template<typename SequenceType>
std::pair<long, long>
orf_in_frame(const SequenceType &seq, const std::set<std::string>
&starts, const std::set<std::string> &stops, short frame, bool incorrectEnds)
{
	long max_start = -1;
    long max_end = -1;

    long frame_start = -1;
    long frame_end = -1;

    long pos = frame;
    long len = seq.size();
    std::string codon = "XXX";
    while (pos+2 < len)
    {
        codon[0] = std::toupper(seq[pos]);
        codon[1] = std::toupper(seq[++pos]);
        codon[2] = std::toupper(seq[++pos]);
        if (frame_start == -1)
        {
            if (starts.find(codon) != starts.end())
                frame_start = pos-2;
        }
        else
        {
            if (stops.find(codon) != stops.end())
			{
                frame_end = pos;
				if ((frame_end-frame_start) > (max_end - max_start))
				{
					max_end = frame_end;
					max_start = frame_start;
				}
				frame_start = -1;
			}
        }
		++pos;
    }
	if ((incorrectEnds) && (frame_start != -1) && (((pos-1)-frame_start) > (max_end - max_start)))
	{
		max_end = pos-1;
		max_start = frame_start;
	}
    return std::pair<long, long>(max_start, max_end);
}


/**
 * \brief Calcualates the longest Orf position.
 * \param seq The sequence in which the longest orft should
 * \param starts The set of start codons
 * \param stops The set of stop codonS
 * \param reverse 0=only + strand, 1=+/- strand, 2=only - strand
 */
template<typename SequenceType>
std::pair<long, long>
longestOrf(const SequenceType &seq, const std::set<std::string> &starts,
const std::set<std::string> &stops, int reverse=0, bool incorrectEnds=true)
{
	long max_length = -1;
    long max_start = -1;
    long max_end = -1;

    if (reverse < 2)
    {
		for (short i=0; i<3; ++i)
		{
			auto p = orf_in_frame(seq, starts, stops, i, incorrectEnds);
			if ((p.second - p.first) > max_length)
			{
				max_start = p.first;
				max_end = p.second;
				max_length = p.second-p.first;
			}
		}
    }

	if (reverse > 0)
	{
		SequenceType newSeq(seq);
		reverseComplement(newSeq);
		for (short i=0; i<3; ++i)
		{
			auto p = orf_in_frame(newSeq, starts, stops, i, incorrectEnds);
			if ((p.second - p.first) > max_length)
			{
				max_length = p.second-p.first;
				max_start = -1*p.first;
				max_end = -1*p.second;

			}
		}
	}

	return std::pair<long,long>(max_start, max_end);
}



/** @} */ //  Sequence modify functions

//*********************************************************************************
//                     Sequence Analysis Functions                                *
//*********************************************************************************

/** @defgroup seqAnalysisFunction Sequence analysis functions
 *  Group with basic sequence analysis functions
 *  @{
 */



/**
 * \brief Calculates the Kyte-Doolittle hydropathy scores for a sequence.
 * @param seq The sequence to analyze.
 * @param windowSize The window size to use (has to be an odd size)
 * @return The calculated scores.
 * \relates Sequence
 */
template<typename SequenceType>
std::vector<float>
calcHydropathy(const SequenceType &seq, int windowSize)
{
	static std::vector<float> aminoAcidScores = hydropathyScores();
	std::vector<float> scores(seq.size()-windowSize+1);
	float sum = 0;
	for (int i=0; i<windowSize; ++i)
		sum += aminoAcidScores[seq[i]];
	scores[0] = sum/windowSize;
	for (size_t i=windowSize; i<seq.size(); ++i)
	{
		sum += aminoAcidScores[seq[i]];
		sum -= aminoAcidScores[seq[i-windowSize]];
		scores[i-windowSize+1] = sum/(windowSize);
	}
	return scores;
}


/**
 * \brief Calculates the CpG o/e value of a DNA sequence.
 *
 * \details The formula used to calculate the CpG o/e value is:  \f$\frac{f_{CG}}{f_C*f_G}\f$, where  \f$ f_x\f$ is the frequency of \f$ x\f$.
 * @param seq The sequence to analyse.
 * @param removeN If set to true the 'N' residues will be 'removed'.
 * @param threshold The minimum percentage of length to use.
 * @return The CpG o/e value.
 */
template<typename SequenceType>
float
calcCpGoe(const SequenceType &seq, bool removeN=true, float threshold=100)
{
	size_t seqLength = seq.size();
	size_t dinucCG = 0, singleC=0, singleG=0, singleCounter=0, dinucCounter=0;
	char residue, prevResidue;
	if (!removeN)
	{
		prevResidue = residue = toupper(seq[0]);
		if (residue == 'C')
			++singleC;
		else
		{
			if (residue == 'G')
				++singleG;
		}
		for (size_t i=1; i<seqLength; ++i)
		{
			residue = toupper(seq[i]);
			switch (residue)
			{
				case 'C':
					++singleC;
				break;
				case 'G':
					++singleG;
					if (prevResidue == 'C')
						++dinucCG;
				break;
			}
			prevResidue = residue;
		}
		singleCounter = seqLength;
		dinucCounter = seqLength-1;
	}
	else
	{
		prevResidue = residue = toupper(seq[0]);
		if (residue == 'C')
			++singleC;
		else
		{
			if (residue == 'G')
				++singleG;
		}
		if (residue != 'N')
			++singleCounter;
		if (toupper(seq[0]) == 'C')
			++singleC;
		else
		{
			if (toupper(seq[0]) == 'G')
				++singleG;
		}
		for (size_t i=1; i<seqLength; ++i)
		{
			residue = toupper(seq[i]);
			switch (residue)
			{
				case 'C':
					++singleC;
				break;
				case 'G':
					++singleG;
					if (prevResidue == 'C')
						++dinucCG;
				break;
			}
			if (residue != 'N')
			{
				++singleCounter;
				if (prevResidue != 'N')
					++dinucCounter;
			}
			prevResidue = residue;
		}
	}
	if (singleCounter < 1.0*seqLength*threshold/100.0)
		return -1;
	if ((singleC > 0) && (singleG > 0))
        return ((1.0*dinucCG/dinucCounter)/((1.0*singleC/singleCounter)*(1.0*singleG/singleCounter)));
    else
        return 0;

}

/**
 * \brief Calculates the CG content of a sequence.
 * @param seq The sequence to analyze.
 * @param removeN If set to true the 'N' residues will be 'removed'.
 * @return The fraction of CG.
 */
template<typename SequenceType>
TwoValues
calcCGcontent(const SequenceType &seq, bool removeN=true)
{
	size_t CGcount = 0, totalCount= 0;
	char residue;
	if (removeN)
	{
		for (size_t i=0; i<seq.size(); ++i)
		{
			residue = toupper(seq[i]);
			if ((residue == 'C') || (residue == 'G'))
				++CGcount;
			if (residue != 'N')
				++ totalCount;
		}
	}
	else
	{
		for (size_t i=0; i<seq.size(); ++i)
		{
			residue = toupper(seq[i]);
			if ((residue == 'C') || (residue == 'G'))
				++CGcount;
		}
		totalCount = seq.size();
	}

	return TwoValues(CGcount, totalCount);
}

/** @} */ //  Sequence analysis functions


} // BioSeqDataLib


#endif /* SEQFUNCTIONS_HPP_ */
