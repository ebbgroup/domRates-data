/*
 * DSM.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: ckeme_01
 */

#include "DSM.hpp"

namespace BioSeqDataLib
{


using namespace std;

DSM::DSM() : nDomains_(0), threshold_(0), useNegative_(false)
{
	// TODO Auto-generated constructor stub

}

DSM::~DSM()
{
	// TODO Auto-generated destructor stub
}

void
DSM::read(const fs::path &inFile)
{
	ifstream inF(inFile.string(), ios::in | ios::binary);
	if (inF.bad())
		throw std::runtime_error("Error: A problem occurred opening 'inFile'.");

	// read description
	int descLen;
	inF.read((char*)&descLen, sizeof(int));
	description_.resize(descLen+1);
	inF.read((char*)&description_[0], sizeof(char)*descLen);
	description_[descLen]='\0';

	// read threshold
	inF.read((char*)&threshold_, sizeof(short));

	// read domain information
	inF.read((char*)&nDomains_, sizeof(int));
	int totNameLength;
	inF.read((char*)&totNameLength, sizeof(int));
	string namesS;
	namesS.resize(totNameLength);
	inF.read((char*)&namesS[0], sizeof(char)*totNameLength);
	auto tokens = split(namesS,";");

	for (int i=0; i<nDomains_; ++i)
		acc2id_.emplace(tokens[i], i);

	size_t nVals;
	inF.read((char*)&nVals, sizeof(size_t));
	rowIDs_.resize(nDomains_);
	colIDs_.resize(nVals);
	values_.resize(nVals);
	inF.read((char*)&rowIDs_[0], sizeof(int) * nDomains_);
	inF.read((char*)&colIDs_[0], sizeof(int) * nVals);
	inF.read((char*)&values_[0], sizeof(short) * nVals);
	inF.close();
}


short
DSM::val(const string &d1, const string &d2) const
{
	size_t id1, id2;
	try
	{
		id1 = acc2id_.at(d1);
	}
	catch (std::out_of_range &e)
	{
		throw std::out_of_range("DSM: Error - " + d1 + " not contained" );
	}

	try
	{
		id2 = acc2id_.at(d2);
	}
	catch (std::out_of_range &e)
	{
		throw std::out_of_range("DSM: Error - " + d2 + " not contained" );
	}

	int posR, posRN, posC;

	// only one half of the matrix is stored, therefore need for checking order.
	if (id1 < id2)
	{
		posR = rowIDs_[id1];
		posRN =  (++id1 == rowIDs_.size()) ? colIDs_.size() : rowIDs_[id1];
		posC = id2;
	}
	else
	{
		posR = rowIDs_[id2];
		posRN =  (++id2 == rowIDs_.size()) ? colIDs_.size() : rowIDs_[id2];
		posC = id1;
	}

	// check the colIDs if value is contained.
	for (int i=posR; i<posRN; ++i)
	{
		if (colIDs_[i] == posC)
		{
			if (!useNegative_)
				return values_[i];
			else
				return 2*values_[i]-100;
		}
	}
	if (!useNegative_)
		return 0;
	else
		return -100;
}



} /* namespace BioSeqDataLib */
