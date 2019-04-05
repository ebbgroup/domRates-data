/*
 * OrthologySet.cpp
 *
 *  Created on: 18 Aug 2014
 *      Author: Carsten Kemena
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
#include "OrthologySet.hpp"

namespace BioSeqDataLib
{

using namespace std;
OrthologySet::OrthologySet()
{
	// TODO Auto-generated constructor stub

}

OrthologySet::~OrthologySet()
{
	// TODO Auto-generated destructor stub
}

void
OrthologySet::read(const string &orthoFile)
{
	AP::Input orthoS(orthoFile);
	//openInFile(orthoFile, orthoS);
	string line;
	getline(orthoS, line);
	orthoS.seekg(0);
	if (line[0] == '#')
		readProteinortho_(orthoS);
	else
		readOrthoMCL_(orthoS);
	orthoS.close();
}

void
OrthologySet::readOrthoMCL_(AP::Input &orthoS)
{

	//ifstream orthoS;
	//openInFile(orthoFile, orthoS);
	string line;
	while (getline(orthoS, line))
	{
		vector<string> tokens = split(line, " ");
		tokens[0].resize(tokens[0].size()-1);
		size_t nTokens = tokens.size();
		for (size_t i=1; i<nTokens; ++i)
		{
			string &token = tokens[i];
			size_t pos = token.find("|");
			std::string speciesID = token.substr(0,pos);
			speciesSet_.insert(speciesID);
			groups_[tokens[0]][speciesID].push_back(token.substr(pos+1));
			id2group_[speciesID][token.substr(pos+1)] = tokens[0];
		}
		//2ants1000: pbar|PB12891-RA pbar|PB14305-RA pbar|PB14364-RA pbar|PB17240-RA sinv|SINV15920-PA sinv|SINV10134-PA sinv|SINV26408-PA sinv|SINV23932-PA sinv|SINV26415-PA sinv|SINV10298-PA sinv|SINV15232-PA sinv|SINV12123-PA sinv|SINV22168-PA sinv|SINV22696-PA sinv|SINV11930-PA sinv|SI
	}
}

void
OrthologySet::readProteinortho_(AP::Input &orthoS)
{
	string line;
	getline(orthoS, line);
	auto header = split(line, "\t");
	size_t lineN=0;
	for (size_t i=3; i<header.size(); ++i)
		speciesSet_.insert(header[i]);
	while (getline(orthoS, line))
	{
		++lineN;
		vector<string> tokens = split(line, "\t");
		size_t nTokens = tokens.size();
		for (size_t i=3; i<nTokens; ++i)
		{
			if (tokens[i] == "*")
				continue;
			auto seqNames=split(tokens[i], ",");
			for (auto &seqName : seqNames)
			{
				groups_[to_string(lineN)][header[i]].push_back(seqName);
				id2group_[header[i]][seqName] = to_string(lineN);
			}
		}
	}
}


const OrthoGroup *
OrthologySet::find(const std::string &species, const std::string &seqID) const
{
	std::map<std::string, std::map<std::string, std::string> >::const_iterator it;
	if ((it=id2group_.find(species)) != id2group_.end())
	{
		std::map<std::string, std::string>::const_iterator it2;
		if ((it2=it->second.find(seqID)) != it->second.end())
			return &groups_.at(it2->second);
	}
	return nullptr;
}

} /* namespace BioSeqDataLib */
