/*
 * Feature.cpp
 *
 *  Created on: 3 Jul 2014
 *	Copyright: 2014
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

#include "Feature.hpp"

namespace BioSeqDataLib {

using namespace std;

Feature::Feature(const std::string &line)
{
	vector<string> tokens;
	vector<std::string> attributeSplit;
	tokens=split(line, "\t");
	attributeSplit= split(tokens[tokens.size()-1], ";");

	seqId_ = tokens[0];
	source_ = tokens[1];
	type_ = tokens[2];
	start_ = stoul(tokens[3])-1;
	end_ = stoul(tokens[4])-1;
	strand_ = tokens[6][0];

	if (tokens[5][0]=='.')
		score_=std::numeric_limits<double>::max();
	else
		score_ = stod(tokens[5]);
	if (tokens[7][0]=='.')
		phase_=-1;
	else
		phase_=stoi(tokens[7]);

	for (size_t i=0; i<attributeSplit.size(); ++i)
	{
		size_t pos = attributeSplit[i].find('=');
		if ((pos==0) || (pos==(attributeSplit[i].size()-1)))
			continue;
		string key = attributeSplit[i].substr(0,pos);
		if (key == "Parent")
		{
			vector<string> parents;
			split(attributeSplit[i].substr(pos+1), ",", parents);
			for (size_t j=0; j<parents.size(); ++j)
				parents_.push_back(parents[j]);
		}
		else
			attributes_[key] = attributeSplit[i].substr(pos+1);
	}

}


bool
operator<(const Feature &a, const Feature &b)
{
	if (b.seqId_ == a.seqId_)
	{
		//if ((!b.parents_.empty()) && (b.parents_.at(0) == a.attributes_.at("ID")))
		//	return true;
		if (b.start_ == a.start_)
				return a.end_ > b.end_;
		else
			return a.start_ < b.start_;
	}
	else
		return a.seqId_ < b.seqId_;
}


bool
operator>(const Feature &a, const Feature &b)
{
	if (a.seqId_ == b.seqId_)
	{
	//	if ((!a.parents_.empty()) && (a.parents_.at(0) == b.attributes_.at("ID")))
	//		return true;
		if (a.start_ == b.start_)
			return a.end_ < b.end_;
		else
			return a.start_ > b.start_;
	}
	else
		return a.seqId_ > b.seqId_;
}

/**
 * \brief Simple reading of a single feature type
 * @param inFile The gff file.
 * @param featType The feature type to use.
 * @return All features grouped by sequence.
 * \relates Feature
 */
std::map<std::string, std::vector<Feature> >
simpleGffRead(const std::string &inFile, const std::string &featType, bool ignoreStrand)
{
	std::map<std::string, std::vector<Feature> > features;
	std::ifstream inS(inFile);
	string line;

	if (ignoreStrand)
	{
		while (getline(inS, line))
		{
			if (line.empty())
				continue;
			if (line.compare("##FASTA") == 0)
				break;
			if (line[0] == '#')
				continue;

			Feature feature(line);
			if (feature.type_ == featType)
				features[feature.seqId_].emplace_back(std::move(feature));
		}
	}
	else
	{
		while (getline(inS, line))
		{
			if (line.empty())
				continue;
			if (line.compare(0,7,"##FASTA") == 0)
				break;
			if (line[0] == '#')
				continue;
			Feature feature(line);

			if (feature.type_ == featType)
				features[feature.seqId_ + (feature.strand_=='+' ? "+" : "-" )].emplace_back(std::move(feature));
		}
	}

	for (auto &pair : features)
		sort(pair.second.begin(), pair.second.end());

	return features;
}

} /* namespace BioSeqDataLib */
