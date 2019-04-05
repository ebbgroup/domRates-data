/*
 * FeatureSet.cpp
 *
 *  Created on: 29 Nov 2013
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

#include "FeatureSet.hpp"
#include <iostream>
using namespace std;

namespace BioSeqDataLib
{

FeatureSet::FeatureSet()
{
	// TODO Auto-generated constructor stub

}

FeatureSet::~FeatureSet()
{
	// TODO Auto-generated destructor stub
}



void
FeatureSet::readGff(const std::string &inFile)
{
	ifstream gffFile;
	gffFile.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
	gffFile.open(inFile);
	gffFile.exceptions(std::ifstream::badbit);
	std::string line;
	getline(gffFile, line);
	if (line != "##gff-version 3")
	{
		const std::string errorMessage = "The file " + inFile + "does not contain the gff version header 3";
		throw FormatException(errorMessage);
	}

	//pbar_scf7180000350264	maker	CDS	169473	169514	.	+	0	ID=PB18752-RA:cds;Parent=PB18752-RA
	while (getline(gffFile, line))
	{
		if ((line.empty()) || (line[0] == '#'))
			continue;

		Feature feature(line);
		id2position_[feature.attribute("ID")] = feature.seqId_;
		typeCounter[feature.type_] += 1;
		features_[feature.seqId_].push_back(std::move(feature));

	}

	for (auto &tmp : features_)
		sort(tmp.second.begin(), tmp.second.end());
	/*
	##gff-version 3
	pbar_scf7180000350264	maker	gene	165232	175342	.	+	.	ID=PB18752;Name=PB18752;score=173.64
	pbar_scf7180000350264	maker	mRNA	165232	175342	.	+	.	ID=PB18752-RA;Parent=PB18752;Name=PB18752-RA;_AED=0.39;_eAED=0.41;_QI=0|0|0|1|1|1|12|0|770;score=173.64;Alias=fgenesh_masked-pbar_scf7180000350264-abinit-gene-1.45-mRNA-1,PB18752-RA
	pbar_scf7180000350264	maker	exon	165232	165313	.	+	.	ID=PB18752-RA:exon:0;Parent=PB18752-RA
	pbar_scf7180000350264	maker	exon	168731	168972	.	+	.	ID=PB18752-RA:exon:1;Parent=PB18752-RA
	*/
}


} /* namespace BioSeqDataLib */
