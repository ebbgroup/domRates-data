/*
 * FeatureTest.hpp
 *
 *  Created on: 1 Dec 2013
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


#ifndef FEATURETEST_HPP_
#define FEATURETEST_HPP_

#include <limits>
#include <vector>

#include "../../src/sequence/Sequence.hpp"
#include "../../src/annotation/FeatureSet.hpp"


BOOST_AUTO_TEST_SUITE(FeatureSet_Test)

BOOST_AUTO_TEST_CASE( Feature_Constructor_Test )
{
	std::ifstream gffStream("../tests/annotation/data/Test.gff");
	std::string line;
	std::getline(gffStream, line);
	std::getline(gffStream, line);
	BioSeqDataLib::Feature gene(line);
	std::getline(gffStream, line);
	BioSeqDataLib::Feature feature(line);

	//pbar_scf7180000350264	maker	mRNA	165232	175342	.	+	.	ID=PB18752-RA;Parent=PB18752;Name=PB18752-RA;_AED=0.39;_eAED=0.41;_QI=0|0|0|1|1|1|12|0|770;score=173.64;Alias=fgenesh_masked-pbar_scf7180000350264-abinit-gene-1.45-mRNA-1,PB18752-RA
	BOOST_CHECK_EQUAL(feature.seqId_, "pbar_scf7180000350264");
	BOOST_CHECK_EQUAL(feature.source_, "maker");
	BOOST_CHECK_EQUAL(feature.type_, "mRNA");
	BOOST_CHECK_EQUAL(feature.start_, 165231);
	BOOST_CHECK_EQUAL(feature.end_, 175341);
	BOOST_CHECK_EQUAL(feature.strand_, '+');
	BOOST_CHECK_CLOSE(feature.score_, std::numeric_limits<double>::max(), 0.0001);
	BOOST_CHECK_EQUAL(feature.phase_, -1);
	BOOST_CHECK_EQUAL(feature.parents_[0], "PB18752");
	BOOST_CHECK_EQUAL(feature.attribute("_AED"), "0.39");


	feature.attribute("_AED", "AHA");
	BOOST_CHECK_EQUAL(feature.attribute("_AED"), "AHA");

	std::getline(gffStream, line);
	BioSeqDataLib::Feature feature2(line);
	BOOST_CHECK_EQUAL(feature<feature2, true);
	BOOST_CHECK_EQUAL(feature>feature2, false);
	//BOOST_CHECK_EQUAL(feature<=feature, true);
	//BOOST_CHECK_EQUAL(feature>=feature, true);
	BOOST_CHECK_EQUAL(feature>feature2, false);
	BOOST_CHECK_EQUAL(feature<feature2, true);
	//BOOST_CHECK_EQUAL(gene>feature, true);
	BOOST_CHECK_EQUAL(feature<gene, false);
	BOOST_CHECK_EQUAL(gene>feature, false);
	//BOOST_CHECK_EQUAL(feature>gene, true);
	BOOST_CHECK_EQUAL(feature==feature2, false);
	BOOST_CHECK_EQUAL(feature==feature, true);
	BOOST_CHECK_EQUAL(feature!=feature2, true);
	BOOST_CHECK_EQUAL(feature2!=feature2, false);
}

/*
BOOST_AUTO_TEST_CASE( FeatureSet_Test )
{
	BioSeqDataLib::FeatureSet set;
	set.readGff("../tests/Annotation/data/Test.gff");
	auto x = set["PB18752-RA"];
	BOOST_CHECK_EQUAL(x->start_, 165231);
	++x;
	BOOST_CHECK_EQUAL(x->attribute("ID"), "PB18752-RA:cds1");

	const BioSeqDataLib::FeatureSet &set2 = set;
	auto y = set2["augustus_masked-pbar_scf7180000348142-processed-gene-0.3-mRNA-1:exon:1389"];
	BOOST_CHECK_EQUAL(y->start_, 182);

	BOOST_CHECK_EQUAL(set.counts("mRNA"), 4);
}*/


BOOST_AUTO_TEST_CASE( SimpleRead_Test )
{
	auto features = BioSeqDataLib::simpleGffRead("../tests/annotation/data/Test.gff", "mRNA", true);
	BOOST_CHECK_EQUAL(features["pbar_scf7180000350264"].size(), 2);
	BOOST_CHECK_EQUAL(features["pbar_scf7180000350264"][0].attribute("ID"), "PB18753-RA");
	BOOST_CHECK_EQUAL(features["pbar_scf7180000348142"][0].attribute("ID"), "test1");
	BOOST_CHECK_EQUAL(features.size(), 3);

	features = BioSeqDataLib::simpleGffRead("../tests/annotation/data/Test.gff", "mRNA", false);
	BOOST_CHECK_EQUAL(features["pbar_scf7180000350264+"].size(), 2);
	BOOST_CHECK_EQUAL(features["pbar_scf7180000348142-"][0].attribute("ID"), "test2");
	BOOST_CHECK_EQUAL(features.size(), 4);

}

BOOST_AUTO_TEST_SUITE_END()

#endif /* FEATURETEST_HPP_ */
