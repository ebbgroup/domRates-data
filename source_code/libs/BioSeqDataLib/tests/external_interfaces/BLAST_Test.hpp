/*
 * Sequence_Test.cpp
 *
 *  Created on: 15 Oct 2013
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


#include <boost/test/unit_test.hpp>
#include <iostream>


#include "../../src/external_interfaces/runDBsearch.hpp"
#include "../../src/annotation/BlastHit.hpp"
#include "../../src/external_interfaces/runPhylo.hpp"
#include "../../src/external_interfaces/runDBsearch.hpp"
#include "../../src/external_interfaces/domainProgs.hpp"
#include "../../src/utility/Settings.hpp"
#include "../../src/domain/DomainArrangementSet.hpp"
#include "../../src/sequence/SequenceSet.hpp"

BOOST_AUTO_TEST_SUITE(DBSearch_Test)


BOOST_AUTO_TEST_CASE( Blast_Test )
{
	BioSeqDataLib::makeBlastDB("../tests/sequence/data/seqSet.fasta", true);
	BioSeqDataLib::runblast("blastp", "../tests/sequence/data/seqSet.fasta", "../tests/sequence/data/seqSet.fasta", "aha", "log");
	BOOST_CHECK_EQUAL(1,1);
}


BOOST_AUTO_TEST_CASE( Alignment_Tests )
{
	BioSeqDataLib::runAlignment("../tests/external_interfaces/data/hem.fasta", BioSeqDataLib::AlnProgram::tcoffee, "tcoffee.fasta", 1);
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<>> set("tcoffee.fasta");
	BOOST_CHECK_EQUAL(set.size(), 4);
	set.clear();
	BioSeqDataLib::runAlignment("../tests/external_interfaces/data/hem.fasta", BioSeqDataLib::AlnProgram::mafft, "mafft.aln", 1, BioSeqDataLib::AlnFormat::clustal);
	// No read clustal format yet
	BioSeqDataLib::runAlignment("../tests/external_interfaces/data/hem.fasta", BioSeqDataLib::AlnProgram::mafft, "mafft.fa", 1, BioSeqDataLib::AlnFormat::fasta);
	set.read("mafft.fa");
	BOOST_CHECK_EQUAL(set.size(), 4);
	set.clear();
	BioSeqDataLib::runAlignment("../tests/external_interfaces/data/hem.fasta", BioSeqDataLib::AlnProgram::clustalo, "clustalo.fasta", 1);
	set.read("clustalo.fasta");
	BOOST_CHECK_EQUAL(set.size(), 4);
	set.clear();
	BioSeqDataLib::runAlignment("../tests/external_interfaces/data/hem.fasta", BioSeqDataLib::AlnProgram::msaprobs, "msaProbs.fasta", 1);
	set.read("msaProbs.fasta");
	BOOST_CHECK_EQUAL(set.size(), 4);
	set.clear();
}




BOOST_AUTO_TEST_SUITE_END()
