/*
 * Alignment_Test.hpp
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
#ifndef ALIGNMENT_TEST_HPP_
#define ALIGNMENT_TEST_HPP_


#include <boost/test/unit_test.hpp>
#include <iostream>


#include "../../src/sequence/Alignment.hpp"
#include "../../src/sequence/AlignmentFunctions.hpp"
#include "../../src/sequence/Sequence.hpp"


BOOST_AUTO_TEST_SUITE(Alignment_Test)


BOOST_AUTO_TEST_CASE( AlignmentRead_Test )
{
	BioSeqDataLib::Alignment<BioSeqDataLib::Sequence<> > aln1("../tests/sequence/data/seqSet2.msf");
	BOOST_CHECK_EQUAL(aln1[0].seq(), "MGTVLSLSPSYR-------GAATVGHYTAVQNSKNAKDKNLKRHSIISVLPWKRIVAVSAKKKNSKKVQPNSSYQNNITHLNNENLKKSLSCANLSTFAQ");

	BioSeqDataLib::Alignment<BioSeqDataLib::Sequence<> > aln;
	aln.read("../tests/sequence/data/seqSet2.msf");
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set = aln.seqSet();
	BOOST_CHECK_EQUAL(set[0].seq(), "MGTVLSLSPSYRGAATVGHYTAVQNSKNAKDKNLKRHSIISVLPWKRIVAVSAKKKNSKKVQPNSSYQNNITHLNNENLKKSLSCANLSTFAQ");
}

BOOST_AUTO_TEST_CASE( Alignment_write_Test )
{
	BioSeqDataLib::Alignment<BioSeqDataLib::Sequence<> > aln("../tests/sequence/data/seqSet2.msf");
	aln.write("test_fasta_out_aln.fa","fasta");
	std::ifstream ifs1_1("../tests/sequence/data/seqSet2.fa");
	std::ifstream ifs2_1("test_fasta_out_aln.fa");
	std::istream_iterator<char> b1_1(ifs1_1), e1_1;
	std::istream_iterator<char> b2_1(ifs2_1), e2_1;
	BOOST_CHECK_EQUAL_COLLECTIONS(b1_1, e1_1, b2_1, e2_1);

	aln.write("test_msf_out_aln.msf","msf");
	std::ifstream ifs1_2("../tests/sequence/data/msf_out.msf");
	std::ifstream ifs2_2("test_msf_out_aln.msf");
	std::istream_iterator<char> b1_2(ifs1_2), e1_2;
	std::istream_iterator<char> b2_2(ifs2_2), e2_2;
	BOOST_CHECK_EQUAL_COLLECTIONS(b1_2, e1_2, b2_2, e2_2);

	aln.write("test_phylip_out_aln.phylips","phylips");
	std::ifstream ifs1_3("../tests/sequence/data/phylips_out.phy");
	std::ifstream ifs2_3("test_phylip_out_aln.phylips");
	std::istream_iterator<char> b1_3(ifs1_3), e1_3;
	std::istream_iterator<char> b2_3(ifs2_3), e2_3;
	BOOST_CHECK_EQUAL_COLLECTIONS(b1_3, e1_3, b2_3, e2_3);
}



BOOST_AUTO_TEST_CASE( AlignmentFunction_Test )
{
	BioSeqDataLib::Alignment<BioSeqDataLib::Sequence<> > aln;
	aln.read("../tests/sequence/data/seqSet2.msf");
	auto err = BioSeqDataLib::alnCheck(aln);
	BOOST_CHECK_EQUAL(err.first(), 0);

	BioSeqDataLib::Alignment<BioSeqDataLib::Sequence<> > aln1;
	aln1.read("../tests/sequence/data/seqSet.fasta");
	auto err1 = BioSeqDataLib::alnCheck(aln1);
	BOOST_CHECK_EQUAL(err1.first(), 1);
	BOOST_CHECK_EQUAL(err1.second(), 1);

	BioSeqDataLib::Alignment<BioSeqDataLib::Sequence<> > aln2;
	aln2.read("../tests/sequence/data/seqSet3.msf");
	auto err2 = BioSeqDataLib::alnCheck(aln2);
	BOOST_CHECK_EQUAL(err2.first(), 2);
	BOOST_CHECK_EQUAL(err2.second(), 5);


	BioSeqDataLib::shrinkToSubAln(aln2, 0, 19);
	BOOST_CHECK_EQUAL(aln2[0].seq(), "MGTVL-LSPSYR-------G");
	BioSeqDataLib::shrinkToSubAln(aln2, 5, 9);
	BOOST_CHECK_EQUAL(aln2[0].seq(), "-LSPS");

	BioSeqDataLib::Alignment<BioSeqDataLib::Sequence<> > aln3;
	aln3.read("../tests/sequence/data/seqSet3.msf");
	BioSeqDataLib::shrinkToSubAln(aln3, 0, 19, "seq1");
	BOOST_CHECK_EQUAL(aln3[0].seq(), "MGTVL-LSPSYR-------GAATVGHYT");
	BioSeqDataLib::shrinkToSubAln(aln3, 2, 9, "seq2");
	BOOST_CHECK_EQUAL(aln3[1].seq(), "TVL-LSPSY");

	BioSeqDataLib::Alignment<BioSeqDataLib::Sequence<> > alnIdentity;
	alnIdentity.push_back(BioSeqDataLib::Sequence<>("Seq1","ACGt", "", "test sequence"));
	alnIdentity.push_back(BioSeqDataLib::Sequence<>("Seq2","AcGt", "", "test sequence"));
	alnIdentity.push_back(BioSeqDataLib::Sequence<>("Seq3","A-et", "", "test sequence"));
	BioSeqDataLib::TwoValues vals = BioSeqDataLib::identity(alnIdentity);
	BOOST_CHECK_EQUAL(vals.first(), 8);
	BOOST_CHECK_EQUAL(vals.second(), 10);
}

BOOST_AUTO_TEST_SUITE_END()



#endif /* ALIGNMENT_TEST_HPP_ */
