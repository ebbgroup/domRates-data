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
#include <boost/test/output_test_stream.hpp>
#include <iostream>


#include "../../src/sequence/Sequence.hpp"
#include "../../src/sequence/SequenceSet.hpp"
#include "../../src/sequence/SeqFunctions.hpp"
#include "../../src/domain/PfamDomain.hpp"
#include "../../src/domain/DomainArrangement.hpp"


BOOST_AUTO_TEST_SUITE(Sequence_Test)


BOOST_AUTO_TEST_CASE( Sequence_Test )
{
	BioSeqDataLib::Sequence<> seq1("seq1", "ACGTCT", "B7", "test sequence");
	BOOST_CHECK_EQUAL(seq1.name(), "seq1");
	BOOST_CHECK_EQUAL(seq1.seq(), "ACGTCT");
	BOOST_CHECK_EQUAL(seq1.accession(), "B7");
	BOOST_CHECK_EQUAL(seq1.comment(), "test sequence");
	BOOST_CHECK_EQUAL(seq1.size(), 6);
	BioSeqDataLib::Sequence<> seq2("seq2", "TTTTTT", "", "test sequence");
	seq1=seq2;
	BOOST_CHECK_EQUAL(seq1, seq2);
	seq1=std::move(seq2);
	BOOST_CHECK_EQUAL(seq1.name(), "seq2");
	//BOOST_CHECK_EQUAL(seq2.name(), "");

	seq1.append("ATTT");
	BOOST_CHECK_EQUAL(seq1.seq(), "TTTTTTATTT");

	seq1.append('C');
	BOOST_CHECK_EQUAL(seq1.seq(), "TTTTTTATTTC");

	seq1.resize(20);
	BOOST_CHECK_EQUAL(seq1.size(), 20);
	seq1.resize(3);
	BOOST_CHECK_EQUAL(seq1.seq(), "TTT");
}


BOOST_AUTO_TEST_CASE( Sequence_Functions_Test)
{
	BioSeqDataLib::Sequence<> seq1("seq1", "ACgtCATac", "", "test sequence");
	reverseComplement(seq1);
	BOOST_CHECK_EQUAL(seq1.seq(), "gtATGacGT");

	BioSeqDataLib::toUpper(seq1);
	BOOST_CHECK_EQUAL(seq1.seq(), "GTATGACGT");

	BioSeqDataLib::toLower(seq1);
	BOOST_CHECK_EQUAL(seq1.seq(), "gtatgacgt");

	BioSeqDataLib::replace(seq1, 'a', 'g');
	BOOST_CHECK_EQUAL(seq1.seq(), "gtgtggcgt");

	BioSeqDataLib::Sequence<> subSeq = BioSeqDataLib::subseq(seq1, 2, 3);
	BOOST_CHECK_EQUAL(seq1.seq(), "gtgtggcgt");
	BOOST_CHECK_EQUAL(subSeq.seq(), "gt");
	BOOST_CHECK_EQUAL(subSeq.name(), "seq1");
	BOOST_CHECK_EQUAL(subSeq.comment(), "");

	BioSeqDataLib::Sequence<> seq2("seq1", "XATGTTTAATTAATTAATTATGXXXXXXTTTAATTAATTAATT", "", "test sequence");
	auto pair = longestOrf(seq2, {"ATG"}, {"TAA"}, 1);
	BOOST_CHECK_EQUAL(pair.first, 19);
	BOOST_CHECK_EQUAL(pair.second, 36);

	BioSeqDataLib::Sequence<> seq3("seq1", "AATTAATTAATTAAACATX", "", "test sequence");
	pair = longestOrf(seq3, {"ATG"}, {"TAA"}, 1);
	BOOST_CHECK_EQUAL(pair.first, -1);
	BOOST_CHECK_EQUAL(pair.second, -12);

}

BOOST_AUTO_TEST_CASE(Sequence_analysis_Test)
{
	// hydropathy score
	BioSeqDataLib::Sequence<> seq("seq1", "AWraFA", "", "test sequence");
	std::vector<float> scores= BioSeqDataLib::calcHydropathy(seq, 3);
	BOOST_CHECK_EQUAL(scores.size(), 4);
	BOOST_CHECK_CLOSE( -1.2, scores[0], 0.0001);
	BOOST_CHECK_CLOSE( -1.2, scores[1], 0.0001);
	BOOST_CHECK_CLOSE( 0.0333333333, scores[2], 0.0001);
	BOOST_CHECK_CLOSE( 2.1333333333, scores[3], 0.0001);

	scores= BioSeqDataLib::calcHydropathy(seq, 1);
	BOOST_CHECK_EQUAL(scores.size(), 6);
	BOOST_CHECK_CLOSE( 1.8, scores[3], 0.0001);
	BOOST_CHECK_CLOSE( 1.8, scores[0], 0.0001);
	BOOST_CHECK_CLOSE( 1.8, scores[5], 0.0001);

	// CpG o/e
	BioSeqDataLib::Sequence<> seq2("seq2", "AGCTCTACG", "", "test sequence");
	float cpgoe = BioSeqDataLib::calcCpGoe(seq2, true);
	BOOST_CHECK_EQUAL(cpgoe, 1.6875);
	cpgoe = BioSeqDataLib::calcCpGoe(seq2, false);
	BOOST_CHECK_EQUAL(cpgoe, 1.6875);
	BioSeqDataLib::Sequence<> seq3("seq3", "NGCTCNACG", "", "test sequence");
	cpgoe = BioSeqDataLib::calcCpGoe(seq3.seq(), false);
	BOOST_CHECK_CLOSE(cpgoe, 1.6875, 0.0001);
	cpgoe = BioSeqDataLib::calcCpGoe(seq3, true, 0.5);
	BOOST_CHECK_CLOSE(cpgoe, 1.633333, 0.0001);

	// CG content
	BioSeqDataLib::TwoValues cgFrac = BioSeqDataLib::calcCGcontent(seq2, true);
	BOOST_CHECK_EQUAL(cgFrac.first(), 5);
	BOOST_CHECK_EQUAL(cgFrac.second(), 9);
	cgFrac = BioSeqDataLib::calcCGcontent(seq2, false);
	BOOST_CHECK_EQUAL(cgFrac.first(), 5);
	BOOST_CHECK_EQUAL(cgFrac.second(), 9);
	cgFrac = BioSeqDataLib::calcCGcontent(seq3.seq(), false);
	BOOST_CHECK_EQUAL(cgFrac.first(), 5);
	BOOST_CHECK_EQUAL(cgFrac.second(), 9);
	cgFrac = BioSeqDataLib::calcCGcontent(seq3, true);
	BOOST_CHECK_EQUAL(cgFrac.first(), 5);
	BOOST_CHECK_EQUAL(cgFrac.second(), 7);

}


BOOST_AUTO_TEST_CASE( Sequence_extended_Test)
{
	BioSeqDataLib::Sequence<BioSeqDataLib::DomainArrangementWrapper<BioSeqDataLib::PfamDomain> > seq1("seq1", "ACgtCATac", "", "test sequence");
	reverseComplement(seq1);
	BOOST_CHECK_EQUAL(seq1.seq(), "gtATGacGT");
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::PfamDomain> &arrangement = seq1.domainArrangement();
	BOOST_CHECK_EQUAL(arrangement.size(), 0);

	BioSeqDataLib::Sequence<> seq2("seq2", "GCTGCCGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAG", "", "test sequence");
	BioSeqDataLib::translate(seq2);
	BOOST_CHECK_EQUAL(seq2.seq(), "AAAALLLLLLRRRRRRKKNNMDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV***");

	BioSeqDataLib::Sequence<> seq3("seq2", "GCUGCCGCAGCGUUAUUGCUUCUCCUACUGCGUCGCCGACGGAGAAGGAAAAAGAAUAACAUGGAUGACUUUUUCUGUUGCCCUCCCCCACCGCAACAGUCUUCCUCAUCGAGUAGCGAAGAGACUACCACAACGGGUGGCGGAGGGUGGCAUCACUAUUACAUUAUCAUAGUUGUCGUAGUGUAAUGAUAG", "", "test sequence");
	BioSeqDataLib::translate(seq3);
	BOOST_CHECK_EQUAL(seq3.seq(), "AAAALLLLLLRRRRRRKKNNMDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV***");
	BOOST_CHECK_EQUAL(seq3.name(), seq3.name());
	BOOST_CHECK_EQUAL(seq3.comment(), seq3.comment());

	BioSeqDataLib::Sequence<BioSeqDataLib::DomainArrangementWrapper<BioSeqDataLib::PfamDomain> > seq4("seq1", "cCN", "", "test sequence");
	BioSeqDataLib::translate(seq4);
	BOOST_CHECK_EQUAL(seq4.seq(), "P");
}

struct cout_redirect {
    cout_redirect( std::streambuf * new_buffer )
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

private:
    std::streambuf * old;
};

BOOST_AUTO_TEST_CASE( Sequence_output_Test)
{
	BioSeqDataLib::Sequence<> seq1("seq1", "HELLO", "", "");

	boost::test_tools::output_test_stream output;
	{
		cout_redirect guard( output.rdbuf( ) );
		std::cout << seq1 << std::endl;
	}
	BOOST_CHECK( output.is_equal( ">seq1\nHELLO\n") );
}


BOOST_AUTO_TEST_SUITE_END()
