#ifndef DOMAIN_ARRANGEMENT_SET_TEST_HPP_
#define DOMAIN_ARRANGEMENT_SET_TEST_HPP_


#include <boost/test/unit_test.hpp>
#include <iostream>


#include "../../src/sequence/Sequence.hpp"
#include "../../src/sequence/SequenceSet.hpp"
#include "../../src/sequence/SeqFunctions.hpp"
#include "../../src/DomainModule.hpp"
#include "../../src/domain/DomainArrangementSet.hpp"


BOOST_AUTO_TEST_SUITE(DomainArrangementSet_Test)


/**********************************************************
 *               READ FUNCTION TESTS                      *
 **********************************************************/


BOOST_AUTO_TEST_CASE( DomainArrangement_readHmmScan_Test )
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
	arrangementSet.read("../tests/domain/data/BB20012.hmmScan");

	BOOST_CHECK_EQUAL(arrangementSet.size(), 27);
	//GTP_EFTU             PF00009.22   188 IF2G_METTH           -            408   5.3e-44  149.9   0.2   1   1   7.1e-47   8.1e-44  149.3   0.2     3   187     6   198     4   199 0.95 Elongation factor Tu GTP binding domain
	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::PfamDomain> &set = arrangementSet["IF2G_METTH"];
	BOOST_CHECK_EQUAL(set.size(), 17);
	const BioSeqDataLib::PfamDomain &dom = set[1];
	BOOST_CHECK_EQUAL(dom.name(), "GTP_EFTU");
	BOOST_CHECK_EQUAL(dom.accession(), "PF00009");
	BOOST_CHECK_EQUAL(dom.start(), 5);
	BOOST_CHECK_EQUAL(dom.end(), 197);
	BOOST_CHECK_EQUAL(dom.envStart(), 3);
	BOOST_CHECK_EQUAL(dom.envEnd(), 198);
	BOOST_CHECK_EQUAL(dom.hmmStart(), 2);
	BOOST_CHECK_EQUAL(dom.hmmEnd(), 186);
	BOOST_CHECK_EQUAL(dom.hmmLength(), 188);
	BOOST_CHECK_EQUAL(dom.evalue(), 8.1e-44);
	BOOST_CHECK_CLOSE(dom.bitscore(), 149.3, 0.01);
	BOOST_CHECK_EQUAL(dom.type(), "");


	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::DomainExt> arrangementSet2;
	arrangementSet2.read("../tests/domain/data/BB20012.hmmScan");
	BOOST_CHECK_EQUAL(arrangementSet2.size(), 27);
	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::DomainExt> &set2 = arrangementSet2["IF2G_METTH"];
	BOOST_CHECK_EQUAL(set2.size(), 17);
	const BioSeqDataLib::DomainExt &dom2 = set2[1];
	BOOST_CHECK_EQUAL(dom2.name(), "GTP_EFTU");
	BOOST_CHECK_EQUAL(dom2.accession(), "PF00009");
	BOOST_CHECK_EQUAL(dom2.start(), 5);
	BOOST_CHECK_EQUAL(dom2.end(), 197);
	BOOST_CHECK_EQUAL(dom2.envStart(), 3);
	BOOST_CHECK_EQUAL(dom2.envEnd(), 198);
	BOOST_CHECK_EQUAL(dom2.hmmStart(), 2);
	BOOST_CHECK_EQUAL(dom2.hmmEnd(), 186);
	BOOST_CHECK_EQUAL(dom2.hmmLength(), 188);
	BOOST_CHECK_EQUAL(dom2.evalue(), 8.1e-44);
	BOOST_CHECK_CLOSE(dom2.bitscore(), 149.3, 0.01);



	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::Domain> arrangementSet3;
	arrangementSet3.read("../tests/domain/data/BB20012.hmmScan");
	BOOST_CHECK_EQUAL(arrangementSet3.size(), 27);
	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> &set3 = arrangementSet3["IF2G_METTH"];
	BOOST_CHECK_EQUAL(set3.size(), 17);
	const BioSeqDataLib::Domain &dom3 = set3[1];
	BOOST_CHECK_EQUAL(dom3.accession(), "PF00009");
	BOOST_CHECK_EQUAL(dom3.start(), 5);
	BOOST_CHECK_EQUAL(dom3.end(), 197);
	BOOST_CHECK_EQUAL(dom3.evalue(), 8.1e-44);


	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<BioSeqDataLib::DomainArrangementWrapper<BioSeqDataLib::PfamDomain> > > seqSet;
	seqSet.read("../tests/domain/data/BB20012.fa");
	arrangementSet.moveTo(seqSet);

	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::PfamDomain> &arrangement = seqSet[0].domainArrangement();
	BOOST_CHECK_EQUAL(arrangement.size(), 32);
	BOOST_CHECK_EQUAL(arrangement[1].name(), "FeoB_N");
}



BOOST_AUTO_TEST_CASE( DomainArrangement_readASS_Test )
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::SFDomain> arrangementSet;
	arrangementSet.read("../tests/domain/data/test.ass");

	BOOST_CHECK_EQUAL(arrangementSet.size(), 2);

	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::SFDomain> &set = arrangementSet["gnl|Cobs_1.4|Cobs_00021-mRNA-1"];
	BOOST_CHECK_EQUAL(set.size(), 2);
	const BioSeqDataLib::SFDomain &dom = set[1];
	BOOST_CHECK_EQUAL(dom.name(), "0045850");
	BOOST_CHECK_EQUAL(dom.accession(), "55074");
	BOOST_CHECK_EQUAL(dom.start(), 1130);
	BOOST_CHECK_EQUAL(dom.end(), 1228);
	BOOST_CHECK_EQUAL(dom.envStart(), 1130);
	BOOST_CHECK_EQUAL(dom.envEnd(), 1228);
	BOOST_CHECK_EQUAL(dom.hmmStart(), 1);
	BOOST_CHECK_EQUAL(dom.hmmEnd(), 0);
	BOOST_CHECK_EQUAL(dom.hmmLength(), 0);
	BOOST_CHECK_EQUAL(dom.evalue(), 2.35e-15);
	BOOST_CHECK_CLOSE(dom.bitscore(), -1, 0.01);
	BOOST_CHECK_EQUAL(dom.modelID(), "0045850");
	BOOST_CHECK_EQUAL(dom.scopDomID(), "55074");
	BOOST_CHECK_EQUAL(dom.scopID(), "39423");
	//gnl|Cobs_1.4|Cobs_00021-mRNA-1	0045850	1131-1229	2.35e-15	2	QSRDKVGVMFASVPNFTEFYSED--VNKGMECIRLLNEIIADFDELLDETAFHCIEKIKTVGATYMAASGLNPSQTDNNGDDM----EHLCKLVDYAVAMRHRLE	9.09e-05	39423	55074

	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::DomainExt> arrangementSet2;
	arrangementSet2.read("../tests/domain/data/test.ass");
	BOOST_CHECK_EQUAL(arrangementSet2.size(), 2);
	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::DomainExt> &set2 = arrangementSet2["gnl|Cobs_1.4|Cobs_00021-mRNA-1"];
	BOOST_CHECK_EQUAL(set2.size(), 2);
	const BioSeqDataLib::DomainExt &dom2 = set2[1];
	BOOST_CHECK_EQUAL(dom2.name(), "0045850");
	BOOST_CHECK_EQUAL(dom2.accession(), "55074");
	BOOST_CHECK_EQUAL(dom2.start(), 1130);
	BOOST_CHECK_EQUAL(dom2.end(), 1228);
	BOOST_CHECK_EQUAL(dom2.envEnd(), 1228);
	BOOST_CHECK_EQUAL(dom2.envStart(), 1130);
	BOOST_CHECK_EQUAL(dom2.hmmStart(), 1);
	BOOST_CHECK_EQUAL(dom2.hmmEnd(), 0);
	BOOST_CHECK_EQUAL(dom2.hmmLength(), 0);
	BOOST_CHECK_EQUAL(dom2.evalue(), 2.35e-15);
	BOOST_CHECK_CLOSE(dom2.bitscore(), -1, 0.01);

	//gnl|Cobs_1.4|Cobs_00021-mRNA-1	0045850	1131-1229	2.35e-15	2	QSRDKVGVMFASVPNFTEFYSED--VNKGMECIRLLNEIIADFDELLDETAFHCIEKIKTVGATYMAASGLNPSQTDNNGDDM----EHLCKLVDYAVAMRHRLE	9.09e-05	39423	55074

	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::Domain> arrangementSet3;
	arrangementSet3.read("../tests/domain/data/test.ass");
	BOOST_CHECK_EQUAL(arrangementSet3.size(), 2);
	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> &set3 = arrangementSet3["gnl|Cobs_1.4|Cobs_00021-mRNA-1"];
	BOOST_CHECK_EQUAL(set3.size(), 2);
	const BioSeqDataLib::Domain &dom3 = set2[1];
	BOOST_CHECK_EQUAL(dom3.accession(), "55074");
	BOOST_CHECK_EQUAL(dom3.start(), 1130);
	BOOST_CHECK_EQUAL(dom3.end(), 1228);
	BOOST_CHECK_EQUAL(dom3.evalue(), 2.35e-15);
}


BOOST_AUTO_TEST_CASE( DomainArrangement_readXDOM_Test )
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::DomainExt> arrangementSet;
	arrangementSet.read("../tests/domain/data/test.xdom");

	BOOST_CHECK_EQUAL(arrangementSet.size(), 3);

	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::DomainExt> &da = arrangementSet["seq1"];
	BOOST_CHECK_EQUAL(da.size(), 2);

	const BioSeqDataLib::DomainExt &dom = da[1];
	BOOST_CHECK_EQUAL(dom.accession(), "PF00002");
	BOOST_CHECK_EQUAL(dom.start(), 19);
	BOOST_CHECK_EQUAL(dom.end(), 39);



}

BOOST_AUTO_TEST_CASE( DomainArrangement_readTsv_Test )
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::Domain> arrangementSet;
	arrangementSet.read("../tests/domain/data/test.tsv");

	BOOST_CHECK_EQUAL(arrangementSet.size(), 2);

	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> &da = arrangementSet["EFX64876"];
	BOOST_CHECK_EQUAL(da.size(), 7);
	//EFX64876	c8b7a5d76f40fa853647cdb7223c235c	104	PANTHER	PTHR24152		1	104	0.0	T	24-02-2015
	const BioSeqDataLib::Domain &dom = da[0];
	BOOST_CHECK_EQUAL(dom.accession(), "PTHR24152");
	BOOST_CHECK_EQUAL(dom.start(), 1);
	BOOST_CHECK_EQUAL(dom.end(), 104);
	BOOST_CHECK_EQUAL(true, dom.db() ==  BioSeqDataLib::DomainDB::panther);

	std::vector<BioSeqDataLib::DomainDB> dbImportance = {BioSeqDataLib::DomainDB::pfam};
	arrangementSet.solveDbOverlaps(dbImportance, 10, 0.10);
	BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> &da2 = arrangementSet["EFX64876"];
	BOOST_CHECK_EQUAL(da2.size(), 1);
}

BOOST_AUTO_TEST_CASE( DomainArrangement_readDAMA_Test )
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::Domain> arrangementSet;
	arrangementSet.read("../tests/domain/data/DAMAOut.txt");
	BOOST_CHECK_EQUAL(arrangementSet.size(), 9);

	auto &arrangement = arrangementSet["sp|P00504|AATC_CHICK"];
	BOOST_CHECK_EQUAL(arrangement.size(), 1);
	BOOST_CHECK_EQUAL(arrangement[0].accession(), "PF00155");
	auto &arrangement2 = arrangementSet["sp|Q5ZIM6|AATF_CHICK"];
	BOOST_CHECK_EQUAL(arrangement2[0].accession(), "PF13339");
	BOOST_CHECK_EQUAL(arrangement2[1].accession(), "PF08164");

}

BOOST_AUTO_TEST_CASE( DomainArrangement_readRADIANT_Test )
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::Domain> arrangementSet;
	arrangementSet.read("../tests/domain/data/RadiantOut.txt");
	BOOST_CHECK_EQUAL(arrangementSet.size(), 6);
	auto &arrangement = arrangementSet["ENSP00000362111.4"];
	BOOST_CHECK_EQUAL(arrangement.size(), 1);
	BOOST_CHECK_EQUAL(arrangement[0].accession(), "PF00335");
	auto &arrangement2 = arrangementSet["ENSP00000363117.3"];
	BOOST_CHECK_EQUAL(arrangement2[0].accession(), "PF00018");
	BOOST_CHECK_EQUAL(arrangement2[1].accession(), "PF00017");
	BOOST_CHECK_EQUAL(arrangement2[0].start(), 90);
	BOOST_CHECK_EQUAL(arrangement2[0].end(), 90);
}


BOOST_AUTO_TEST_CASE( DomainArrangement_readPfamScan_Test )
{
	// check read with PfamDomain
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
	arrangementSet.read("../tests/domain/data/BB20012.pfamScan");
	BOOST_CHECK_EQUAL(arrangementSet.size(), 27);
	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::PfamDomain> &set = arrangementSet["IF2G_ARCFU"];
	BOOST_CHECK_EQUAL(set.size(), 3);
	const BioSeqDataLib::PfamDomain &dom = set[1];
	BOOST_CHECK_EQUAL(dom.name(), "GTP_EFTU_D2");
	BOOST_CHECK_EQUAL(dom.accession(), "PF03144");
	BOOST_CHECK_EQUAL(dom.start(), 248);
	BOOST_CHECK_EQUAL(dom.end(), 326);
	BOOST_CHECK_EQUAL(dom.envStart(), 248);
	BOOST_CHECK_EQUAL(dom.envEnd(), 327);
	BOOST_CHECK_EQUAL(dom.hmmStart(), 0);
	BOOST_CHECK_EQUAL(dom.hmmEnd(), 72);
	BOOST_CHECK_EQUAL(dom.hmmLength(), 74);
	BOOST_CHECK_CLOSE(dom.evalue(), 0.0000000073, 0.000001);
	BOOST_CHECK_CLOSE(dom.bitscore(), 35.6, 0.01);
	BOOST_CHECK_EQUAL(dom.significance(), 1);
	BOOST_CHECK_EQUAL(dom.clan(), "CL0023");
	BOOST_CHECK_EQUAL(dom.type(), "Domain");
	//IF2G_ARCFU    249    327    249    328 PF03144.20  GTP_EFTU_D2       Domain     1    73    74     35.6   7.3e-09   1 CL0023
	BOOST_CHECK_EQUAL(set[2].clan(), "");

	// check read with Domain
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::Domain> arrangementSet2;
	arrangementSet2.read("../tests/domain/data/BB20012.pfamScan");
	BOOST_CHECK_EQUAL(arrangementSet2.size(), 27);
	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::Domain> &set2 = arrangementSet2["IF2G_ARCFU"];
	BOOST_CHECK_EQUAL(set2.size(), 3);
	const BioSeqDataLib::Domain &dom2 = set2[1];
	BOOST_CHECK_EQUAL(dom2.accession(), "PF03144");
	BOOST_CHECK_EQUAL(dom2.start(), 248);
	BOOST_CHECK_EQUAL(dom2.end(), 326);
	BOOST_CHECK_CLOSE(dom2.evalue(), 0.0000000073, 0.000001);

	// check read with DomainExt
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::DomainExt> arrangementSet3;
	arrangementSet3.read("../tests/domain/data/BB20012.pfamScan");
	BOOST_CHECK_EQUAL(arrangementSet3.size(), 27);
	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::DomainExt> &set3 = arrangementSet3["IF2G_ARCFU"];
	BOOST_CHECK_EQUAL(set.size(), 3);
	const BioSeqDataLib::DomainExt &dom3 = set3[1];
	BOOST_CHECK_EQUAL(dom3.name(), "GTP_EFTU_D2");
	BOOST_CHECK_EQUAL(dom3.accession(), "PF03144");
	BOOST_CHECK_EQUAL(dom3.start(), 248);
	BOOST_CHECK_EQUAL(dom3.end(), 326);
	BOOST_CHECK_EQUAL(dom3.envStart(), 248);
	BOOST_CHECK_EQUAL(dom3.envEnd(), 327);
	BOOST_CHECK_EQUAL(dom3.hmmStart(), 0);
	BOOST_CHECK_EQUAL(dom3.hmmEnd(), 72);
	BOOST_CHECK_EQUAL(dom3.hmmLength(), 74);
	BOOST_CHECK_CLOSE(dom3.evalue(), 0.0000000073, 0.000001);
	BOOST_CHECK_CLOSE(dom3.bitscore(), 35.6, 0.01);
}


BOOST_AUTO_TEST_CASE(ERROR_WHILE_READING_TEST)
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
	BOOST_CHECK_THROW(
		arrangementSet.read("../tests/domain/data/BB20012.pfamScan_format_error"), BioSeqDataLib::FormatException);

	BOOST_CHECK_THROW(arrangementSet.read("../tests/domain/data/does_not_exist"), std::exception);

	BOOST_CHECK_THROW(arrangementSet.read("../tests/domain/data/noKnownFormat.txt"), BioSeqDataLib::FormatException);
}


/**********************************************************
 *               WRITE FUNCTION TESTS                      *
 **********************************************************/

BOOST_AUTO_TEST_CASE(write_pfam_Test)
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
	arrangementSet.read("../tests/domain/data/BB20012.pfamScan");
	arrangementSet.write("test.txt", "pfam");
	auto nArrangements = arrangementSet.size();
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> testSet;
	testSet.read("test.txt");
	BOOST_CHECK_EQUAL(nArrangements, testSet.size());

	BioSeqDataLib::DomainArrangement<BioSeqDataLib::PfamDomain>  &set2 = arrangementSet["IF2G_ARCFU"];
	auto it = set2.begin();
	++it;
	set2.erase(it);
	BOOST_CHECK_EQUAL(set2.size(), 2);
}



BOOST_AUTO_TEST_CASE(write_xdom_Test)
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
	arrangementSet.read("../tests/domain/data/test.xdom");

	// check without sequence length
	arrangementSet.write("test.txt", "xdom");
	std::ifstream ifs1("../tests/domain/data/test_out.xdom");
    std::ifstream ifs2("test.txt");
    std::istream_iterator<char> b1(ifs1), e1;
    std::istream_iterator<char> b2(ifs2), e2;
    BOOST_CHECK_EQUAL_COLLECTIONS(b1, e1, b2, e2);


	// check with sequence length
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > seqSet;
	seqSet.push_back(BioSeqDataLib::Sequence<>("seq1", "AAAAAAAAAAAAAAAAAAAA", "", ""));
	seqSet.push_back(BioSeqDataLib::Sequence<>("seq2", "AAAAAAAAAAAAAAAAAAAAAAAA", "", ""));
	seqSet.push_back(BioSeqDataLib::Sequence<>("seq3", "AAAAAAAAAAAAAAAAAAAAAAA","", ""));

	arrangementSet.write("test.txt", "xdom", seqSet);
	std::ifstream ifs1_2("../tests/domain/data/test_out2.xdom");
    std::ifstream ifs2_2("test.txt");
    std::istream_iterator<char> b1_2(ifs1_2), e1_2;
    std::istream_iterator<char> b2_2(ifs2_2), e2_2;

    BOOST_CHECK_EQUAL_COLLECTIONS(b1_2, e1_2, b2_2, e2_2);
}


BOOST_AUTO_TEST_CASE(write_exception_Test)
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
	arrangementSet.read("../tests/domain/data/test.xdom");

	BOOST_CHECK_THROW(arrangementSet.write("test.txt", "not_a_format"), BioSeqDataLib::FormatException);

	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > seqSet;
	BOOST_CHECK_THROW(arrangementSet.write("test.txt", "not_a_format", seqSet), BioSeqDataLib::FormatException);

}


/**********************************************************
 *               OTHER FUNCTION TESTS                      *
 **********************************************************/





BOOST_AUTO_TEST_CASE( DomainArrangement_several_Test )
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
	arrangementSet.read("../tests/domain/data/BB20012.pfamScan");

	BOOST_CHECK_EQUAL(arrangementSet.size(), 27);

	std::vector<std::string> accessionIds = {"PF00009", "PF03144", "PF03143"};
	std::vector<std::string> geneIds = arrangementSet.search(accessionIds);
	BOOST_CHECK_EQUAL(geneIds.size(), 6);
	BOOST_CHECK_EQUAL(geneIds[0], "EFT1_PASMU");
	BOOST_CHECK_EQUAL(geneIds[1], "EFTU_BRELN");
	BOOST_CHECK_EQUAL(geneIds[2], "EFTU_LACLA");
	BOOST_CHECK_EQUAL(geneIds[3], "EFTU_LACPL");
	BOOST_CHECK_EQUAL(geneIds[4], "EFTU_SPIPL");
	BOOST_CHECK_EQUAL(geneIds[5], "EFTU_YEAST");


	std::vector<std::string> accessionIds2 = {"PF09106","PF09106"};
	geneIds = arrangementSet.contain(accessionIds2);
	BOOST_CHECK_EQUAL(geneIds.size(), 3);
	BOOST_CHECK_EQUAL(geneIds[0], "SELB_DESBA");
	BOOST_CHECK_EQUAL(geneIds[1], "SELB_ECOLI");
	BOOST_CHECK_EQUAL(geneIds[2], "SELB_MOOTH");


	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::DomainExt> arrangementSetDomain;
	arrangementSetDomain.read("../tests/domain/data/BB20012.pfamScan");
	BOOST_CHECK_EQUAL(arrangementSetDomain.size(), 27);
	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::DomainExt> &setDomain = arrangementSetDomain["IF2G_ARCFU"];

	BOOST_CHECK_EQUAL(setDomain.size(), 3);
	const BioSeqDataLib::DomainExt &dom1 = setDomain[1];
	BOOST_CHECK_EQUAL(dom1.name(), "GTP_EFTU_D2");
	BOOST_CHECK_EQUAL(dom1.accession(), "PF03144");
	BOOST_CHECK_EQUAL(dom1.start(), 248);
	BOOST_CHECK_EQUAL(dom1.end(), 326);
	BOOST_CHECK_EQUAL(dom1.envStart(), 248);
	BOOST_CHECK_EQUAL(dom1.envEnd(), 327);
	BOOST_CHECK_EQUAL(dom1.hmmStart(), 0);
	BOOST_CHECK_EQUAL(dom1.hmmEnd(), 72);
	BOOST_CHECK_EQUAL(dom1.hmmLength(), 74);
	BOOST_CHECK_CLOSE(dom1.evalue(), 0.0000000073, 0.000001);
	BOOST_CHECK_CLOSE(dom1.bitscore(), 35.6, 0.01);

}


 /*
  * Check that the begin and end functions of are working properly.
  * Iterators themsleves do not need to be tested as they are simply hidden STL containers.
  */
 BOOST_AUTO_TEST_CASE( Iterator_Test )
 {
		BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
		arrangementSet.read("../tests/domain/data/BB20012.pfamScan");

		// iterator and reverse_iterator test
		auto it = arrangementSet.begin();
		BOOST_CHECK_EQUAL(it->first, "1ega_A");
		it = --arrangementSet.end();
		BOOST_CHECK_EQUAL(it->first, "SELB_MOOTH");
		auto rit = arrangementSet.rbegin();
		BOOST_CHECK_EQUAL(rit->first, "SELB_MOOTH");
		rit = --arrangementSet.rend();
		BOOST_CHECK_EQUAL(rit->first, "1ega_A");

		// const_iterator and const_reverse_iterator test
		const auto &c_arrangementSet = arrangementSet;
		auto c_it = c_arrangementSet.begin();
		BOOST_CHECK_EQUAL(c_it->first, "1ega_A");
		c_it = --c_arrangementSet.end();
		BOOST_CHECK_EQUAL(c_it->first, "SELB_MOOTH");
		auto c_rit = c_arrangementSet.rbegin();
		BOOST_CHECK_EQUAL(c_rit->first, "SELB_MOOTH");
		c_rit = --c_arrangementSet.rend();
		BOOST_CHECK_EQUAL(c_rit->first, "1ega_A");
 }


 BOOST_AUTO_TEST_CASE( Find_Test )
 {
		BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
		arrangementSet.read("../tests/domain/data/BB20012.pfamScan");

		auto it = arrangementSet.find("SELB_MOOTH");
		BOOST_CHECK_EQUAL(it->first, "SELB_MOOTH");

		const auto &c_aSet = arrangementSet;
		auto c_it = c_aSet.find("SELB_MOOTH");
		BOOST_CHECK_EQUAL(c_it->first, "SELB_MOOTH");
 }


// Testing the count functions
BOOST_AUTO_TEST_CASE( Domain_count_Test )
{
	BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
	arrangementSet.read("../tests/domain/data/BB20012.pfamScan");

	// checking DomainArrangement counts
	const BioSeqDataLib::DomainArrangement<BioSeqDataLib::PfamDomain> &set = arrangementSet["IF2G_ARCFU"];
	BOOST_CHECK_EQUAL(set.size(), 3);
	std::map<std::string, int> daCounts = domainCounts(arrangementSet["1mky_A"]);
	BOOST_CHECK_EQUAL(2, daCounts.size());
	BOOST_CHECK_EQUAL(daCounts["PF01926"], 2);
	BOOST_CHECK_EQUAL(daCounts["PF14714"], 1);

	daCounts = clanCounts(arrangementSet["1mky_A"]);
	BOOST_CHECK_EQUAL(2, daCounts.size());
	BOOST_CHECK_EQUAL(daCounts["CL0023"], 2);
	BOOST_CHECK_EQUAL(daCounts["PF14714"], 1);

	// checking DomainArrangementSet counts
	std::map<std::string, int> setCounts = domainCounts(arrangementSet);
	BOOST_CHECK_EQUAL(12, setCounts.size());
	BOOST_CHECK_EQUAL(setCounts["PF09173"], 14);
	setCounts = clanCounts(arrangementSet);

	BOOST_CHECK_EQUAL(9, setCounts.size());
	BOOST_CHECK_EQUAL(setCounts["CL0023"], 52);
	BOOST_CHECK_EQUAL(setCounts["PF09106"], 3);
	setCounts = clanCounts(arrangementSet);

}

BOOST_AUTO_TEST_CASE( format_string_Test )
{
	BOOST_CHECK_EQUAL(BioSeqDataLib::getFormatString(BioSeqDataLib::DomainFileFormat::unknown), "unknown");
	BOOST_CHECK_EQUAL(BioSeqDataLib::getFormatString(BioSeqDataLib::DomainFileFormat::pfam), "Pfam");
	BOOST_CHECK_EQUAL(BioSeqDataLib::getFormatString(BioSeqDataLib::DomainFileFormat::hmmscan_domtbl), "HMMscan domtblout");
	BOOST_CHECK_EQUAL(BioSeqDataLib::getFormatString(BioSeqDataLib::DomainFileFormat::xdom), "XDOM");
	BOOST_CHECK_EQUAL(BioSeqDataLib::getFormatString(BioSeqDataLib::DomainFileFormat::ass), "SUPERFAMILY");
	BOOST_CHECK_EQUAL(BioSeqDataLib::getFormatString(BioSeqDataLib::DomainFileFormat::interpro_tsv), "InterPro TSV");
	BOOST_CHECK_EQUAL(BioSeqDataLib::getFormatString(BioSeqDataLib::DomainFileFormat::dama), "DAMA");
	BOOST_CHECK_EQUAL(BioSeqDataLib::getFormatString(BioSeqDataLib::DomainFileFormat::radiant), "RADIANT");
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* DOMAIN_ARRANGEMENT_SET_TEST_HPP_ */
