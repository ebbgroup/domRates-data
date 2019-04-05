/*
 * SequenceSet_Test.cpp
 *
 *  Created on: 20 Oct 2013
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
#include "../../src/sequence/SeqSetFunctions.hpp"

BOOST_AUTO_TEST_SUITE(SequenceSet_Test)


BOOST_AUTO_TEST_CASE( SequenceSet_simple_Test)
{
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set;
	set.emplace_back("seq1", "ACGTTTT", "", "");
	BOOST_CHECK_EQUAL(set[0].name(), "seq1");
	BOOST_CHECK_EQUAL(set[0].seq(), "ACGTTTT");
}



BOOST_AUTO_TEST_CASE( Sequence_read_Test)
{
	// test fasta reading
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set("../tests/sequence/data/seqSet.fasta", std::vector<std::string>(1, "gnl|em|HSFAU"));
	BOOST_CHECK_EQUAL(set[0].name(), "gnl|em|HSFAU");
	BOOST_CHECK_EQUAL(set[0].comment(), "(X65923) H.sapiens fau mRNA");
	BOOST_CHECK_EQUAL(set.size(), 1);

	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > setb("../tests/sequence/data/seqSet.fasta", std::vector<std::string>(1, "gnl|em|HSFAU"), true);
	BOOST_CHECK_EQUAL(setb[0].name(), "gnl|em|HSFAU1");
	BOOST_CHECK_EQUAL(setb[0].comment(), "(X65921) H.sapiens fau 1 gene");
	BOOST_CHECK_EQUAL(setb.size(), 1);

	// test genbank reading
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set1("../tests/sequence/data/seqSet.genbank", std::vector<std::string>(1, "HSFAU"));
	BOOST_CHECK_EQUAL(set1[0].name(), "HSFAU");
	BOOST_CHECK_EQUAL(set1[0].accession(), "X65923");
	BOOST_CHECK_EQUAL(set1[0].comment(), "H.sapiens fau mRNA");
	BOOST_CHECK_EQUAL(set1[0].seq(), "ttcctctttctcgactccatcttcgcggtagctgggaccgccgttcagtcgccaatatgcagctctttgtccgcgcccaggagctacacaccttcgaggtgaccggccaggaaacggtcgcccagatcaaggctcatgtagcctcactggagggcattgccccggaagatcaagtcgtgctcctggcaggcgcgcccctggaggatgaggccactctgggccagtgcggggtggaggccctgactaccctggaagtagcaggccgcatgcttggaggtaaagttcatggttccctggcccgtgctggaaaagtgagaggtcagactcctaaggtggccaaacaggagaagaagaagaagaagacaggtcgggctaagcggcggatgcagtacaaccggcgctttgtcaacgttgtgcccacctttggcaagaagaagggccccaatgccaactcttaagtcttttgtaattctggctttctctaataaaaaagccacttagttcagtcaaaaaaaaaa");
	BOOST_CHECK_EQUAL(set1.size(), 1);

	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set1b("../tests/sequence/data/seqSet.genbank", std::vector<std::string>(1, "HSFAU"), true);
	BOOST_CHECK_EQUAL(set1b[0].name(), "HSFAU1");
	BOOST_CHECK_EQUAL(set1b[0].accession(), "X65921");
	BOOST_CHECK_EQUAL(set1b[0].comment(), "H.sapiens fau 1 gene");
	BOOST_CHECK_EQUAL(set1b[0].seq(), "ctaccattttccctctcgattctatatgtacactcgggacaagttctcctgatcgaaaacggcaaaactaaggccccaagtaggaatgccttagttttcggggttaacaatgattaacactgagcctcacacccacgcgatgccctcagctcctcgctcagcgctctcaccaacagccgtagcccgcagccccgctggacaccggttctccatccccgcagcgtagcccggaacatggtagctgccatctttacctgctacgccagccttctgtgcgcgcaactgtctggtcccgccccgtcctgcgcgagctgctgcccaggcaggttcgccggtgcgagcgtaaaggggcggagctaggactgccttgggcggtacaaatagcagggaaccgcgcggtcgctcagcagtgacgtgacacgcagcccacggtctgtactgacgcgccctcgcttcttcctctttctcgactccatcttcgcggtagctgggaccgccgttcaggtaagaatggggccttggctggatccgaagggcttgtagcaggttggctgcggggtcagaaggcgcggggggaaccgaagaacggggcctgctccg");
	BOOST_CHECK_EQUAL(set1b.size(), 1);

	// test swissprot reading
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set2("../tests/sequence/data/seqSet.swissprot", std::vector<std::string>(1, "104K_THEPA"));
	BOOST_CHECK_EQUAL(set2[0].name(), "104K_THEPA");
	BOOST_CHECK_EQUAL(set2[0].accession(),"P15711;");
	BOOST_CHECK_EQUAL(set2[0].comment(), "104 kDa microneme-rhoptry antigen.");
	BOOST_CHECK_EQUAL(set2[0].seq(), "MKFLILLFNILCLFPVLAADNHGVGPQGASGVDPITFDINSNQTGPAFLTAVEMAGVKYLQVQHGSNVNIHRLVEGNVVIWENASTPLYTGAIVTNNDGPYMAYVEVLGDPNLQFFIKSGDAWVTLSEHEYLAKLQEIRQAVHIESVFSLNMAFQLENNKYEVETHAKNGANMVTFIPRNGHICKMVYHKNVRIYKATGNDTVTSVVGFFRGLRLLLINVFSIDDNGMMSNRYFQHVDDKYVPISQKNYETGIVKLKDYKHAYHPVDLDIKDIDYTMFHLADATYHEPCFKIIPNTGFCITKLFDGDQVLYESFNPLIHCINEVHIYDRNNGSIICLHLNYSPPSYKAYLVLKDTGWEATTHPLLEEKIEELQDQRACELDVNFISDKDLYVAALTNADLNYTMVTPRPHRDVIRVSDGSEVLWYYEGLDNFLVCAWIYVSDGVASLVHLRIKDRIPANNDIYVLKGDLYWTRITKIQFTQEIKRLVKKSKKKLAPITEEDSDKHDEPPEGPGASGLPPKAPGDKEGSEGHKGPSKGSDSSKEGKKPGSGKKPGPAREHKPSKIPTLSKKPSGPKDPKHPRDPKEPRKSKSPRTASPTRRPSPKLPQLSKLPKSTSPRSPPPPTRPSSPERPEGTKIIKTSKPPSPKPPFDPSFKEKFYDDYSKAASRSKETKTTVVLDESFESILKETLPETPGTPFTTPRPVPPKRPRTPESPFEPPKDPDSPSTSPSEFFTPPESKRTRFHETPADTPLPDVTAELFKEPDVTAETKSPDEAMKRPRSPSEYEDTSPGDYPSLPMKRHRLERLRLTTTEMETDPGRMAKDASGKPVKLKRSKSFDDLTTVELAPEPKASRIVVDDEGTEADDEETHPPEERQKTEVRRRRPPKKPSKSPRPSKPKKPKKPDSAYIPSILAILVVSLIVGIL");
	BOOST_CHECK_EQUAL(set2.size(), 1);

	// test stockholm reading
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set3("../tests/sequence/data/seqSet.stockholm", std::vector<std::string>(1, "AP001509.1"));
	BOOST_CHECK_EQUAL(set3[0].name(), "AP001509.1");
	BOOST_CHECK_EQUAL(set3[0].comment(), "");
	BOOST_CHECK_EQUAL(set3[0].seq(), "UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCAAUAUGG-GAUGAGGGUCUCUAC-AGGUA-CCGUAAA-UACCUAGCUACGAAAAGAAUGCAGUUAAUGU");
	BOOST_CHECK_EQUAL(set3.size(), 1);

	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> >::iterator it=set3.begin();
	BOOST_CHECK_EQUAL(set3[0], *it);

	// test msf reading
	std::vector<std::string> seq_names(1, "CD5R_HUMAN");
	seq_names.push_back("CD5R_MOUSE");
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set4("../tests/sequence/data/seqSet.msf", seq_names);
	BOOST_CHECK_EQUAL(set4.size(), 2);
	BOOST_CHECK_EQUAL(set4[0].name(), "CD5R_HUMAN");
	BOOST_CHECK_EQUAL(set4[0].comment(), "");
	BOOST_CHECK_EQUAL(set4[0].seq(), "MGTVLSLSPSYRKATLFEDGAATVGHYTAVQNSKNAKDKNLKRHSIISVLPWKRIVAVSAKKKNSKKVQPNSSYQNNITHLNNENLKKSLSCANLSTFAQPPPAQPPAPPASQLSGSQTGGSSSVKKAPHPAVTSAGTPKRVIVQASTSELLRCLGEFLCRRCYRLKHLSPTDPVLWLRSVDRSLLLQGWQDQGFITPANVVFLYMLCRDVISSEVGSDHELQAVLLTCLYLSYSYMGNEISYPLKPFLVESCKEAFWDRCLSVINLMSSKMLQINADPHYFTQVFSDLKNESGQEDKKRLLLGLDR");

	// test phylip reading
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set5("../tests/sequence/data/seqSet.phylip_s");
	BOOST_CHECK_EQUAL(set5.size(), 3);
	BOOST_CHECK_EQUAL(set5[1].name(), "ALEU_HORVU");
	BOOST_CHECK_EQUAL(set5[2].seq(), "------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK------------FHFKSWMSKHRKTY-STEEYHHRLQTFASNWRKINAHN----NGNHTFKMALNQFSDMSFAEIKHKYLWSEPQNCSAT--KSNYLRGT--GPYPPSVDWRKKGNFVSPVKNQGACGSCWTFSTTGALESAIAIATGKMLSLAEQQLVDCAQDFNNY--------GCQGGLPSQAFEYILYNKGIMGEDTYPYQGKDGY-CKFQPGKAIGFVKDVANITIYDEEAMVEAVALYNPVSFAFEVTQDFMMYRTGIYSSTSCHKTPDKVNHAVLAVGYGEKNGI-----PYWIVKNSWGPQWGMNGYFLIERGKNMCGLAACASYPIPLV");

	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set6("../tests/sequence/data/seqSet.phylip_s", std::vector<std::string>(1, "ALEU_HORVU"));
	BOOST_CHECK_EQUAL(set6.size(), 1);
	BOOST_CHECK_EQUAL(set6[0].seq(), "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALRFARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTN----RKGLPYRLGINRFSDMSWEEFQATRL-GAAQTCSATLAGNHLMRDA--AALPETKDWREDG-IVSPVKNQAHCGSCWTFSTTGALEAAYTQATGKNISLSEQQLVDCAGGFNNF--------GCNGGLPSQAFEYIKYNGGIDTEESYPYKGVNGV-CHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVIDGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGV-----PYWLIKNSWGADWGDNGYFKMEMGKNMCAIATCASYPVVAA");

	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set7("../tests/sequence/data/seqSet.phylip_i");
	BOOST_CHECK_EQUAL(set7.size(), 3);
	BOOST_CHECK_EQUAL(set7[1].name(), "ALEU_HORVU");
	BOOST_CHECK_EQUAL(set7[2].seq(), "------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK------------FHFKSWMSKHRKTY-STEEYHHRLQTFASNWRKINAHN----NGNHTFKMALNQFSDMSFAEIKHKYLWSEPQNCSAT--KSNYLRGT--GPYPPSVDWRKKGNFVSPVKNQGACGSCWTFSTTGALESAIAIATGKMLSLAEQQLVDCAQDFNNY--------GCQGGLPSQAFEYILYNKGIMGEDTYPYQGKDGY-CKFQPGKAIGFVKDVANITIYDEEAMVEAVALYNPVSFAFEVTQDFMMYRTGIYSSTSCHKTPDKVNHAVLAVGYGEKNGI-----PYWIVKNSWGPQWGMNGYFLIERGKNMCGLAACASYPIPLV");

	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set8("../tests/sequence/data/seqSet.phylip_i", std::vector<std::string>(1, "ALEU_HORVU"));
	BOOST_CHECK_EQUAL(set8.size(), 1);
	BOOST_CHECK_EQUAL(set8[0].seq(), "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALRFARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTN----RKGLPYRLGINRFSDMSWEEFQATRL-GAAQTCSATLAGNHLMRDA--AALPETKDWREDG-IVSPVKNQAHCGSCWTFSTTGALEAAYTQATGKNISLSEQQLVDCAGGFNNF--------GCNGGLPSQAFEYIKYNGGIDTEESYPYKGVNGV-CHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVIDGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGV-----PYWLIKNSWGADWGDNGYFKMEMGKNMCAIATCASYPVVAA");

	// test amps
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set9("../tests/sequence/data/seqSet.amps");
	BOOST_CHECK_EQUAL(set9.size(), 6);
	BOOST_CHECK_EQUAL(set9[1].name(), "A0231");
	BOOST_CHECK_EQUAL(set9[2].seq(), "NAP--LV");

	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set10("../tests/sequence/data/seqSet.amps", std::vector<std::string>(1, "Four_Alpha"));
	BOOST_CHECK_EQUAL(set10.size(), 1);
	BOOST_CHECK_EQUAL(set10[0].seq(), "DAPPWLV");




}

BOOST_AUTO_TEST_CASE( Sequence_write_Test)
{
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence <>> set1;
	set1.push_back(BioSeqDataLib::Sequence<>("seq1", "ACGT", "", ""));
	set1.push_back(BioSeqDataLib::Sequence<>("seq2", "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALRFARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTN", "", ""));
	set1.write("test.fa", "fasta");

	std::ifstream ifs1("../tests/sequence/data/fasta_out.fa");
	std::ifstream ifs2("test.fa");
	std::istream_iterator<char> b1(ifs1), e1;
	std::istream_iterator<char> b2(ifs2), e2;
	BOOST_CHECK_EQUAL_COLLECTIONS(b1, e1, b2, e2);
}



BOOST_AUTO_TEST_CASE( Sequence_changing_Test)
{
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence <>> set1, set2;
	set1.push_back(BioSeqDataLib::Sequence<>("seq1", "ACGT", "", ""));
	BioSeqDataLib::Sequence<> seq("seq2", "ACGT", "", "");
	set1.push_back(seq);
	BOOST_CHECK_EQUAL(set1.size(), 2);

	set2.transfer(set1);
	BOOST_CHECK_EQUAL(set2.size(), 2);
	BOOST_CHECK_EQUAL(set1.size(), 0);
	BOOST_CHECK_EQUAL(set1.empty(), true);
	BOOST_CHECK_EQUAL(set2[0].name(), "seq1");
	BOOST_CHECK_EQUAL(set2[1].name(), "seq2");

	BOOST_CHECK_EQUAL(set2["seq2"].name(), "seq2");

	BOOST_CHECK_EQUAL(set2.find("seq2")->name(), "seq2");

	BioSeqDataLib::Sequence<> tmpSeq("seq3", "ACGT", "", "");
	set1.push_back(tmpSeq);
	BOOST_CHECK_EQUAL(set1.size(), 1);

	set2.erase(set2.find("seq2"));
	BOOST_CHECK_EQUAL(set2.size(), 1);

}

BOOST_AUTO_TEST_CASE( Sequence_analysis_Test)
{
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence <>> set;
	set.push_back(BioSeqDataLib::Sequence<>("seq1", "ACGT", "", ""));
	set.push_back(BioSeqDataLib::Sequence<>("seq2", "gCG", "", ""));
	set.push_back(BioSeqDataLib::Sequence<>("seq3", "aCGT", "", ""));

	BioSeqDataLib::SimilarityMatrix<float> simMat;
	simMat.setBlosum62();
	int identity = BioSeqDataLib::avgId(set, simMat)*1000;
	BOOST_CHECK_EQUAL(identity, 777);


	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence <>> set1;
	set1.push_back(BioSeqDataLib::Sequence<>("seq1", "ACGTa", "", ""));
	set1.push_back(BioSeqDataLib::Sequence<>("seq2", "gCG", "", ""));
	set1.push_back(BioSeqDataLib::Sequence<>("seq3", "aCGT", "", ""));
	int val = BioSeqDataLib::avgLength(set1);
	BOOST_CHECK_EQUAL(val, 4);

	BOOST_CHECK_EQUAL(BioSeqDataLib::maxLength(set1), 5);
	BOOST_CHECK_EQUAL(BioSeqDataLib::minLength(set1), 3);

}


BOOST_AUTO_TEST_CASE( SequenceSet_iterator_Test)
{
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> > set3("../tests/sequence/data/seqSet.stockholm");
	BOOST_CHECK_EQUAL(set3.size(), 2);

	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> >::iterator it=set3.begin();
	BOOST_CHECK_EQUAL(set3[0], *it);
	BOOST_CHECK_EQUAL(set3[0].name(), it->name());
	++it;
	BOOST_CHECK_EQUAL(set3[1], *it);
	BOOST_CHECK_EQUAL(set3[1].name(), it->name());

	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> >::const_iterator cit=set3.cend();
	--cit;
	BOOST_CHECK_EQUAL(set3[1], *cit);
	BOOST_CHECK_EQUAL(set3[1].name(), it->name());
	--cit;
	BOOST_CHECK_EQUAL(set3[0], *cit);
	BOOST_CHECK_EQUAL(set3[0].name(), cit->name());


	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence<> >::reverse_iterator itr, itr_end =set3.rend();
	size_t pos=2;
	for (itr=set3.rbegin(); itr!=itr_end; ++itr)
		BOOST_CHECK_EQUAL(set3[--pos], *itr);
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

BOOST_AUTO_TEST_CASE( SequenceSet_output_Test)
{
	BioSeqDataLib::SequenceSet<BioSeqDataLib::Sequence <>> set1;
	set1.push_back(BioSeqDataLib::Sequence<>("seq1", "ACGT", "", "comment"));
	set1.push_back(BioSeqDataLib::Sequence<>("seq2", "MAHAR", "", ""));

	boost::test_tools::output_test_stream output;
	{
		cout_redirect guard( output.rdbuf( ) );
		set1.write("fasta");
	}
	BOOST_CHECK( output.is_equal( ">seq1 comment\nACGT\n>seq2\nMAHAR\n") );
}




BOOST_AUTO_TEST_SUITE_END()
