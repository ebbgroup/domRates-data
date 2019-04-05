/*
 * DomainArrangementSet.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: ckeme_01
 */


#include "DomainArrangementSet.hpp"

namespace BioSeqDataLib
{


std::string
getFormatString(const DomainFileFormat &format)
{
	switch (format)
	{
		case pfam:
			return "Pfam";
			break;
		case hmmscan_domtbl:
			return "HMMscan domtblout";
			break;
		case xdom:
			return "XDOM";
			break;
		case ass:
			return "SUPERFAMILY";
			break;
		case interpro_tsv:
			return "InterPro TSV";
			break;
		case dama:
			return "DAMA";
			break;
		case radiant:
			return "RADIANT";
			break;
		default:
			return "unknown";
	}
}

//<seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>
template<>
void
DomainArrangementSet<Domain>::pfamTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<Domain> &da)
{
	da.emplace_back(acc, stoul(tokens[1])-1, stoul(tokens[2])-1, stod(tokens[12]));
}

template<>
void
DomainArrangementSet<DomainExt>::pfamTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<DomainExt> &da)
{
	da.emplace_back(acc, tokens[6], std::stoul(tokens[1])-1, std::stoul(tokens[2])-1, std::stoul(tokens[3])-1, std::stoul(tokens[4])-1, std::stoul(tokens[8])-1, std::stoul(tokens[9])-1, std::stoul(tokens[10]), std::stod(tokens[11]), std::stod(tokens[12]));
}


template<>
void
DomainArrangementSet<PfamDomain>::pfamTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<PfamDomain> &da)
{
	//accession, std::string name, size_t seqStart, size_t seqEnd, size_t envStart, size_t envEnd, size_t hmmStart, size_t hmmEnd, size_t hmm_length, double bit_score, double evalue, double significance, std::string clan, std::string type
	std::string clan = (tokens[14] == "No_clan") ? "" : tokens[14];
	da.emplace_back(acc, tokens[6], std::stoul(tokens[1])-1, std::stoul(tokens[2])-1, std::stoul(tokens[3])-1, std::stoul(tokens[4])-1, std::stoul(tokens[8])-1, std::stoul(tokens[9])-1, std::stoul(tokens[10]), std::stod(tokens[11]), std::stod(tokens[12]), std::stod(tokens[13]), clan, tokens[7]);
}


template<>
void
DomainArrangementSet<Domain>::radiantTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<Domain> &da)
{
	da.emplace_back(acc, stoul(tokens[1])-1, stoul(tokens[2])-1, -1);
}

/*
template<>
void
DomainArrangementSet<DomainExt>::radiantTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<DomainExt> &da)
{
	da.emplace_back(acc, tokens[6], std::stoul(tokens[1])-1, std::stoul(tokens[2])-1, std::stoul(tokens[1])-1, std::stoul(tokens[2])-1, std::stoul(tokens[8])-1, std::stoul(tokens[9])-1, std::stoul(tokens[10]), std::stod(tokens[11]), std::stod(tokens[12]));
}


template<>
void
DomainArrangementSet<PfamDomain>::radiantTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<PfamDomain> &da)
{

	//accession, std::string name, size_t seqStart, size_t seqEnd, size_t envStart, size_t envEnd, size_t hmmStart, size_t hmmEnd, size_t hmm_length, double bit_score, double evalue, double significance, std::string clan, std::string type
	std::string clan = (tokens[14] == "No_clan") ? "" : tokens[14];
	da.emplace_back(acc, tokens[6], std::stoul(tokens[1])-1, std::stoul(tokens[2])-1, std::stoul(tokens[3])-1, std::stoul(tokens[4])-1, std::stoul(tokens[8])-1, std::stoul(tokens[9])-1, std::stoul(tokens[10]), std::stod(tokens[11]), std::stod(tokens[12]), std::stod(tokens[13]), clan, tokens[7]);
}*/


/*
 * hmmscan tblout format
 * #                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
 * # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
 * # ------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
 * # GTP_EFTU             PF00009.22   188 1kjz_A               -            400   5.9e-43  146.4   0.1   1   1   1.8e-45   1.2e-42  145.5   0.1     3   187     5   193     3   194 0.95 Elongation factor Tu GTP binding domain
 */

template<>
void
DomainArrangementSet<Domain>::hmmTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<Domain> &da)
{
	da.emplace_back(acc, stoul(tokens[17])-1, stoul(tokens[18])-1, stod(tokens[12]));
}

template<>
void
DomainArrangementSet<DomainExt>::hmmTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<DomainExt> &da)
{
	//accession, name seqStart seqEnd envStart envEnd hmmStart hmmEnd hmm_length bit_score evalue)
	da.emplace_back(acc, tokens[0], std::stoul(tokens[17])-1, std::stoul(tokens[18])-1, std::stoul(tokens[19])-1, std::stoul(tokens[20])-1, std::stoul(tokens[15])-1, std::stoul(tokens[16])-1, std::stoul(tokens[2]), std::stod(tokens[13]), std::stod(tokens[12]));
}

template<>
void
DomainArrangementSet<PfamDomain>::hmmTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<PfamDomain> &da)
{
	//std::string accession, std::string name, size_t seqStart, size_t seqEnd, size_t envStart, size_t envEnd, size_t hmmStart, size_t hmmEnd, size_t hmm_length, double bit_score, double evalue, double significance, std::string clan, std::string type
	da.emplace_back(acc, tokens[0], std::stoul(tokens[17])-1, std::stoul(tokens[18])-1, std::stoul(tokens[19])-1, std::stoul(tokens[20])-1, std::stoul(tokens[15])-1, std::stoul(tokens[16])-1, std::stoul(tokens[2]), std::stod(tokens[13]), std::stod(tokens[12]), 0, "", "");
}


template<>
void
DomainArrangementSet<SFDomain>::hmmTokens2Domain_(std::string &&acc, const std::vector<std::string> &tokens, DomainArrangement<SFDomain> &da)
{
	da.emplace_back(acc, tokens[0], std::stoul(tokens[17])-1, std::stoul(tokens[18])-1, std::stoul(tokens[19])-1, std::stoul(tokens[20])-1, std::stoul(tokens[15])-1, std::stoul(tokens[16])-1, std::stoul(tokens[2]), std::stod(tokens[13]), std::stod(tokens[12]), "");
}


template<>
void
DomainArrangementSet<Domain>::damaTokens2Domain_(const std::vector<std::string> &tokens, DomainArrangement<Domain> &da)
{
	da.emplace_back(tokens[4], stoul(tokens[1])-1, stoul(tokens[2])-1, stod(tokens[0]));
}

template<>
void
DomainArrangementSet<Domain>::assTokens2Domain_(const std::vector<std::string> &tokens, DomainArrangement<Domain> &da)
{
	auto positions = split(tokens[2], "-");
	da.emplace_back(std::move(tokens[1]), stoul(positions[0])-1, stoul(positions[1])-1, stod(tokens[6]));
}

template<>
void
DomainArrangementSet<DomainExt>::assTokens2Domain_(const std::vector<std::string> &tokens, DomainArrangement<DomainExt> &da)
{
	auto positions = split(tokens[2], "-");
	da.emplace_back(tokens[8], tokens[1], std::stoul(positions[0])-1,	std::stoul(positions[1])-1, std::stoul(positions[0])-1, std::stoul(positions[1])-1, std::stoul(tokens[4])-1, 0, 0, -1, std::stod(tokens[3]));
}

template<>
void
DomainArrangementSet<SFDomain>::assTokens2Domain_(const std::vector<std::string> &tokens, DomainArrangement<SFDomain> &da)
{
	auto positions = split(tokens[2], "-");
	da.emplace_back(tokens[8], tokens[1], std::stoul(positions[0])-1,	std::stoul(positions[1])-1, std::stoul(positions[0])-1, std::stoul(positions[1])-1, std::stoul(tokens[4])-1, 0, 0, -1, std::stod(tokens[3]), tokens[7]);
}



template<>
void
DomainArrangementSet<Domain>::interProTSVTokens2Domain_(const std::vector<std::string> &tokens, DomainArrangement<Domain> &da)
{
	//EFX70648	2c281c89dcf49d2e45ab9cf625a5ac16	1040	PRINTS	PR00419	Adrenodoxin reductase family signature	337	351	1.299999245E-14	T	24-02-2015
	std::map<std::string, DomainDB> str2db = {{"Gene3D", DomainDB::gene3d}, {"Hamap", DomainDB::hamap}, {"PANTHER", DomainDB::panther}, {"Pfam", DomainDB::pfam}, {"PIRSF", DomainDB::pirsf}, {"PRINTS", DomainDB::prints}, {"ProDom", DomainDB::prodom}, {"ProSitePatterns", DomainDB::prosite}, {"ProSiteProfiles", DomainDB::prosite}, {"SMART", DomainDB::smart}, {"SUPERFAMILY", DomainDB::superfamily}, {"TIGRFAM", DomainDB::tigrfams}};
	da.emplace_back(tokens[4], std::stoul(tokens[6]), std::stoul(tokens[7]), (tokens[8][0] == '-') ? 999 : stod(tokens[8]), str2db[tokens[3]]);
}



template<>
void
DomainArrangementSet<PfamDomain>::_writePfamScanOutput(std::ofstream &outFile)
{
	outFile << "# pfam_scan.pl\n";
	outFile << "# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>\n\n";
	size_t len = 0;
	for (const auto &pair : arrangements_)
	{
		if (len < pair.first.length())
			len = pair.first.length();
	}

	outFile.precision(1);
	for (const auto &pair : arrangements_)
	{
		for (const auto &domain : pair.second)
		{
			outFile.width(len);
			outFile << std::left << pair.first;
			outFile.width(7);
			outFile << std::right << domain.start()+1;
			outFile.width(7);
			outFile << domain.end()+1;
			outFile.width(7);
			outFile << domain.envStart()+1;
			outFile.width(7);
			outFile << domain.envEnd()+1;
			outFile << " " << domain.accession() << " ";// << "." << std::setw(4) << std::left << domain.version();
			outFile << std::setw(18) << std::left << domain.name() <<  domain.type();
			outFile << std::setw(6) << std::right << domain.hmmStart()+1 << std::setw(6) << std::right << domain.hmmEnd()+1 << std::setw(6) << std::right << domain.hmmLength();
			outFile << std::setw(9) << std::fixed << std::right << domain.bitscore();
			outFile << std::setw(10) << std::scientific << std::right << domain.evalue();
			outFile.unsetf ( std::ios::floatfield);
			outFile << std::setw(4) << std::right << domain.significance();
			if (domain.clan().empty())
				outFile << " No_Clan";
			else
				outFile << " " << domain.clan();
			outFile << "\n";
		}
	}
}








/**
 *
 * \brief Counts the occurrence of each clan in a domainArrangmentSet.
 * \param If useDomain is true, the domain accession number is used if the domain is not assigned to any clan.
 * @return The domain and the number of their occurrence.
 * \relates DomainArrangement
 */
std::map<std::string, int>
clanCounts(const DomainArrangementSet<PfamDomain> &daSet, bool useDomain)
{
	std::map<std::string, int> summary;
	for (const auto da : daSet)
	{
		for (const PfamDomain &domain : da.second)
		{
			if (domain.clan().empty())
			{
				if (useDomain)
					summary[domain.accession()] += 1;
			}
			else
				summary[domain.clan()] += 1;
		}
	}
	return summary;
}

}
