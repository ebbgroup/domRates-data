/*
 * DomainArrangement.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: ckeme_01
 */


#include "DomainArrangement.hpp"

namespace BioSeqDataLib
{

std::map<std::string, int>
clanCounts(const DomainArrangement<PfamDomain> &daSet, bool useDomain)
{
	std::map<std::string, int> summary;
	for (const PfamDomain &domain : daSet)
	{
		if (domain.clan().empty())
		{
			if (useDomain)
				summary[domain.accession()] += 1;
		}
		else
			summary[domain.clan()] += 1;
	}
	return summary;
}




}
