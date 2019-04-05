
#include "domainProgs.hpp"

using namespace std;

namespace BioSeqDataLib
{

    void
    runPfamScan(const fs::path &inputFile, const fs::path &outputFile, const fs::path &database, size_t nThreads,
        const fs::path &logFile)
    {
        std::string log = logFile.empty() ? "/dev/null" : logFile.string();

        string cmd = "pfam_scan.pl -fasta " + inputFile.string() + " -dir " + database.string() + " -cpu "
            + std::to_string(nThreads) + " -outfile " + outputFile.string() + " 2>> " + log;
        if (system(cmd.c_str()))
            throw std::runtime_error("pfam_scan.pl failed to run on " + inputFile.string());
    }

    void
    runPfamScan(const fs::path &inputFile, const fs::path &outputFile, const fs::path &database, float evalue, size_t nThreads,
        const fs::path &logFile)
    {
        std::string log = logFile.empty() ? "/dev/null" : logFile.string();

        string cmd = "pfam_scan.pl -fasta " + inputFile.string() + " -dir " + database.string() + " -cpu "
            + std::to_string(nThreads) + " -outfile " + outputFile.string() + " -e_seq " + to_string(evalue) + " -e_dom "+  to_string(evalue) + " 2>> " + log;
        if (system(cmd.c_str()))
            throw std::runtime_error("pfam_scan.pl failed to run on " + inputFile.string());
    }


    // join <(cut -f 1,2,8 pfamA.txt) <(cut -f 1,2 Pfam-A.clans.tsv) > pfamA_31.txt
    void
    readDomainInfo(const fs::path &pfamInfo, std::map<std::string, DomainInfo> &info)
    {
    	AlgorithmPack::Input in(pfamInfo);
    	std::string line;
    	while (getline(in, line))
    	{
    		auto tokens = split(line, " ");
    		if (tokens.size()==3)
    			info.emplace(std::piecewise_construct, std::make_tuple(tokens[0]), std::make_tuple(tokens[1], tokens[2], "No_clan"));
    		else
    			info.emplace(std::piecewise_construct, std::make_tuple(tokens[0]), std::make_tuple(tokens[1], tokens[2], tokens[3]));
    	}
    }


}
