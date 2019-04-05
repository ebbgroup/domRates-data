/*
 * runPhylo.cpp
 *
 *  Created on: Sep 15, 2015
 *      Author: ckeme_01
 */


#include "runPhylo.hpp"


using namespace std;


namespace BioSeqDataLib
{

/**
 * A class defining different supported alignment programs
 */


string
makeClustalOCmd(const fs::path &inF, const fs::path &outF, int nThreads, AlnFormat format, const fs::path &logFile)
{
	string outFormat;
	switch (format)
	{
		case AlnFormat::fasta:
			outFormat = "fasta";
			break;
		case AlnFormat::phylip:
			outFormat = "phylip";
			break;
		case AlnFormat::clustal:
			outFormat = "clustal";
	}
	string logF;
	if (logFile.empty())
		logF =  " > /dev/null";
	else
		logF = " >> " + logFile.string() +" 2>&1";
	return "clustalo -i "+ inF.string() + " --threads " + to_string(nThreads) + " --force --outfmt=" + outFormat + " -o " + outF.string() + logF;
}

void
runClustalO(const fs::path &inF, const fs::path &outF, int nThreads, const AlnFormat &format, const fs::path &log)
{
	string command = makeClustalOCmd(inF, outF, nThreads, format, log);
	if (system(command.c_str()))
		throw std::runtime_error("ClustalO failed to run on " + inF.string());
}


string
makeTCoffeeCmd(const fs::path &inF, const fs::path &outF, int nThreads, const AlnFormat &format, const fs::path &logFile)
{
	string outFormat;
	switch (format)
	{
		case AlnFormat::fasta:
			outFormat = "fasta_aln";
			break;
		case AlnFormat::phylip:
			outFormat = "phylip";
			break;
		case AlnFormat::clustal:
			outFormat = "clustalw";
	}
	string logF;
	if (logFile.empty())
		logF =  " > /dev/null";
	else
		logF = " >> " + logFile.string() +" 2>&1";
	return "t_coffee -seq "+ inF.string() + " -n_core " + to_string(nThreads) + " -output " + outFormat + " -outfile " + outF.string() + logF;
}


void
runTCoffee(const fs::path &inF, const fs::path &outF, int nThreads, const AlnFormat format, const fs::path &logFile)
{
	string command = makeTCoffeeCmd(inF, outF, nThreads, format, logFile);
	if (system(command.c_str()))
		throw std::runtime_error("T-Coffee failed to run on " + inF.string());
}



string
makeMafftCmd(const fs::path &inF, const fs::path &outF, int nThreads, const AlnFormat format, const fs::path &logFile)
{
	string outFormat;
	switch (format)
	{
		case AlnFormat::fasta:
			outFormat = " ";
			break;
		case AlnFormat::clustal:
			outFormat = " --clustalout ";
			break;
		default:
			throw (std::runtime_error("The choosen output format is no supported by Mafft"));
	}
	string logF;
	if (logFile.empty())
		logF =  "";
	else
		logF = " 2>> " + logFile.string();
	return "mafft --thread " + to_string(nThreads) + outFormat + inF.string() + " > " + outF.string() + logF;
}


void
runMafft(const fs::path &inF, const fs::path &outF, int nThreads, AlnFormat format, const fs::path &logFile)
{
	string command = makeMafftCmd(inF, outF, nThreads, format, logFile);
	if (system(command.c_str()))
		throw std::runtime_error("Mafft failed to run on " + inF.string());
}

string
makeMSAProbsCmd(const fs::path &inF, const fs::path &outF, int nThreads, const AlnFormat format, const fs::path &logFile)
{
	string outFormat;
	switch (format)
	{
		case AlnFormat::fasta:
			outFormat = " ";
			break;
		case AlnFormat::clustal:
			outFormat = " -clustalw ";
			break;
		default:
			throw (std::runtime_error("The choosen output format is no supported by MSAProbs"));
	}
	string logF;
	if (logFile.empty())
		logF =  " >> /dev/null";
	else
		logF = " >> " + logFile.string() +" 2>&1";
	return "msaprobs -num_threads " + to_string(nThreads) + " -o " + outF.string() + outFormat + inF.string() + logF;
}


void
runMSAProbs(const fs::path &inF, const fs::path &outF, int nThreads, AlnFormat format, const fs::path &logFile)
{
	string command = makeMSAProbsCmd(inF, outF, nThreads, format, logFile);
	if (system(command.c_str()))
		throw std::runtime_error("MSAProbs failed to run on " + inF.string());
}

void
runAlignment(const fs::path &inF, AlnProgram program, const fs::path &outF, int nThreads, AlnFormat format, const fs::path &logFile)
{
	switch (program)
	{
		case AlnProgram::clustalo:
			runClustalO(inF, outF, nThreads, format, logFile);
			break;
		case AlnProgram::mafft:
			runMafft(inF, outF, nThreads, format, logFile);
			break;
		case AlnProgram::tcoffee:
			runTCoffee(inF, outF, nThreads, format, logFile);
			break;
		case AlnProgram::msaprobs:
			runMSAProbs(inF, outF, nThreads, format, logFile);
			break;
		default:
			throw std::runtime_error("Unsupported alignment program");
	}
}


string
makeTrimAlCmd(const fs::path &inF, const fs::path &outF, const string &parameter, const fs::path &logFile)
{
	string logF;
	if (logFile.empty())
		logF =  " >> /dev/null";
	else
		logF = " >> " + logFile.string() +" 2>&1";
	return "trimal -in " +  inF.string() + " -out " + outF.string() + " " + parameter + logF;
}


void
runTrimAl(const fs::path &inF, const fs::path &outF, const string &parameter, const fs::path &logFile)
{
	string command = makeTrimAlCmd(inF, outF, parameter, logFile);
	if (system(command.c_str()))
		throw std::runtime_error("TrimAl failed to run on " + inF.string());
}


string
makeProttestCmd(const fs::path &inF, const fs::path &outF, int nThreads, const string &parameter, const fs::path &logFile)
{
	string logF;
	if (logFile.empty())
		logF =  " >> /dev/null";
	else
		logF = " >> " + logFile.string() + " 2>&1";
	return "prottest -i " + inF.string() +  " -threads " + to_string(nThreads) + " -o " + outF.string() + " " + parameter + logF;
}


std::string
runProttest(const fs::path &inF, const fs::path &outF, int nThreads, const string &parameter, const fs::path &logFile)
{
	string command = makeProttestCmd(inF, outF, nThreads, parameter, logFile);
	if (system(command.c_str()))
		throw std::runtime_error("Prottest failed to run on " + inF.string());

	ifstream inS(outF.string());
	string line;
	std::smatch m;
	std::regex e ("Best model according to.*: (.*)");
	while (getline(inS, line))
	{
		if (regex_match(line, m, e))
		{
			string bestModel= m[1];
			removeSpaces(bestModel);
			return bestModel;
		}
	}
	return "";
}


string
makeRunRaxmlCmd(const fs::path &inF, const fs::path &outF, const fs::path &wDir, const std::string &model, int nThreads, const std::string &type, const fs::path &logFile)
{
	vector<string> tokens = split(model, "+");
	std::set<string> params;
	for (auto &token : tokens)
		params.insert(token);
	auto itEnd = params.end();
	bool gFound = true;
	string raxmlFormat = type;

	// Gamma needed for raxml
	raxmlFormat += "GAMMA";
	gFound = true;

	if ((params.find("I") != itEnd) && (!gFound))
		raxmlFormat += "I";
	raxmlFormat += tokens[0];
	if (params.find("F") != itEnd)
		raxmlFormat += "F";

	string logF;
	if (logFile.empty())
		logF =  " >> /dev/null";
	else
		logF = " >> " + logFile.string() + " 2>&1";
	return "RaxML -m " + raxmlFormat + " -p 12345  -s " + inF.string() + " -w " + wDir.string() + " -n " + outF.string() + " -T " + to_string(nThreads) + " -x 20384 -# autoMRE -fa " + logF;
}


void
runRaxml(const fs::path &inF, const fs::path &outF, const fs::path &wDir, const std::string &model, int nThreads, const std::string &type, const fs::path &logFile)
{
	string command = makeRunRaxmlCmd(inF, outF, wDir, model, nThreads, type, logFile);
	if (system(command.c_str()))
		throw std::runtime_error("RaxML failed to run on " + inF.string());
}

}
