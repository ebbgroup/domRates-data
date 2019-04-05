/*
 * runPhylo.hpp
 *
 *  Created on: Sep 15, 2015
 *      Author: ckeme_01
 */

/**
 * \file runPhylo.hpp
 * \brief Contains programs to calculate phylogenetic trees
 */

#ifndef SRC_EXTERNAL_RUNPHYLO_HPP_
#define SRC_EXTERNAL_RUNPHYLO_HPP_

// C++ header
#include <string>
#include <stdexcept>
#include <fstream>
#include <regex>
#include <set>
#include <vector>
#include <iostream>

#include <boost/filesystem.hpp>

#include "../utility/stringHelpers.hpp"

namespace BioSeqDataLib
{

/**
 * \brief List of alignment programs.
 * Programs are not included in BSDL. They need to be installed separately and need to be in the PATH variable.
 */
enum class AlnProgram
{
	/*! T-Coffee program */
	tcoffee,
	/*! Mafft program */
	mafft,
	/*! ClustalO program */
	clustalo,
	/*! MSAProbs program */
	msaprobs
};

/**
 * \brief The output formats for alignment construction.
 * Be aware that not all programs support all of the formats listed here.
 */
enum class AlnFormat
{
	/*! The fasta format */
	fasta,
	/*! The pylip format */
	phylip,
	/*! The clustal format */
	clustal
};

namespace fs=boost::filesystem;

/**
 * \brief Run trimAl, a alignment trimming program.
 * \param inF The input file.
 * \param outF The output file.
 * \param parameter The parameters to used. Default is optimization for Maximum Likelihood phylogentic reconstruction.
 * \param log The logfile
 */
void
runTrimAl(const fs::path &inF, const fs::path &outF, const std::string &parameter=" -automated1", const fs::path &log="/dev/null");

/**
 * \brief Runs prottest to determine model for phylogenetic reconstruction.
 * \param inF The input file
 * \param outF The output file
 * \param nThreads The number of threads to use
 * \param The parameter to use. By default it uses some predefined settings
 * \param log The logfile \param log The logfile
 */
std::string
runProttest(const fs::path &inF, const fs::path &outF, int nThreads, const std::string &parameter=" -AICC  -I  -G -F -IG  -JTT  -LG  -DCMut -Dayhoff -WAG -Blosum62 -VT ", const fs::path &log="/dev/null");

/**
 * \brief Run runRaxml
 * \param inF The input file.
 * \param outF The output file.
 * \param wDir The working directory.
 * \param model The model to use.
 * \param nThreads The number of threads to use.
 * \param type The type to use
 * \param log The logfile to use
 * Always used protgamma and it ignores -I as suggested by the author of RAXML.
 */
void
runRaxml(const fs::path &inF, const fs::path &outF, const fs::path &wDir, const std::string &model, int nThreads, const std::string &type="PROT", const fs::path &log="/dev/null");

/**
 *
 */
//std::string
//makeRunRaxmlCmd(const fs::path &inF, const fs::path &outF, const fs::path &wDir, const std::string &model, int nThreads, const std::string &type, const fs::path &log);

/**
 * \brief Calls an alignment program.
 * \param inF The input file.
 * \param program The alignment program to use.
 * \param outF The output File of the program.
 * \param format The output format. Be aware that not all programs support all formats.
 * \param log The log file to use.
 */
void
runAlignment(const fs::path &inF, AlnProgram program, const fs::path &outF, int nThreads=1, AlnFormat format=AlnFormat::fasta, const fs::path &log="/dev/null");

}

#endif /* SRC_EXTERNAL_RUNPHYLO_HPP_ */
