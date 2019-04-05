/*
 * domainProgs.hpp
 *
 *  Created on: October 25, 2016
 *      Author: Carsten Kemena
 *	 Copyright: 2016
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

/**
 * @file domainProgs.hpp
 * \brief This files contains functions to run different domain related programs (e.g. annotation)
 */

#ifndef DOMAIN_PROGS_HPP_
#define DOMAIN_PROGS_HPP_


#include <string>
#include <iostream>
#include <boost/filesystem.hpp>

#include "../utility/stringHelpers.hpp"
#include "../external/Input.hpp"
#include "../domain/DomainArrangementSet.hpp"

namespace fs=boost::filesystem;

namespace BioSeqDataLib
{

     /**
      * \brief Runs pfam_scan.
      * @param inputFile  The input files with the sequences to annotate.
      * @param outputFile The output file.
      * @param database   The database to use.
      * @param n_threads  The number of threads to use.
      * @param logFile    The logfile to use.
      */
     void
     runPfamScan(const fs::path &inputFile, const fs::path &outputFile, const fs::path &database, size_t nThreads=1,
         const fs::path &logFile="");

     /**
      * \brief Runs pfam_scan.
      * @param inputFile  The input files with the sequences to annotate.
      * @param outputFile The output file.
      * @param database   The database to use.
      * @param evalue     The evalue that should be used for e_dom and e_seq
      * @param n_threads  The number of threads to use.
      * @param logFile    The logfile to use.
      */
     void
     runPfamScan(const fs::path &inputFile, const fs::path &outputFile, const fs::path &database, float evalue,
         size_t nThreads=1, const fs::path &logFile="");


    /**
     * \brief Runs pfam_scan and reads results.
     * @param inputFile  The input files with the sequences to annotate.
     * @param outputFile The output file. If an empy string is provided the file will be deletet after reading.
     * @param database   The database to use.
     * @param domSet     The domain set to store the data in.
     * @param n_threads  The number of threads to use.
     * @param logFile    The logfile to use.
     */
    template<typename DomainType>
    void
    runPfamScan(const fs::path &inputFile, fs::path outputFile, const fs::path &database, DomainArrangementSet<DomainType> &domSet,
        size_t nThreads=1, const std::string &logFile="")
    {
        bool makeTmp = outputFile.empty();
        if (makeTmp)
            outputFile = boost::filesystem::unique_path();
        runPfamScan(inputFile, outputFile, database, nThreads, logFile);
        domSet.read(outputFile);
        if (makeTmp)
            fs::remove(outputFile);
    }


    /**
     * \brief Runs pfam_scan and reads results.
     * @param inputFile  The input files with the sequences to annotate.
     * @param outputFile The output file. If an empy string is provided the file will be deletet after reading.
     * @param database   The database to use.
     * @param evalue     The evalue that should be used for e_dom and e_seq
     * @param domSet     The domain set to store the data in.
     * @param n_threads  The number of threads to use.
     * @param logFile    The logfile to use.
     */
    template<typename DomainType>
    void
    runPfamScan(const fs::path &inputFile, fs::path outputFile, const fs::path &database, float evalue, DomainArrangementSet<DomainType> &domSet,
        size_t nThreads=1, const std::string &logFile="")
    {
        bool makeTmp = outputFile.empty();
        if (makeTmp)
            outputFile = boost::filesystem::unique_path();
        runPfamScan(inputFile, outputFile, database, evalue, nThreads, logFile);
        domSet.read(outputFile);
        if (makeTmp)
            fs::remove(outputFile);
    }


    struct DomainInfo
    {
    	std::string name;
    	std::string type;
    	std::string clan;

    	DomainInfo(const std::string &n, const std::string &t, const std::string &c) : name(n), type(t), clan(c)
    	{}


    };


    void
    readDomainInfo(const fs::path &pfamInfo, std::map<std::string, DomainInfo> &info);


}





#endif /* DOMAIN_PROGS_HPP_ */
