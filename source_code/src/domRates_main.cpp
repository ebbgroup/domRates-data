/*
 * DomRates is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DomRates is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DomRates.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <vector>
#include <map>
#include <numeric>

// boost header
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// BioSeqDataLib header
#include "../libs/BioSeqDataLib/src/utility/Matrix.hpp"
#include "../libs/BioSeqDataLib/src/DomainModule.hpp"
#include "../libs/BioSeqDataLib/src/phylogeny/PhylogeneticTree.hpp"
#include "../libs/BioSeqDataLib/src/phylogeny/fitch.hpp"
#include "../libs/BioSeqDataLib/src/phylogeny/dollo.hpp"

//DomRates header
#include "domRates.hpp"
#include "helperStructs.hpp"
#include "version.hpp"


namespace BSDL = BioSeqDataLib;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using std::string;
using std::vector;
using std::set;
using std::map;
using std::tuple;
using std::cout;
using std::ofstream;
using std::cerr;
using std::endl;
using std::fixed;
using std::setprecision;

int
main(int argc, char *argv[]) {

    string treeFile;
    string outgroup;
    string ending;
    string lca;
    string domrates_param_str;
    fs::path annotationDirectory;
    fs::path outFile;
    fs::path addOut;
    bool detailed;
    unsigned int nthreads;

    // argument parsing
    std::string DomRatesVersion(std::string(STR(MAJOR_VERSION)) + "." + std::string(STR(MINOR_VERSION)) + "." + std::string(STR(PATCH_VERSION)) );
    po::options_description allOpts(
            "DomRates version " + DomRatesVersion + " (C) 2016-2019 Elias Dohmen\nThis program comes with ABSOLUTELY NO WARRANTY;\n\nAllowed options are displayed below.");
    po::options_description general("General options");
    general.add_options()
            ("help,h", "Produces this help message")
            ("version,v", "Shows program version")
            ("tree,t", po::value<string>(&treeFile)->required(), "The phylogenetic tree in newick format.")
            ("annotationDirectory,a", po::value<fs::path>(&annotationDirectory)->required(),
             "A directory with all domain annotation files for species in the tree. Note that species names in tree and the regarding domain annotation file names have to be the same.")
            ("outgroup,g", po::value<string>(&outgroup)->required(),
             "The name of the outgroup as it is labeled in the tree.")
            ("ending,e", po::value<string>(&ending)->default_value(".dom"),
             "The filename extension of your domain annotation files.")
            ("out,o", po::value<fs::path>(&outFile),
             "The output file. If no output file is chosen, results will be printed to console.")
            ("statistics,s", po::value<fs::path>(&addOut),
             "File to store additional information (such as number of events per node in the tree). Additional information is just stored in file, if specified.")
            ("node,n", po::value<string>(&lca),
             "If two species names devided by ':' are provided, all arrangements involved in rearrangement events at the node representing the last common ancestor of both species will be listed in the statistics file.\n"
             "Just usable if statistics file (-s parameter) is set."
             "Example for use: '-n Drosophila_melanogaster:Caenorhabditis_elegans'")
            ("detailed,d", po::value<bool>(&detailed)->default_value(false)->zero_tokens(),
             "If this parameter is set, the output files also contain statistics about identical arrangements that have not changed. i.e. the arrangement stays conserved, and complex solutions, i.e. the rearrangement event leading to the new arrangement cannot be determined. (This can heavily increase file size.)")
            ("threads,p", po::value<unsigned int>(&nthreads)->default_value(1),
             "Number of parallel threads to use for computation.");

    allOpts.add(general);

    try {
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(allOpts).run(), vm);
        if (vm.count("help")) {
            cout << allOpts << "\n";
            return EXIT_SUCCESS;
        }

        if (vm.count("version")) {
            cout << "DomRates version " + DomRatesVersion + " (C) 2016-2019 Elias Dohmen" << "\n";
            return EXIT_SUCCESS;
        }

        if (vm.count("node") && (!vm.count("statistics"))) {
            throw po::error(
                    string("Missing option: node (-n) was used without statistics (-s) specified. Please specify -s option."));
        }

        po::notify(vm);
    }
    catch (boost::program_options::error &e) {
        cerr << "An error occurred parsing the commandline: \n";
        cerr << e.what() << "\n";
        cerr << "Please use -h/--help for more information.\n";
        return EXIT_FAILURE;
    }

    try {
        domrates_param_str = "domRates -t " + treeFile + " -a " + annotationDirectory.string() + " -g " + outgroup + " -e " + ending + " -o " + outFile.string() + " -s " + addOut.string() + " -n " + lca + " -d " + std::to_string(detailed) + " -p " + std::to_string(nthreads);
        analyseDomRates(treeFile, annotationDirectory, outgroup, ending, outFile, addOut, lca, detailed, nthreads, domrates_param_str);
    }
    catch ( const std::exception& e ) {
        cerr << "An error occured during the DomRates run: \n";
        cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

