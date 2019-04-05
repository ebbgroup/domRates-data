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

#ifndef SRC_DOMRATES_HPP
#define SRC_DOMRATES_HPP

#include <string>
#include <vector>
#include <map>
#include <utility>

// boost header
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// BioSeqDataLib header
#include "../libs/BioSeqDataLib/src/utility/Matrix.hpp"
#include "../libs/BioSeqDataLib/src/DomainModule.hpp"
#include "../libs/BioSeqDataLib/src/external/Output.hpp"
#include "../libs/BioSeqDataLib/src/phylogeny/PhylogeneticTree.hpp"
#include "../libs/BioSeqDataLib/src/phylogeny/fitch.hpp"
#include "../libs/BioSeqDataLib/src/phylogeny/dollo.hpp"

//DomRates header
#include "helperStructs.hpp"

namespace fs = boost::filesystem;
namespace BSDL = BioSeqDataLib;

/**
 * @brief alters a given filename and adds the given appendix to it
 * @details to a given filename the appendix is added before the file extension and if no extension is given
 * ".txt" is additionally added (e.g. outfile.txt --> outfile_appdx.txt; outfile --> outfile_appdx.txt;
 * outfile.out --> outfile_appdx.out)
 *
 * @param inp input filename which should be altered
 */
fs::path alter_filename(const fs::path  &inp, const std::string &appendix);

/**
 * @brief checks if a domain arrangement can be explained by a fusion event
 * @details a fusion event is defined as a bigger arrangement formed by exactly two subarrangements (A + B-C -> A-B-C)
 *
 * @param domorder existing domain arrangements and their index position
 * @param posorder index positions and their related domain arrangement
 * @param pNode presence/absence (1/-1) states for all domain arrangements at the parental node (order matching with domorder/posorder)
 * @param posi index position of the arrangement to be checked
 */
std::pair<unsigned int, std::vector<std::pair<std::string,std::string> > > check_fusion(const std::map<std::vector<std::string>,unsigned int> & domorder, const std::map<unsigned int, std::vector<std::string>> & posorder, const std::vector<int> & pNode, const unsigned int &posi);

/**
 * @brief checks if a domain arrangement can be explained by a fission or a terminal loss event
 * @details a fission is defined as a split of a bigger arrangement into two subarrangements
 * with both split products still present (A-B-C-D -> A-B | C-D),
 * if one of the split products is not present it is considered a terminal loss (A-B-C-D -> A-B)
 *
 * @param domorder existing domain arrangements and their index position
 * @param posorder index positions and their related domain arrangement
 * @param[out] fission_pairs all fission subarrangements; pair.first: index second subarrangement pair.second: index parent arrangement
 * @param[out] termLoss_pairs all index positions of possible terminal loss parent arrangements
 * @param cNode presence/absence (1/-1) states for all domain arrangements at the current node (order matching with domorder/posorder)
 * @param pNode presence/absence (1/-1) states for all domain arrangements at the parental node (order matching with domorder/posorder)
 * @param posi index position of the arrangement to be checked
 */
std::vector<unsigned int> check_fission_termLoss_event(const std::map<std::vector<std::string>,unsigned int> & domorder, const std::map<unsigned int, std::vector<std::string>> & posorder, std::vector<std::pair<unsigned int, unsigned int> > & fission_pairs, std::vector<unsigned int> & termLoss_pairs, const std::vector<int> & cNode, const std::vector<int> & pNode, const unsigned int &posi);

/**
 * @brief checks if a domain arrangement can be explained by a terminal emergence event
 * @details a terminal emergence is defined by addition of exactly one domain (E),
 * which did not exist in an ancestor, at either the C- or N-terminus of a domain arrangement (A-B -> A-B-E)
 *
 * @param posorder index positions and their related domain arrangement
 * @param singleDom_parent_vec presence/absence (1/-1) states for all single domains at the parent node (order matching with single_domainorder)
 * @param single_domainorder existing domain arrangements and their index position
 * @param parent_states presence/absence (1/-1) states for all domain arrangements at the parent node
 * @param posi index position of the arrangement to be checked
 */
std::pair<unsigned int, std::string> check_termGain(const std::map< unsigned int, std::vector<std::string>> & posorder, const std::vector<int> & singleDom_parent_vec, const std::map<std::string,unsigned int> & single_domainorder, const std::vector<int> & parent_states, const unsigned int &posi);

/**
 * @brief finds the last common ancestor of two given species in a tree
 *
 * @param lca two species names separated by a ":" of which the last common ancestor should be found
 * @param nTree phylogenetic tree in which the species can be found
 */
unsigned int findLCA(const std::string &lca, const BSDL::PhylogeneticTree<std::vector<int> > &nTree);

/**
 * @brief saves presence/absence (1/-1) states for domain arrangements and single domains for every leaf in a given phylogenetic tree
 *
 * @param[in|out] nTree tree to store all domain arrangements
 * @param[in|out] singleDomTree tree to store all single domains
 * @param annotationDirectory directory containing domain annotation data (e.g. PfamScan output files) for all species in the tree
 * @param outgroup the species/group to be used as outgroup (it should be located closest to root (regarding the hirarchy levels in the tree, not branch length))
 * @param ending file extension that has to be added to species names in the tree to read the related annotation file
 */
std::pair<posOrderMaps, eventMaps> saveDomData(BSDL::PhylogeneticTree<std::vector<int> > &nTree, BSDL::PhylogeneticTree<std::vector<int> > &singleDomTree, const fs::path &annotationDirectory, const std::string &outgroup, const std::string &ending);

/**
 * @brief infers six domain rearrangement event types and their frequency per node in a given phylogentic tree
 * @details the six infered event types are 1) fusion 2) fission 3) terminal loss 4) terminal emergence 5) single domain loss 6) single domain emergence
 * additionally, four solution types are distinguished 1) exact solution 2) non-ambiguous solution 3) ambiguous solution 4) complex solution
 * the different solution types are defined by the amount and types of event types that can explain a new domain arrangement
 *
 * @param pomaps data structure storing maps matching domain arrangements/single domains to their related index positions in the data set
 * @param emaps data structure storing information of matching nodeIDs to nodes and reconstructed events per node
 * @param nthreads number of threads that run event reconstruction in parallel
 */
std::pair<solutionTypes, eventTypes> eventReconstruction(const posOrderMaps &pomaps, eventMaps &emaps, const unsigned int &nthreads);

/**
 * @brief creates human readable output of different statistics and writes them to specified output files
 * @details dependent on input parameters different output is created containing differing levels of details and
 * either stored in output files or printed to the console
 *
 * @param nTree phylogenetic tree used for reconstruction of events
 * @param emaps data structure storing information of reconstructed events per node
 * @param sTypes data structure storing number of reconstructed solution types
 * @param eTypes data structure storing number of reconstructed event types
 * @param outFile name for an output file the frequency of different event and solution types is written to
 * @param addOut name for an output file containing additional statistics of the event reconstruction
 * @param lca two species names separated by ":" for whose last common ancestor reconstruction details are written to output
 * @param lca_id node-ID of last common ancestor of the species defined in the lca parameter
 * @param detailed if set, output contains information about maintained arrangements, which haven't been rearranged
 */
void summary(BSDL::PhylogeneticTree<std::vector<int> > &nTree, const eventMaps &emaps, const solutionTypes &sTypes, const eventTypes &eTypes, const fs::path &outFile, const fs::path &addOut, const std::string &lca, const unsigned int &lca_id, const bool &detailed, const std::string &domrates_param_str);

/**
 * @brief wrapper function called by the main script coordinating all necessary steps for full DomRates analysis
 * @details this function provides error handling, reading in trees and calling all necessary functions with according parameters
 *
 * @param annotationDirectory directory containing domain annotation data (e.g. PfamScan output files) for all species in the tree
 * @param outgroup the species/group to be used as outgroup (it should be located closest to root (regarding the hirarchy levels in the tree, not branch length))
 * @param ending file extension that has to be added to species names in the tree to read the related annotation file
 * @param outFile name for an output file the frequency of different event and solution types is written to
 * @param addOut name for an output file containing additional statistics of the event reconstruction
 * @param lca two species names separated by ":" for whose last common ancestor reconstruction details are written to output
 * @param detailed if set, output contains information about maintained arrangements, which haven't been rearranged
 * @param nthreads number of threads that run event reconstruction in parallel
 */
void analyseDomRates(const std::string &treeFile, const fs::path &annotationDirectory, const std::string &outgroup, const std::string &ending, const fs::path &outFile, const fs::path &addOut, const std::string &lca, const bool &detailed, const unsigned int &nthreads, const std::string &domrates_param_str);


#endif //SRC_DOMRATES_HPP
