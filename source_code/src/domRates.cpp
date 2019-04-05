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
#include <stdexcept>
#include <regex>
#include <ctime>

// boost header
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// BioSeqDataLib header
#include "../libs/BioSeqDataLib/src/DomainModule.hpp"
#include "../libs/BioSeqDataLib/src/external/Output.hpp"
#include "../libs/BioSeqDataLib/src/phylogeny/PhylogeneticTree.hpp"
#include "../libs/BioSeqDataLib/src/phylogeny/Tree.hpp"
#include "../libs/BioSeqDataLib/src/phylogeny/fitch.hpp"
#include "../libs/BioSeqDataLib/src/phylogeny/dollo.hpp"

//DomRates header
#include "domRates.hpp"
#include "helperStructs.hpp"
#include "version.hpp"

namespace BSDL = BioSeqDataLib;
namespace fs = boost::filesystem;

    using std::string;
    using std::vector;
    using std::set;
    using std::map;
    using std::cout;
    using std::ofstream;
    using std::endl;
    using std::fixed;
    using std::setprecision;

fs::path alter_filename(const fs::path  &inp, const string &appendix) {
    boost::filesystem::path nam = inp;
    boost::filesystem::path cpof = nam;

    if (!nam.extension().empty()) {
        string new_filename = nam.stem().string() + appendix + nam.extension().string();

        if (!cpof.remove_leaf().empty()) {
            nam.remove_leaf() += "/" + new_filename;
        }
        else {
            nam.remove_leaf() += new_filename;
        }
    }
    else {
        if (!cpof.remove_leaf().empty()) {
            nam.remove_leaf() += "/" + nam.stem().string() + appendix + ".txt";
        }
        else {
            nam +=  appendix + ".txt";
        }
    }

    return nam;
}

std::pair<unsigned int, vector<std::pair<string,string> > >
check_fusion(const map<vector<string>,unsigned int> & domorder, const map<unsigned int,vector<string>> & posorder, const vector<int> & pNode, const unsigned int &posi)
{
    unsigned int fus_even = 0;

    vector<string> new_rrgmnt = posorder.find(posi)->second;
    vector<std::pair<string,string>> all_fusions;

    for(unsigned int sep=1;sep<new_rrgmnt.size();++sep)
    {
        vector<string> partONE(new_rrgmnt.begin(),new_rrgmnt.begin()+sep);
        vector<string> partTWO(new_rrgmnt.begin()+sep,new_rrgmnt.end());

        bool foundONE = false;
        bool foundTWO = false;

        if (domorder.count(partONE))
        {
            unsigned int arr_pos = domorder.at(partONE);
            if (pNode[arr_pos] == 1) foundONE=true;
        }
        if (domorder.count(partTWO))
        {
            unsigned int arr_pos = domorder.at(partTWO);
            if (pNode[arr_pos] == 1) foundTWO=true;
        }

        if (foundONE and foundTWO)
        {
            ++fus_even;
            string str_part1 = boost::algorithm::join(partONE, " ");
            string str_part2 = boost::algorithm::join(partTWO, " ");
            all_fusions.emplace_back(str_part1, str_part2);
        }
    }

    return std::pair<unsigned int, vector<std::pair<string,string>>>(fus_even, all_fusions);
}

vector<unsigned int>
check_fission_termLoss_event(const map<vector<string>,unsigned int> & domorder, const map<unsigned int,vector<string>> & posorder, vector<std::pair<unsigned int, unsigned int> > & fission_pairs, vector<unsigned int> & termLoss_pairs, const vector<int> & cNode, const vector<int> & pNode, const unsigned int &posi)
{
    unsigned int fission_events = 0;
    unsigned int termLoss_events = 0;
    vector<string> act_rrngmnt = posorder.find(posi)->second;
    auto size_act_rrngmnt = act_rrngmnt.size();

    for(unsigned int parent_present_pos=0; parent_present_pos<pNode.size(); ++parent_present_pos)
    {
        if (pNode[parent_present_pos] == 1)
        {
            if(posorder.find(parent_present_pos)->second.size() > size_act_rrngmnt)
            {
                vector<string> subarrangement_front(posorder.find(parent_present_pos)->second.begin(),
                                                    posorder.find(parent_present_pos)->second.begin() + size_act_rrngmnt);
                vector<string> subarrangement_back(posorder.find(parent_present_pos)->second.end() - size_act_rrngmnt,
                                                   posorder.find(parent_present_pos)->second.end());

                if (act_rrngmnt == subarrangement_front)
                {
                    vector<string> second_part_front(posorder.find(parent_present_pos)->second.end() - (posorder.find(parent_present_pos)->second.size() - size_act_rrngmnt),
                                                     posorder.find(parent_present_pos)->second.end());

                    // does second split product exist in child node?
                    if(domorder.count(second_part_front))
                    {
                        if (cNode[domorder.at(second_part_front)] == 1)
                        {
                            ++fission_events;
                            // pair.first = second part of the fission; pair.second = parental full arrangement
                            fission_pairs.emplace_back(domorder.at(second_part_front), parent_present_pos);
                        }
                        else
                        {
                            ++termLoss_events;
                            termLoss_pairs.push_back(parent_present_pos);
                        }
                    }
                    else
                    {
                        ++termLoss_events;
                        termLoss_pairs.push_back(parent_present_pos);
                    }
                }
                else if (act_rrngmnt == subarrangement_back)
                {
                    vector<string> second_part_back(posorder.find(parent_present_pos)->second.begin(),
                                                    posorder.find(parent_present_pos)->second.end() - size_act_rrngmnt);

                    // does second split product exist in child node?
                    if(domorder.count(second_part_back))
                    {
                        if (cNode[domorder.at(second_part_back)] == 1)
                        {
                            ++fission_events;
                            // pair.first = second part of the fission; pair.second = parental full arrangement
                            fission_pairs.emplace_back(domorder.at(second_part_back), parent_present_pos);
                        }
                        else
                        {
                            ++termLoss_events;
                            termLoss_pairs.push_back(parent_present_pos);
                        }
                    }
                    else
                    {
                        ++termLoss_events;
                        termLoss_pairs.push_back(parent_present_pos);
                    }
                }
            }
        }
    }

    return {fission_events, termLoss_events};
}


std::pair<unsigned int, string>
check_termGain(const map< unsigned int,vector<string>> & posorder, const vector<int> & singleDom_parent_vec, const map<string,unsigned int> & single_domainorder, const vector<int> & parent_states, const unsigned int &posi)
{
    unsigned int term_gain_eve = 0;
    auto size_parent_arrangement = posorder.find(posi)->second.size()-1;

    for(unsigned int k = 0; k<parent_states.size(); ++k)
    {
        if(parent_states[k] == 1)
        {
            if(posorder.find(k)->second.size() == size_parent_arrangement)
            {
                vector<string> subarrangement_front(posorder.find(posi)->second.begin(),
                                                    posorder.find(posi)->second.begin() + size_parent_arrangement);
                vector<string> subarrangement_back(posorder.find(posi)->second.end() - size_parent_arrangement,
                                                   posorder.find(posi)->second.end());
                if(subarrangement_front == posorder.find(k)->second)
                {
                    if(singleDom_parent_vec[single_domainorder.find(posorder.find(posi)->second.back())->second] == -1) return std::pair<unsigned int, string>(++term_gain_eve, posorder.find(posi)->second.back());
                }
                if(subarrangement_back == posorder.find(k)->second)
                {
                    if(singleDom_parent_vec[single_domainorder.find(posorder.find(posi)->second.front())->second] == -1) return std::pair<unsigned int, string>(++term_gain_eve, posorder.find(posi)->second.front());
                }
            }
        }
    }

    return std::pair<unsigned int, string>(term_gain_eve, "");
}

unsigned int
findLCA(const string &lca, const BSDL::PhylogeneticTree<vector<int> > & nTree)
{

    string spec1, spec2;
    vector<string> specvec;
    boost::algorithm::split(specvec, lca, boost::is_any_of(":"));
    if (specvec.size() == 2) {
        spec1 = specvec[0];
        spec2 = specvec[1];
    }
    else {
        throw std::runtime_error("Error (-n option): Please check if exactly two species provided and separated by ':' (e.g.: -n Drosophila_melanogaster:Caenorhabditis_elegans)");
    }

    const BSDL::TreeNodePhylo<vector<int> >* spen1 = nullptr;
    const BSDL::TreeNodePhylo<vector<int> >* spen2 = nullptr;
    unsigned int lca_specs_f = 0;

    for(auto lcaNode=nTree.preorderBegin(); lcaNode!= nTree.preorderEnd(); ++lcaNode) {
        if (lcaNode->name == spec1 || lcaNode->name == spec2) {
            if (lca_specs_f == 0) {
                spen1 = &*lcaNode;
                ++lca_specs_f;
            }
            else if (lca_specs_f == 1) {
                spen2 = &*lcaNode;
                ++lca_specs_f;
            }
            else {
                break;
            }
        }
    }

    // check if both species have been found in tree
    if (spen1 == nullptr || spen2 == nullptr) {
        throw std::runtime_error("Error (-n option): Couldn't find both species in the provided tree.");
    }

    set<unsigned int> pot_lcas;

    do {
        if (spen1->parent() != nullptr) {
            if (pot_lcas.count(spen1->parent()->id)) {
                return spen1->parent()->id;
            }
            else {
                pot_lcas.insert(spen1->parent()->id);
                spen1 = spen1->parent();
            }
        }
        if (spen2->parent() != nullptr) {
            if (pot_lcas.count(spen2->parent()->id)) {
                return spen2->parent()->id;
            }
            else {
                pot_lcas.insert(spen2->parent()->id);
                spen2 = spen2->parent();
            }
        }
    } while (spen1->parent() != nullptr && spen2->parent() != nullptr);
    throw std::runtime_error("Error (-n option): Couln't determine last common ancestor of both species.");
}

std::pair<posOrderMaps, eventMaps>
saveDomData(BSDL::PhylogeneticTree<vector<int> > & nTree, BSDL::PhylogeneticTree<vector<int> > & singleDomTree, const fs::path &annotationDirectory, const string &outgroup, const string &ending)
{
    eventMaps treeEvents;
    posOrderMaps pomaps;
    unsigned int counter = 0;
    unsigned int counter2 = 0;

    auto bNode = singleDomTree.preorderBegin();

    for (auto aNode=nTree.preorderBegin(); aNode!=nTree.preorderEnd(); ++aNode)
    {

        treeEvents.treemap[aNode->id] = &*aNode;
        treeEvents.id_to_tree[bNode->id] = &*bNode;
        treeEvents.identities_node[aNode->id] = 0;
        treeEvents.events_per_node[aNode->id].assign(6, 0);

        if(aNode->isLeaf())
        {
            if (!pomaps.domainorder.empty()) {
                aNode->data.assign(pomaps.domainorder.size(), -1);
                bNode->data.assign(pomaps.single_domainorder.size(), -1);
            }

            // check if outgroup exists in tree and is located closest to root
            if (aNode->parent()->parent() == nullptr && aNode->name != outgroup) {
                throw std::runtime_error("Error (-g / --outgroup): Please check if an outgroup with this name exists in the tree and if its branch is closest to the root");
            }

            BSDL::DomainArrangementSet<BSDL::Domain> arrangementSet;
            fs::path nafile = aNode->name + ending;
            fs::path filepan = annotationDirectory / nafile;

            arrangementSet.read(filepan);
            arrangementSet.solveDbOverlaps({BSDL::DomainDB::pfam, BSDL::DomainDB::superfamily, BSDL::DomainDB::gene3d, BSDL::DomainDB::unknown},10,0.1);

            for (auto & arrangement : arrangementSet)
            {
                auto collapsed_arrangement = arrangement.second;
                collapsed_arrangement.collapse();
                vector<string> nDomVec;
                string sca = collapsed_arrangement.str();
                boost::split(nDomVec, sca, boost::is_any_of("- "), boost::token_compress_on);

                if (!pomaps.domainorder.count(nDomVec)) {
                    pomaps.posorder[counter] = nDomVec;
                    pomaps.domainorder[nDomVec] = counter++;
                    aNode->data.resize(pomaps.domainorder.size(),1);
                }
                else {
                    aNode->data[pomaps.domainorder[nDomVec]] = 1;
                }


                for(auto & single_dom : collapsed_arrangement)
                {
                    string ssd = single_dom.accession();
                    if (!pomaps.single_domainorder.count(ssd)) {
                        pomaps.single_posorder[counter2] = ssd;
                        pomaps.single_domainorder[ssd] = counter2++;
                        bNode->data.resize(pomaps.single_domainorder.size(),1);
                    }
                    else {
                        bNode->data[pomaps.single_domainorder[single_dom.accession()]] = 1;
                    }
                }
            }
        }
        ++bNode;
    }

    bNode = singleDomTree.preorderBegin();
    for (auto aNode=nTree.preorderBegin(); aNode!=nTree.preorderEnd(); ++aNode) {
        if(aNode-> isLeaf()) {
            aNode->data.resize(pomaps.domainorder.size(), -1);
            bNode->data.resize(pomaps.single_domainorder.size(), -1);
        }
        ++bNode;
    }

    return std::pair<posOrderMaps, eventMaps>(pomaps, treeEvents);
}


std::pair<solutionTypes, eventTypes>
eventReconstruction(const posOrderMaps &pomaps, eventMaps &emaps, const unsigned int &nthreads) {
    solutionTypes sTypes;
    eventTypes eTypes;

    float &fusion = eTypes.fusion;
    float &fission = eTypes.fission;
    float &termGain = eTypes.termGain;
    float &termLoss = eTypes.termLoss;
    float &singleDomGain = eTypes.singleDomGain;
    float &singleDomLoss = eTypes.singleDomLoss;
    unsigned int &identities_total = eTypes.identities_total;

    float &exact_solution = sTypes.exact_solution;
    float &non_ambiguous_solution = sTypes.non_ambiguous_solution;
    float &ambiguous_solution = sTypes.ambiguous_solution;
    float &complex_solution = sTypes.complex_solution;

    vector<string> &event_listing = emaps.event_listing;
    vector<string> &identities_listing = emaps.identities_listing;
    vector<string> &complex_listing = emaps.complex_listing;

    #pragma omp declare reduction (merge : std::vector<string> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

    #pragma omp parallel num_threads(nthreads)
    {
        #pragma omp for reduction(+:fusion, fission, termGain, termLoss, singleDomGain, singleDomLoss, exact_solution, non_ambiguous_solution, ambiguous_solution, complex_solution, identities_total) reduction(merge: event_listing, complex_listing, identities_listing)
        for (unsigned int ind = 0; ind < emaps.treemap.size(); ++ind) {

            BSDL::TreeNodePhylo<vector<int> >* actNode = emaps.treemap.at(ind);

            // don't count events if current node is outgroup or root or first successor of root
            // first successor of root would count all differences between outgroup and rest of the tree
            if (actNode->parent() != nullptr and actNode->parent()->parent() != nullptr) {

                vector<int> parentNode_data = actNode->parent()->data;
                vector<int> actNode_data = actNode->data;
                set<unsigned int> certain_fission;
                // helper_nonambig_fissions[pos second fission part] = pair< pos first fission part, vector<fission part positions> >
                map<unsigned int, std::pair< unsigned int, vector<std::pair<unsigned int, unsigned int> > > > helper_nonambig_fissions;

                // create set with single domains for parent and current node
                vector<int> singleDom_actNode_vec = emaps.id_to_tree.at(actNode->id)->data;
                vector<int> singleDom_parent_vec = emaps.id_to_tree.at(actNode->parent()->id)->data;

                for (unsigned int i = 0; i < pomaps.domainorder.size(); ++i) {
                    int parental_state = parentNode_data[i];
                    int child_state = actNode_data[i];

                    unsigned int fusion_event = 0;
                    unsigned int fission_event = 0;
                    unsigned int termGain_event = 0;
                    unsigned int termLoss_event = 0;
                    unsigned int singleDomGain_event = 0;
                    unsigned int singleDomLoss_event = 0;

                    vector<std::pair<string, string> > fusions_mark;
                    string termGain_dom;
                    vector<std::pair<unsigned int, unsigned int> > fission_pairs;
                    vector<unsigned int> termLoss_pairs;

                    // maintained arrangements check
                    if (parental_state == 1 && child_state == 1) {
                        identities_total++;

                        string solstr = std::to_string(actNode->id) + "\tmaintained\tmaintained\t" + boost::algorithm::join(pomaps.posorder.at(i), " ") + "\t" + boost::algorithm::join(pomaps.posorder.at(i), " ");
                        identities_listing.push_back(solstr);

                        #pragma omp atomic
                        emaps.identities_node[actNode->id]++;
                    }
                    else if (parental_state != child_state and !certain_fission.count(i)) {

                        bool lost_arr_run = false;

                        if (parental_state < child_state) {

                            if (pomaps.posorder.at(i).size() == 1 and
                                singleDom_parent_vec[pomaps.single_domainorder.at(pomaps.posorder.at(i).at(0))] ==
                                -1) {
                                ++singleDomGain_event; // single domain emergence
                            }
                            else {
                                if (pomaps.posorder.at(i).size() > 1) {
                                    std::pair<unsigned int, vector<std::pair<string, string> > > fusion_res = check_fusion(
                                            pomaps.domainorder, pomaps.posorder, parentNode_data, i); // fusion

                                    fusion_event = fusion_res.first;
                                    fusions_mark = fusion_res.second;

                                    if (fusion_event == 0) {
                                        std::pair<unsigned int, string> terGai_res = check_termGain(
                                                pomaps.posorder,
                                                singleDom_parent_vec,
                                                pomaps.single_domainorder,
                                                parentNode_data,
                                                i); // terminal emergence
                                        termGain_event = terGai_res.first;
                                        termGain_dom = terGai_res.second;
                                    }
                                }

                                vector<unsigned int> fis_answ = check_fission_termLoss_event(
                                        pomaps.domainorder, pomaps.posorder, fission_pairs, termLoss_pairs,
                                        actNode_data, parentNode_data, i); // fission and terminal loss
                                fission_event = fis_answ[0];
                                termLoss_event = fis_answ[1];
                            }
                        } else if (parental_state > child_state) {
                            lost_arr_run = true;
                            if (pomaps.posorder.at(i).size() == 1 and
                                singleDom_actNode_vec[pomaps.single_domainorder.at(pomaps.posorder.at(i).at(0))] ==
                                -1) {
                                lost_arr_run = false;
                                ++singleDomLoss_event; // single domain loss
                            }
                        }

                        bool non_ambig_fission = false;

                        unsigned int num_events = 0;
                        if (fusion_event >= 1) ++num_events;
                        if (fission_event >= 1) ++num_events;
                        if (termGain_event >= 1) ++num_events;
                        if (termLoss_event >= 1) ++num_events;
                        if (singleDomGain_event >= 1) ++num_events;
                        if (singleDomLoss_event >= 1) ++num_events;

                        if (num_events == 1) {
                        // ^ just exact and non-ambiguous solutions v
                            if (fusion_event >= 1) {
                                ++fusion;
                                emaps.events_per_node[actNode->id][0]++;
                                string solution_str;
                                if (fusion_event == 1) {
                                    solution_str = "\texact solution\tfusion\t";
                                } else {
                                    solution_str = "\tnon-ambiguous solution\tfusion\t";
                                }
                                for (auto & fus : fusions_mark) {
                                    string solstr = std::to_string(actNode->id) + solution_str +
                                                    boost::algorithm::join(pomaps.posorder.at(i), " ") +
                                                    "\t" + fus.first + " + " + fus.second;
                                    event_listing.push_back(solstr);
                                }
                            } else if (fission_event >= 1) {
                                // ^ fission events have a second subarrangement that can be checked later and lead to a different solution v
                                string solution_str;
                                bool later_occurence = false;

                                if (fission_event == 1) {
                                    solution_str = "\texact solution\tfission\t";
                                    certain_fission.insert(fission_pairs.at(0).first);
                                } else {
                                    solution_str = "\tnon-ambiguous solution\tfission\t";

                                    for (auto & fis : fission_pairs) {
                                        if ((parentNode_data.at(fis.first) == -1) and (actNode_data.at(fis.first) == 1)) {
                                            later_occurence = true;
                                            non_ambig_fission = true;
                                            // take into account the case if the fission pair is checked later, but can be ambiguous, instead of like in this run a non-ambiguous solution
                                            helper_nonambig_fissions[fis.first] = std::pair<unsigned int, vector<std::pair<unsigned int, unsigned int> > >(i, fission_pairs);
                                        }
                                    }
                                }
                                if (!later_occurence) {
                                    // if the second subarrangement is not checked later in the algorithm this solution is added to output
                                    ++fission;
                                    emaps.events_per_node[actNode->id][1]++;
                                    for (auto & fis : fission_pairs) {
                                        string arrp1 = boost::algorithm::join(pomaps.posorder.at(i), " ");
                                        string arrp2 = boost::algorithm::join(pomaps.posorder.at(fis.first),
                                                                              " ");
                                        string parr = boost::algorithm::join(pomaps.posorder.at(fis.second),
                                                                             " ");

                                        if ((arrp1 + " " + arrp2) == parr) {
                                            string solstr = std::to_string(actNode->id).append(
                                                    solution_str).append(arrp1).append(" | ").append(
                                                    arrp2).append("\t").append(parr);
                                            event_listing.push_back(solstr);
                                        } else {
                                            string solstr = std::to_string(actNode->id).append(
                                                    solution_str).append(arrp2).append(" | ").append(
                                                    arrp1).append("\t").append(parr);
                                            event_listing.push_back(solstr);
                                        }
                                    }
                                }
                            } else if (termLoss_event >= 1) {
                                ++termLoss;
                                emaps.events_per_node[actNode->id][2]++;
                                string solution_str;
                                if (termLoss_event == 1) {
                                    solution_str = "\texact solution\tterminal loss\t";
                                } else {
                                    solution_str = "\tnon-ambiguous solution\tterminal loss\t";
                                }
                                for (auto & tl : termLoss_pairs) {
                                    string solstr = std::to_string(actNode->id) + solution_str +
                                                    boost::algorithm::join(pomaps.posorder.at(i), " ") +
                                                    "\t" +
                                                    boost::algorithm::join(pomaps.posorder.at(tl), " ");
                                    event_listing.push_back(solstr);
                                }
                            } else if (termGain_event >= 1) {
                                ++termGain;
                                emaps.events_per_node[actNode->id][3]++;
                                string solution_str;
                                if (termGain_event == 1) {
                                    solution_str = "\texact solution\tterminal emergence\t";
                                } else {
                                    solution_str = "\tnon-ambiguous solution\tterminal emergence\t";
                                }

                                string domarrstr = boost::algorithm::join(pomaps.posorder.at(i),
                                                                          " ");
                                string prevarr = domarrstr.erase(domarrstr.find(termGain_dom),
                                                                 termGain_dom.size());
                                string solstr = std::to_string(actNode->id).append(solution_str).append(
                                        boost::algorithm::join(pomaps.posorder.at(i), " ")).append("\t").append(prevarr);
                                event_listing.push_back(solstr);
                            } else if (singleDomLoss_event >= 1) {
                                ++singleDomLoss;
                                emaps.events_per_node[actNode->id][4]++;
                                string solution_str;
                                if (singleDomLoss_event == 1) {
                                    solution_str = "\texact solution\tsingle domain loss\t\t";
                                } else {
                                    solution_str = "\tnon-ambiguous solution\tsingle domain loss\t\t";
                                }
                                string solstr = std::to_string(actNode->id) + solution_str +
                                                pomaps.posorder.at(i).at(0);
                                event_listing.push_back(solstr);
                            } else if (singleDomGain_event == 1) {
                                ++singleDomGain;
                                emaps.events_per_node[actNode->id][5]++;

                                string solution_str = "\texact solution\tsingle domain emergence\t";
                                string solstr = std::to_string(actNode->id) + solution_str + pomaps.posorder.at(i).at(0) + "\t";
                                event_listing.push_back(solstr);
                            }

                            if (fusion_event + fission_event + termGain_event + termLoss_event +
                                singleDomGain_event + singleDomLoss_event == 1) {
                                ++exact_solution;
                            } else {
                                if (!non_ambig_fission) {
                                    // non-ambiguous fission part that would be checked later could also become an exact solution
                                    ++non_ambiguous_solution;
                                }
                            }
                        } else if (num_events > 1) {
                        // ^ ambiguous solutions v
                            if (helper_nonambig_fissions.count(i)) {
                                // ^ if in this run the second fission part led to an ambiguous solution take the more precise non-ambiguous solution results from the first fission part
                                ++fission;
                                emaps.events_per_node[actNode->id][1]++;
                                ++non_ambiguous_solution;
                                for (auto fisp = helper_nonambig_fissions.at(i).second.begin(); fisp != helper_nonambig_fissions.at(i).second.end(); ++fisp) {
                                    string arrp1 = boost::algorithm::join(pomaps.posorder.at(helper_nonambig_fissions.at(i).first), " ");
                                    string arrp2 = boost::algorithm::join(pomaps.posorder.at(fisp->first), " ");
                                    string parr = boost::algorithm::join(pomaps.posorder.at(fisp->second), " ");

                                    if ((arrp1 + " " + arrp2) == parr) {
                                        string solstr = std::to_string(actNode->id).append(
                                                "\tnon-ambiguous solution\tfission\t").append(arrp1).append(" | ").append(
                                                arrp2).append("\t").append(parr);
                                        event_listing.push_back(solstr);
                                    } else {
                                        string solstr = std::to_string(actNode->id).append(
                                                "\tnon-ambiguous solution\tfission\t").append(arrp2).append(" | ").append(
                                                arrp1).append("\t").append(parr);
                                        event_listing.push_back(solstr);
                                    }
                                }
                            }
                            else {
                                if (fission_event == 1) {
                                    if (!(parentNode_data.at(fission_pairs.at(0).first) == -1 and
                                          actNode_data.at(fission_pairs.at(0).first) == 1)) {
                                        // skipped if the second part of the fission is still checked as new arrangement and solution and event determination will be done later
                                        ++ambiguous_solution;
                                        if (fusion_event >= 1) {
                                            for (auto & fus : fusions_mark) {
                                                string solstr = std::to_string(actNode->id).append(
                                                        "\tambiguous solution\tfusion\t").append(
                                                        boost::algorithm::join(pomaps.posorder.at(i), " ")).append(
                                                        "\t").append(fus.first).append(" + ").append(fus.second);
                                                complex_listing.push_back(solstr);
                                            }
                                        }
                                        if (fission_event >= 1) {
                                            for (auto & fis : fission_pairs) {
                                                string arrp1 = boost::algorithm::join(pomaps.posorder.at(i), " ");
                                                string arrp2 = boost::algorithm::join(pomaps.posorder.at(fis.first),
                                                                                      " ");
                                                string parr = boost::algorithm::join(pomaps.posorder.at(fis.second),
                                                                                     " ");

                                                if ((arrp1 + " " + arrp2) == parr) {
                                                    string solstr = std::to_string(actNode->id).append(
                                                            "\tambiguous solution\tfission\t").append(arrp1).append(
                                                            " | ").append(
                                                            arrp2).append("\t").append(parr);
                                                    complex_listing.push_back(solstr);
                                                } else {
                                                    string solstr = std::to_string(actNode->id).append(
                                                            "\tambiguous solution\tfission\t").append(arrp2).append(
                                                            " | ").append(
                                                            arrp1).append("\t").append(parr);
                                                    complex_listing.push_back(solstr);
                                                }
                                            }
                                        }
                                        if (termLoss_event >= 1) {
                                            for (auto &tl : termLoss_pairs) {
                                                string solstr = std::to_string(actNode->id).append(
                                                        "\tambiguous solution\tterminal loss\t").append(
                                                        boost::algorithm::join(pomaps.posorder.at(i), " ")).append(
                                                        "\t").append(
                                                        boost::algorithm::join(pomaps.posorder.at(tl), " "));
                                                complex_listing.push_back(solstr);
                                            }
                                        }
                                    }
                                } else if (fission_event >= 1) {
                                    bool later_occurence = false;
                                    for (auto & fis : fission_pairs) {
                                        if (parentNode_data.at(fis.first) == -1 and actNode_data.at(fis.first) == 1) {
                                            later_occurence = true;
                                        }
                                    }
                                    if (!later_occurence) {
                                        // skipped if at least one second part of the fissions is still checked as new arrangement -> in that case solution and event determination will be done in later run
                                        ++ambiguous_solution;
                                        if (fusion_event >= 1) {
                                            for (auto & fus : fusions_mark) {
                                                string solstr = std::to_string(actNode->id).append(
                                                        "\tambiguous solution\tfusion\t").append(
                                                        boost::algorithm::join(pomaps.posorder.at(i), " ")).append(
                                                        "\t").append(fus.first).append(" + ").append(fus.second);
                                                complex_listing.push_back(solstr);
                                            }
                                        }
                                        if (fission_event >= 1) {
                                            for (auto & fis : fission_pairs) {
                                                string arrp1 = boost::algorithm::join(pomaps.posorder.at(i), " ");
                                                string arrp2 = boost::algorithm::join(pomaps.posorder.at(fis.first),
                                                                                      " ");
                                                string parr = boost::algorithm::join(pomaps.posorder.at(fis.second),
                                                                                     " ");

                                                if ((arrp1 + " " + arrp2) == parr) {
                                                    string solstr = std::to_string(actNode->id).append(
                                                            "\tambiguous solution\tfission\t").append(arrp1).append(
                                                            " | ").append(
                                                            arrp2).append("\t").append(parr);
                                                    complex_listing.push_back(solstr);
                                                } else {
                                                    string solstr = std::to_string(actNode->id).append(
                                                            "\tambiguous solution\tfission\t").append(arrp2).append(
                                                            " | ").append(
                                                            arrp1).append("\t").append(parr);
                                                    complex_listing.push_back(solstr);
                                                }
                                            }
                                        }
                                        if (termLoss_event >= 1) {
                                            for (auto & tl : termLoss_pairs) {
                                                string solstr = std::to_string(actNode->id).append(
                                                        "\tambiguous solution\tterminal loss\t").append(
                                                        boost::algorithm::join(pomaps.posorder.at(i), " ")).append(
                                                        "\t").append(
                                                        boost::algorithm::join(pomaps.posorder.at(tl), " "));
                                                complex_listing.push_back(solstr);
                                            }
                                        }
                                    }
                                } else if (fission_event == 0) {
                                    ++ambiguous_solution;
                                    if (fusion_event >= 1) {
                                        for (auto &fus : fusions_mark) {
                                            string solstr = std::to_string(actNode->id).append(
                                                    "\tambiguous solution\tfusion\t").append(
                                                    boost::algorithm::join(pomaps.posorder.at(i), " ")).append(
                                                    "\t").append(fus.first).append(" + ").append(fus.second);
                                            complex_listing.push_back(solstr);
                                        }
                                    }
                                    if (termLoss_event >= 1) {
                                        for (auto &tl : termLoss_pairs) {
                                            string solstr = std::to_string(actNode->id).append(
                                                    "\tambiguous solution\tterminal loss\t").append(
                                                    boost::algorithm::join(pomaps.posorder.at(i), " ")).append(
                                                    "\t").append(
                                                    boost::algorithm::join(pomaps.posorder.at(tl), " "));
                                            complex_listing.push_back(solstr);
                                        }
                                    }
                                }
                            }
                        }
                        else if (num_events == 0 and !lost_arr_run) {
                        // ^ complex solution if no defined event type could explain the new arrangement v
                            ++complex_solution;
                            string solstr = std::to_string(actNode->id).append("\tcomplex solution\t?\t").append(boost::algorithm::join(pomaps.posorder.at(i), " ").append("\t"));
                            complex_listing.push_back(solstr);
                        }
                    }
                }
            }
        }
    }
    return {sTypes, eTypes};
}

void
summary(BSDL::PhylogeneticTree<vector<int> > &nTree, const eventMaps &emaps, const solutionTypes &sTypes, const eventTypes &eTypes, const fs::path &outFile, const fs::path &addOut, const string &lca, const unsigned int &lca_id, const bool &detailed, const string &domrates_param_str)
{

    AlgorithmPack::Output out(outFile);
    time_t now = time(0);
    char* dt = std::ctime(&now);

    // main output - overview for whole tree
    out << "# DomRates version " + string(STR(MAJOR_VERSION)) + "." + string(STR(MINOR_VERSION)) + "." + string(STR(PATCH_VERSION)) + " at " << dt << "# ";
    out << domrates_param_str << '\n';
    out << "# Solution types" << '\n';
    out << "Exact solutions: " << sTypes.exact_solution << '\n';
    out << "Non-ambiguous solutions: " << sTypes.non_ambiguous_solution << '\n';
    out << "Ambiguous solutions: " << sTypes.ambiguous_solution << '\n';
    out << "Complex solutions: " << sTypes.complex_solution << '\n';
    if (detailed) {
        out << "Maintained arrangements total: " << eTypes.identities_total << '\n';
    }
    out << '\n';

    out << "Exact and non-ambiguous solutions: " << sTypes.exact_solution + sTypes.non_ambiguous_solution << '\n';
    out << '\n';

    out << "# Event types" << '\n';
    out << "Fusions: " << eTypes.fusion << '\n';
    out << "Fissions: " << eTypes.fission << '\n';
    out << "Terminal Loss: " << eTypes.termLoss << '\n';
    out << "Terminal Emergences: " << eTypes.termGain << '\n';
    out << "Single Domain Losses: " << eTypes.singleDomLoss << '\n';
    out << "Single Domain Emergences: " << eTypes.singleDomGain << '\n';
    out << '\n';

    out << "# Event rates" << '\n';
    out << "Fusion rate: " << fixed << setprecision(2)
         << eTypes.fusion / (sTypes.exact_solution + sTypes.non_ambiguous_solution) * 100 << "%" << '\n';
    out << "Fission rate: " << fixed << setprecision(2)
         << eTypes.fission / (sTypes.exact_solution + sTypes.non_ambiguous_solution) * 100 << "%" << '\n';
    out << "Terminal Loss rate: " << fixed << setprecision(2)
         << eTypes.termLoss / (sTypes.exact_solution + sTypes.non_ambiguous_solution) * 100 << "%" << '\n';
    out << "Terminal Emergence rate: " << fixed << setprecision(2)
         << eTypes.termGain / (sTypes.exact_solution + sTypes.non_ambiguous_solution) * 100 << "%" << '\n';
    out << "Single Domain Loss rate: " << fixed << setprecision(2)
         << eTypes.singleDomLoss / (sTypes.exact_solution + sTypes.non_ambiguous_solution) * 100 << "%" << '\n';
    out << "Single Domain Emergence rate: " << fixed << setprecision(2)
         << eTypes.singleDomGain / (sTypes.exact_solution + sTypes.non_ambiguous_solution) * 100 << "%" << '\n';
    out << '\n';

    // additional output file with events per node
    if (!addOut.empty()) {
        AlgorithmPack::Output aout(addOut);

        // create ordered set with all events from vectors (vectors can be parallelised during the event calculation, but contain then nodes and events in random order)
        set<string> eventset;
        for (auto &eve : emaps.event_listing) {
            eventset.insert(eve);
        }
        if (detailed) {
            for (auto &ident : emaps.identities_listing) {
                eventset.insert(ident);
            }
            for (auto &comp : emaps.complex_listing) {
                eventset.insert(comp);
            }
        }

        aout << "# Number of events per node." << "\n";
        string headerstr_addOut = "# Node ID\t#Fusions\t#Fissions\t#TerminalLosses\t#TerminalEmergences\t#SingleDomainLosses\t#SingleDomainEmergences";
        if (detailed) {
            headerstr_addOut.append("\t#Maintained");
        }
        aout << headerstr_addOut << "\n";

        for (auto actNode = nTree.preorderBegin(); actNode != nTree.preorderEnd(); ++actNode) {
            string epn_line = std::to_string(actNode->id).append("\t").append(std::to_string(emaps.events_per_node.at(actNode->id).at(0))).append(
                    "\t").append(std::to_string(emaps.events_per_node.at(actNode->id).at(1))).append("\t").append(std::to_string(emaps.events_per_node.at(actNode->id).at(2))).append(
                            "\t").append(std::to_string(emaps.events_per_node.at(actNode->id).at(3))).append("\t").append(std::to_string(emaps.events_per_node.at(actNode->id).at(4))).append(
                                    "\t").append(std::to_string(emaps.events_per_node.at(actNode->id).at(5)));
            if (detailed) {
                epn_line.append("\t").append(std::to_string(emaps.identities_node.at(actNode->id)));
            }
            aout << epn_line << "\n";
        }

        if (!lca.empty()) {
            aout << "# Events per domain arrangement for last common ancestor of " << lca << "." << "\n";

            string headerstr_lca = "# Node-ID\tsolution type\tevent type\tnew arrangement at current node\tarrangement at parental node";
            aout << headerstr_lca << "\n";

            for (auto &eve : eventset) {
                vector<string> tmpvec;
                boost::algorithm::split(tmpvec, eve, boost::is_any_of("\t"));
                if (tmpvec[0] == std::to_string(lca_id)) {
                    aout << eve << "\n";
                }
            }
        }

        fs::path ed_addOut = alter_filename(addOut, "_epd");
        AlgorithmPack::Output asout(ed_addOut);

        string headerstr_edOut = "# Node-ID\tsolution type\tevent type\tnew arrangement at node\tarrangement at parental node";
        asout << headerstr_edOut << "\n";

        for(auto &el : eventset) {
            asout << el << "\n";
        }
    }
}

void
analyseDomRates(const string &treeFile, const fs::path &annotationDirectory, const string &outgroup, const string &ending, const fs::path &outFile, const fs::path &addOut, const string &lca, const bool &detailed, const unsigned int &nthreads, const string &domrates_param_str)
{

    // start initialisation
    BSDL::PhylogeneticTree<vector<int> > nTree;
    BSDL::PhylogeneticTree<vector<int> > singleDomTree;

    // load tree two times (for single domains and for domain arrangements)
    cout << "load tree..." << endl;
    nTree.read(treeFile);
    singleDomTree.read(treeFile);

    // check if tree is strictly bifurcating
    try {
        BSDL::isBifurcatingTree(nTree, true);
    }
    catch (...) {
        throw;
    }

    // find last common ancestor if -n option was specified
    unsigned int lca_id = 0;
    if (!lca.empty()) {
        try {
            lca_id = findLCA(lca, nTree);
        }
        catch (...){
            throw;
        }
    }

    // save all domain arrangements from annotation files in node->data of species tree (&nTree)
    cout << "read all arrangements..." << endl;
    std::pair<posOrderMaps, eventMaps> emaps;
    try {
        emaps = saveDomData(nTree, singleDomTree, annotationDirectory, outgroup, ending);
    }
    catch ( ... ) {
        throw;
    }
    posOrderMaps &pomapping = emaps.first;
    eventMaps &emapping = emaps.second;

    // reconstruction of ancestral domain states
    cout << "reconstructing ancestral states..." << endl;
    BSDL::fitch(nTree);
    BSDL::dollo(singleDomTree);

    cout << "event reconstruction..." << endl;
    std::pair<solutionTypes, eventTypes> setypes = eventReconstruction(pomapping, emapping, nthreads);
    solutionTypes &nSol = setypes.first;
    eventTypes &nEve = setypes.second;

    cout << "write summary..." << endl;
    summary(nTree, emapping, nSol, nEve, outFile, addOut, lca, lca_id, detailed, domrates_param_str);
}