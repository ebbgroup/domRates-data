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


#ifndef DOMRATES_HELPERSTRUCTS_H
#define DOMRATES_HELPERSTRUCTS_H

#include <string>
#include <vector>
#include <map>

// BioSeqDataLib header
#include "../libs/BioSeqDataLib/src/phylogeny/PhylogeneticTree.hpp"

namespace BSDL = BioSeqDataLib;

/**
 * structure to store for a given domain arrangement data set all arrangements and single domains
 * in a given order with index positions for fast access and mapping later from simple vectors of integers representing
 * presence/absence (1/-1) states of all domains/arrangements
 */
struct
posOrderMaps {
    std::map<std::vector<std::string>,unsigned int> domainorder; // domain arrangement, index position of this arrangement
    std::map<unsigned int, std::vector<std::string> > posorder; // index position of arrangement, domain arrangement
    std::map<std::string,unsigned int> single_domainorder; // single domain, index position
    std::map<unsigned int, std::string> single_posorder; // index position, single domain
};

/**
 * structure to store mapping information of nodes and their corresponding IDs
 * as well as reconstructed events for later output
 */
struct
eventMaps {
    // domain arrangement annotation for arrangement Tree
    std::map<unsigned int, BSDL::TreeNodePhylo<std::vector<int> >* > treemap;
    // single domain annotation for singleDomTree
    std::map<unsigned int, BSDL::TreeNodePhylo<std::vector<int> >* > id_to_tree;
    // counts identical domain arrangements between node and its parent node
    std::map<unsigned int,unsigned int> identities_node;
    // counts events per node. id -> (#fusion, #fission, #terminal loss, #terminal gain, #single loss, #single gain)
    std::map<unsigned int, std::vector<unsigned int> > events_per_node;
    // summary of all solved events
    std::vector<std::string> event_listing;
    // summary of all identities (arrangements that do not change from one node to the other)
    std::vector<std::string> identities_listing;
    // summary of all complex solutions (new arrangements without explanation by any event type)
    std::vector<std::string> complex_listing;
};

/**
 * structure to count frequency of reconstructed solution types
 */
struct
solutionTypes {
    float exact_solution = 0.0; // one solution, one event type
    float non_ambiguous_solution = 0.0; // multiple solutions of same event type
    float ambiguous_solution = 0.0; // multiple solutions of different event types
    float complex_solution = 0.0; // no solution in a single step with considered event types
};

/**
 * structure to count frequency of reconstructed event types
 */
struct
eventTypes {
    float fusion  = 0.0;
    float fission = 0.0;
    float termGain= 0.0;
    float termLoss= 0.0;
    float singleDomGain = 0.0;
    float singleDomLoss = 0.0;
    unsigned int identities_total = 0; // arrangement was found in the ancestor and didn't change
};

#endif //DOMRATES_HELPERSTRUCTS_H
