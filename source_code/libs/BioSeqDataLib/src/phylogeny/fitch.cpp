//
// Created by e_dohm01 on 03.08.16.
//


#include "fitch.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

namespace BioSeqDataLib {


    /**
    * \brief Infers the state of all instances (e.g. domainArrangements) by majority rule.
    * @param child_l vector<int> with states of one of the child nodes
    * @param child_r vector<int> with states of the other child nodes
    * @return Parental state based on both child states.
    */
    vector<int>
    parentalState(vector<int> &child_l, vector<int> &child_r)
    {
        /*
         * The states for the instances (e.g. domainArrangements) are given as integers
         * taking only the values +1 (existent), 0 (status unclear) and -1 (not existent)
         * Arithmetics:
         *  1 +  1 =  1
         *  1 +  0 =  1
         * -1 + -1 = -1
         * -1 +  0 = -1
         * -1 +  1 =  0
         *  0 +  0 =  0
         *
         */

        vector<int> parState;

        for (unsigned int i = 0; i < child_l.size(); ++i) {
            if (child_l[i] == child_r[i]) {
                parState.emplace_back(child_l[i]);
            }
            else {
                parState.emplace_back(child_l[i] + child_r[i]);
            }
        }

        return parState;
    }


    /**
    * \brief Collects all children of a given node and sets all calculated states for them.
    * @param r1Node Pointer to the node of a PhylogeneticTree (choose root to infer states for the whole tree).
    * @return The inferred states for the given node as a vector<int>. (All states of its childrens in the PhylogeneticTree have been set)
    * \relates TreeNodePhylo
    */
    vector<int>
    calChildren(TreeNodePhylo <vector<int> > *r1Node) {
        if (r1Node->isLeaf()) {
            return r1Node->data;
        }

        auto cl = r1Node->child(0);
        auto cr = r1Node->child(1);

        cl->data = calChildren(cl);
        cr->data = calChildren(cr);

        return parentalState(cl->data, cr->data);

    }


    /**
    * \brief Calculate in a second iteration over the tree all unknown states.
    * @param r2Node Root node of a phylogenetic tree to iterate over.
    * \relates TreeNodePhylo
    */
    void calUncertainStates(TreeNodePhylo <vector<int>> *r2Node) {
        if (r2Node->isLeaf()) //node is leaf
        {
            return;
        }
        else if (r2Node->parent() == nullptr) //node is root
        {
            //set all uncertain states to 1
            for (int &state: r2Node->data) {
                if (state == 0) {
                    state = 1;
                }
            }
        }
        else {
            //set all instances of 0 to parental state
            for (unsigned int n = 0; n < r2Node->data.size(); ++n) {
                if (r2Node->data[n] == 0) {
                    r2Node->data[n] = r2Node->parent()->data[n];
                }
            }

        }

        auto l_child = r2Node->child(0);
        auto r_child = r2Node->child(1);
        calUncertainStates(l_child);
        calUncertainStates(r_child);

    }

    /**
     * \brief Infers all states of inner nodes of a phylogenetic tree according to a fitch parsimony.
     * @param phyTree The input tree. Has to be a BioSeqDataLib::PhylogeneticTree with a vector<int> in the data field.
     * \relates PhylogeneticTree
     */
    void fitch(PhylogeneticTree <vector<int>> &phyTree) {
        TreeNodePhylo <vector<int>> *rootTreeNode = &phyTree.root();

        rootTreeNode->data = calChildren(rootTreeNode);

        calUncertainStates(rootTreeNode);

    }
}
