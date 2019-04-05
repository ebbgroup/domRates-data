#!/usr/bin/env python
# -*- coding: utf-8 -*-

# visualization of domain rearrangement events on a phylogenetic tree

#    This file is part of DomRates.
#
#    DomRates is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    DomRates is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with DomRates.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import print_function
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, PieChartFace, TextFace, CircleFace
import os, sys
import matplotlib.pyplot as pl
import argparse

__author__ = "Elias Dohmen" 
__version__ = "0.8"
__email__ = "e.dohmen@uni-muenster.de"
__institute__ = "IEB Münster"


def layout_idtree(node):

    nstyle = NodeStyle()
    nstyle["size"] = 1
    node.set_style(nstyle)
    if node.is_leaf():
        # Add node name to leaf nodes
        N = AttrFace("name", fsize=14, fgcolor="black")
        faces.add_face_to_node(N, node, 0, position='aligned')

    if not node.is_root():
        # Creates a sphere face whose size is proportional to node's
        # feature "weight"


        TN = TextFace(str(node.ID), fsize=12)
        faces.add_face_to_node(TN, node, 1, position="branch-right")



def layout_gen_events(node):

    global scal
    global eve

    nstyle = NodeStyle()
    nstyle["size"] = 1
    node.set_style(nstyle)
    if node.is_leaf():
        # Add node name to leaf nodes
        N = AttrFace("name", fsize=14, fgcolor="black")
        faces.add_face_to_node(N, node, 0, position='aligned')

    if not node.is_root():
        # Creates a sphere face whose size is proportional to node's
        # feature "weight"
        tot = float(node.fusion + node.fission + node.termLoss + node.termGain + node.singLoss + node.singGain)
        if tot > 0 :
            perc_fus = 100 * (node.fusion / tot)
            perc_fis = 100 * (node.fission / tot)
            perc_tL = 100 * (node.termLoss / tot)
            perc_tG = 100 * (node.termGain / tot)
            perc_sL = 100 * (node.singLoss / tot)
            perc_sG = 100 * (node.singGain / tot)

            if scal > 0:
                if eve == 0:
                    sf = tot*scal
                else:
                    sf = [node.fusion, node.fission, node.termLoss, node.termGain, node.singLoss, node.singGain][eve-1]*scal
            else:
                sf = 50

            if eve != 0:
                ev_vec = [100]
                col_vec = []
                col_vec.append(["DimGray", "DeepPink", "YellowGreen", "DarkBlue", "Chocolate", "DeepSkyBlue"][eve-1])
            else:
                ev_vec = [perc_fus, perc_fis, perc_tL, perc_tG, perc_sL, perc_sG]
                col_vec = ["DimGray", "DeepPink", "YellowGreen", "DarkBlue", "Chocolate", "DeepSkyBlue"]

            C = PieChartFace( ev_vec, height=sf, width=sf, colors=col_vec)

            if eve == 0:
                TN = TextFace(str(int(tot)), fsize=12)
            else:
                e = [node.fusion, node.fission, node.termLoss, node.termGain, node.singLoss, node.singGain][eve-1]
                TN = TextFace(str(e),fsize=12)

            faces.add_face_to_node(TN, node, 1, position="branch-right")

            # place as a float face over the tree
            faces.add_face_to_node(C, node, 0, position="float")

def main( ) :
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", action="store", dest="treef", help="The tree file")
    parser.add_argument("-s", action="store", dest="statf", help="The statistics file")
    parser.add_argument("-c", action="store", default=0, type= float, dest="scaling", help="If this parameter is set "
                                                                                       "to more than 0, "
                                                                                       "the size of the pie charts "
                                                                                       "correlate with the total "
                                                                                       "number of events at a node "
                                                                                       "(and are scaled by the factor "
                                                                                       "given as a float).")
    parser.add_argument("-e", action="store", type=str, default= "all", dest="event", help="If an event type is specified, just this event type is visualized. Per default all event types are shown on the tree.\n"
                                                                                       "all\n"
                                                                                       "fusions\n"
                                                                                       "fissions\n"
                                                                                       "termLosses\n"
                                                                                       "termEmergences\n"
                                                                                       "singleDomainLosses\n"
                                                                                       "singleDomainEmergences")
    parser.add_argument("-p", action="store", dest="treeshape", default="r", choices=["c", "r"],
        help="shape of the tree, circle (c) or tree format (r)")
    parser.add_argument("-o", action="store", dest="outputname")
    parser.add_argument("-y", action="store", type=str, dest="NodeIDtreeName", default=None, help="Name for output file that shows a tree with all node IDs.")
    parser.add_argument("-l", dest="short_legend", help="Writes the full legend for all events in two levels for short trees", action="store_true")
	
    params = parser.parse_args()

    if (params.event not in ("all", "fusions", "fissions", "termLosses", "termEmergences", "singleDomainLosses", "singleDomainEmergences")):
        print("Error: Please specify a valid event type. For a list of possible options use the --help parameter.")
        sys.exit(1)

    if params.NodeIDtreeName != None:
        id_tree = Tree(params.treef, format=0)

        coun = 0
        for node in id_tree.traverse('preorder'):
            node.add_features(ID=coun)
            coun += 1


        # Create empty TreeStyle
        ts = TreeStyle()

        # Set custom layout function
        ts.layout_fn = layout_idtree
        # Draw tree
        ts.mode = params.treeshape
        ts.complete_branch_lines_when_necessary = True
        ts.extra_branch_line_type = 0
        ts.extra_branch_line_color = "black"
        # ts.optimal_scale_level ="full"
        ts.branch_vertical_margin = 40
        ts.scale = 100

        # We will add node names manually
        ts.show_leaf_name = False
        ts.draw_guiding_lines = True

        if (params.NodeIDtreeName.endswith(".pdf")):
            pathout = params.NodeIDtreeName
        else:
            pathout = params.NodeIDtreeName + ".pdf"
        id_tree.render(pathout, dpi=1200, tree_style=ts)
        pl.close()

    else:

        tree = Tree(params.treef, format=0)

        # Read statistics file
        node_stat_dict = {}
        with open(params.statf, "r") as sf:
            for line in sf:
                # Stop the loop at the second part of statistics file
                if line.startswith("# Number of events per domain.") or line.startswith("# Events per domain arrangement for last common ancestor"):
                    break
                if line[0] not in('#', '\n'):
                    vecline = line.strip().split()
                    id = vecline.pop(0)
                    stats = [int(i) for i in vecline]
                    node_stat_dict[int(id)] = stats

        # determine max. number of events per node for scaling
        fus_max = 0
        fis_max = 0
        termLoss_max = 0
        termGain_max = 0
        singLoss_max = 0
        singGain_max = 0
        tot_max = 0

        # Assign rearrangement events to leaves
        c = 0
        for node in tree.traverse('preorder'):
            node.add_features(fusion=node_stat_dict[c][0])
            if (node_stat_dict[c][0] > fus_max):
                fus_max = node_stat_dict[c][0]
            node.add_features(fission=node_stat_dict[c][1])
            if (node_stat_dict[c][1] > fis_max):
                fis_max = node_stat_dict[c][1]
            node.add_features(termLoss=node_stat_dict[c][2])
            if (node_stat_dict[c][2] > termLoss_max):
                termLoss_max = node_stat_dict[c][2]
            node.add_features(termGain=node_stat_dict[c][3])
            if (node_stat_dict[c][3] > termGain_max):
                termGain_max = node_stat_dict[c][3]
            node.add_features(singLoss=node_stat_dict[c][4])
            if (node_stat_dict[c][4] > singLoss_max):
                singLoss_max = node_stat_dict[c][4]
            node.add_features(singGain=node_stat_dict[c][5])
            if (node_stat_dict[c][5] > singGain_max):
                singGain_max = node_stat_dict[c][5]
            if (sum(node_stat_dict[c]) > tot_max):
                tot_max = sum(node_stat_dict[c])
            c += 1


        global scal
        scal = params.scaling

        global eve
        event_options = {
            "all": 0,
            "fusions": 1,
            "fissions": 2,
            "termLosses": 3,
            "termEmergences": 4,
            "singleDomainLosses": 5,
            "singleDomainEmergences": 6
        }
        eve = event_options[params.event]


        # Create empty TreeStyle
        ts = TreeStyle()

        # Set custom layout function
        ts.layout_fn = layout_gen_events
        # Draw tree
        ts.mode = params.treeshape
        ts.complete_branch_lines_when_necessary = True
        ts.extra_branch_line_type = 0
        ts.extra_branch_line_color = "black"
        #ts.optimal_scale_level ="full"
        ts.branch_vertical_margin = 40
        ts.scale =  100

        # We will add node names manually
        ts.show_leaf_name = False

        # legend creation
        if (params.event == "all"):
            ts.legend.add_face(CircleFace(10, "DimGray"), column=0)
            ts.legend.add_face(TextFace(" Fusion     ", fsize=16, fgcolor='DimGray'), column=1)
            ts.legend.add_face(CircleFace(10, "DeepPink"), column=2)
            ts.legend.add_face(TextFace(' Fission     ', fsize=16, fgcolor='DeepPink'), column=3)
            ts.legend.add_face(CircleFace(10, "YellowGreen"), column=4)
            ts.legend.add_face(TextFace(' Terminal Loss     ', fsize=16, fgcolor='YellowGreen'), column=5)
            if params.short_legend:
                ts.legend.add_face(CircleFace(10, "DarkBlue"), column=0)
                ts.legend.add_face(TextFace(' Terminal Emergence     ', fsize=16, fgcolor='DarkBlue'), column=1)
                ts.legend.add_face(CircleFace(10, "Chocolate"), column=2)
                ts.legend.add_face(TextFace(' Single Domain Loss     ', fsize=16, fgcolor='Chocolate'), column=3)
                ts.legend.add_face(CircleFace(10, "DeepSkyBlue"), column=4)
                ts.legend.add_face(TextFace(' Single Domain Emergence     ', fsize=16, fgcolor='DeepSkyBlue'), column=5)
            else:
                ts.legend.add_face(CircleFace(10, "DarkBlue"), column=6)
                ts.legend.add_face(TextFace(' Terminal Emergence     ', fsize=16, fgcolor='DarkBlue'), column=7)
                ts.legend.add_face(CircleFace(10, "Chocolate"), column=8)
                ts.legend.add_face(TextFace(' Single Domain Loss     ', fsize=16, fgcolor='Chocolate'), column=9)
                ts.legend.add_face(CircleFace(10, "DeepSkyBlue"), column=10)
                ts.legend.add_face(TextFace(' Single Domain Emergence     ', fsize=16, fgcolor='DeepSkyBlue'), column=11)
        elif (params.event == "fusions"):
            ts.legend.add_face(CircleFace(10, "DimGray"), column=0)
            ts.legend.add_face(TextFace(" Fusion     ", fsize=16, fgcolor='DimGray'), column=1)
        elif (params.event == "fissions"):
            ts.legend.add_face(CircleFace(10, "DeepPink"), column=0)
            ts.legend.add_face(TextFace(' Fission     ', fsize=16, fgcolor='DeepPink'), column=1)
        elif (params.event == "termLosses"):
            ts.legend.add_face(CircleFace(10, "YellowGreen"), column=0)
            ts.legend.add_face(TextFace(' Terminal Loss     ', fsize=16, fgcolor='YellowGreen'), column=1)
        elif (params.event == "termEmergences"):
            ts.legend.add_face(CircleFace(10, "DarkBlue"), column=0)
            ts.legend.add_face(TextFace(' Terminal Emergence     ', fsize=16, fgcolor='DarkBlue'), column=1)
        elif (params.event == "singleDomainLosses"):
            ts.legend.add_face(CircleFace(10, "Chocolate"), column=0)
            ts.legend.add_face(TextFace(' Single Domain Loss     ', fsize=16, fgcolor='Chocolate'), column=1)
        elif (params.event == "singleDomainEmergences"):
            ts.legend.add_face(CircleFace(10, "DeepSkyBlue"), column=0)
            ts.legend.add_face(TextFace(' Single Domain Emergence     ', fsize=16, fgcolor='DeepSkyBlue'), column=1)

        ts.legend_position=1
        ts.draw_guiding_lines = True

        if (params.outputname.endswith(".pdf")):
            pathout = params.outputname
        else:
            pathout = params.outputname + ".pdf"
        tree.render(pathout, dpi=1200, tree_style=ts)
        pl.close( )

    sys.exit( 0 )
    
if __name__ == "__main__" :
    main( )
