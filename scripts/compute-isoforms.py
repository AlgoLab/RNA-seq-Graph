#!/usr/bin/env python3
####
#
#
#                              RNA-Graph
#
# Description
#
# Copyright (C) 2011  Gianluca Della Vedova
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
####

import networkx as nx
import pprint as pp


def enumerate_all_paths(graph, origin, destination):
    '''
    enumerate all paths between an origin and a destination
    http://codeidol.com/python/python3/Data-Structures/Graph-Searching/
    implementation from http://groups.google.com/group/networkx-discuss/browse_thread/thread/3f9c08f3dab15077/3b1cd30b37ffa247
    '''

    paths = []
    stack = ([origin], [])
    while stack:
        front, stack = stack
        end = front[-1]
        if end == destination:
            paths.append(front)
        else:
            for successor in graph.successors(end):
                if successor not in front:
                    stack = (front + [successor]), stack

    return paths


def splice_fasta(string):
    rows = []
    while string:
        (head, string) = (string[:60], string[60:])
        rows.append(head)
    return rows


if __name__ == '__main__':
    G=nx.DiGraph()
    G=nx.read_graphml("RNA-seq-graph.graphml")
    # calculate the paths between all OD pairs
    isoforms = []
    for cc in nx.weakly_connected_components(G):
        for origin in cc:
            for destination in cc:
                if origin != destination:
                    paths = enumerate_all_paths(G, origin, destination)
                    for p in paths:
                        isoforms.append(''.join([nx.get_node_attributes(G,'seq')[v] for v in p]))

    for (i,isoform) in enumerate(isoforms):
        print("> isoform {}".format(i))
        [print(r) for r in splice_fasta(isoform)]
