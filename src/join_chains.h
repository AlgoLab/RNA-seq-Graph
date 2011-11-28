/**
* RNA-seq-Graph
* Method for reconstructing the Isoform Graph of a gene from RNA-seq data, without the genome information
*
* Copyright (C) 2011 Stefano Beretta <ste.beretta(-at-)gmail.com>
*
* Distributed under the terms of the GNU Affero General Public License (AGPL)
*
*
* This file is part of RNA-seq-Graph.
*
* RNA-seq-Graph is free software: you can redistribute it and/or modify
* it under the terms of the GNU Affero General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* RNA-seq-Graph is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU Affero General Public License for more details.
*
* A full copy of the GNU Affero General Public License is reported into
* the file COPYING. However more informations ca be found at
* <http://www.gnu.org/licenses/>.
**/
#ifndef JOIN_CHAINS_H
#define JOIN_CHAINS_H

#include <queue>
#include "read_fasta.h"
#include "graph_refinement.h"

struct links_pair{
    unsigned long long D_chain;
    unsigned long long A_chain;
};

struct small_frag{
    links_pair frag_links;
    CharString frag;
    ::std::vector<links_pair> other_links;
};

int get_left_linked_read(string, tables&, int);
int get_right_linked_read(string, tables&, int);

::std::map<unsigned long long, unsigned long long> chain_back_merging(map<unsigned long long, string>&, int);
::std::map<unsigned long long, unsigned long long> chains_unify(map<unsigned long long, string>&, unsigned int);

void link_fragment_chains(tables&, ::std::map<unsigned long long, string> &, int, char*);

void print_graph(::std::vector<table_entry*>, ::std::map<unsigned long long, string>, 
                 ::std::map<unsigned long long, unsigned long long>, char*);

void check_cutted_frags(CharString, ::std::vector<table_entry*> &, 
                        map<unsigned long long, string> &, unsigned int);
#endif
