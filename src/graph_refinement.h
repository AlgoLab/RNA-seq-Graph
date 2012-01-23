/**
* RNA-seq-Graph
* Method for reconstructing the Isoform Graph of a gene from RNA-seq data
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
* You should have received a copy of the GNU Affero General Public License
* along with RNA-seq-Graph.  If not, see <http://www.gnu.org/licenses/>.
**/
#ifndef GRAPH_REFINEMENT_H
#define GRAPH_REFINEMENT_H

#include "join_chains.h"

void add_linking_reads(::std::vector<table_entry*> &, const map<unsigned long long, string>, unsigned int);

unsigned int overlappedStringLength(string, string);
int* computeBackTrackTable(string);
void small_blocks(::std::vector<table_entry*> &, map<unsigned long long, string> &, unsigned int,
                  ::std::map<unsigned long long, unsigned long long>&);
void tiny_blocks(::std::vector<table_entry*> &, map<unsigned long long, string> &, int,
                 ::std::map<unsigned long long, unsigned long long>&, unsigned int);

void linking_refinement(::std::vector<table_entry*> &, map<unsigned long long, string> &, unsigned int, 
                        ::std::map<unsigned long long, unsigned long long>&);

void check_overlapping_nodes(::std::vector<table_entry*> &, map<unsigned long long, string> &, int,
                             ::std::map<unsigned long long, unsigned long long>&, unsigned int, int);

#endif
