/**
* RNA-seq-Graph
* Method for reconstructing the Isoform Graph of a gene from RNA-seq data, without the genome information
*
* Copyright (C) 2011 Stefano Beretta <ste.beretta(-at-)gmail.com>
*
* Distributed under the terms of the GNU General Public License (GPL)
*
*
* This file is part of RNA-seq-Graph.
*
* RNA-seq-Graph is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* RNA-seq-Graph is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* A full copy of the GNU General Public License is reported into
* the file COPYING. However more informations ca be found at
* <http://www.gnu.org/licenses/>.
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
                 ::std::map<unsigned long long, unsigned long long>&);

void linking_refinement(::std::vector<table_entry*> &, map<unsigned long long, string> &, unsigned int, 
                        ::std::map<unsigned long long, unsigned long long>&);

void check_overlapping_nodes(::std::vector<table_entry*> &, map<unsigned long long, string> &, int,
                             ::std::map<unsigned long long, unsigned long long>&, unsigned int, int);

#endif
