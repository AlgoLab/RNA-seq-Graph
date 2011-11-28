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
#ifndef BUILD_CHAINS_H
#define BUILD_CHAINS_H

#include "read_fasta.h"

//RNA-seqs Print Functions
void print_spliced(const tables&);
void print_half_spliced(const tables&);
void print_unspliced(const tables&);
void print_hash_table(const tables&, char);

//Chains Building (half or specified overlap)
void build_unspliced_chains(tables&);
void build_unspliced_chains(tables&, int);

//Chains Printing
void print_unspliced_chains(const tables&);
void print_unspliced_chains(const tables&, int);

//Chains Merging
void merge_unspliced_chains(tables&, map<unsigned long long, string>&);
void print_merged_chains(::std::map<unsigned long long, string> &);
#endif
