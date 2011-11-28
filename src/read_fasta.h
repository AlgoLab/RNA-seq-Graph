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
#ifndef READ_FASTA_H
#define READ_FASTA_H

#include <map>
#include <string>

#include <seqan/sequence.h>
#include <seqan/file.h>

#include "table_entry.h"

using namespace seqan;

//Elements of the hash tables
struct element_table{
    table_entry* p;
    bool unspliced;
    bool half_spliced;
};

typedef ::std::map<unsigned long long, element_table> hash_map;


//Hash tables
struct tables{
    hash_map left_map;
    hash_map right_map;
};

//Convert DNA sequence into binary
unsigned long long fingerprint(const string&);

//Parse Fasta Information
table_entry* parse_fasta(String<Dna5>, string);

//Build table_entry
table_entry* build_entry(String<Dna5>);

//Add entries in the Hash Table
void add_entry(hash_map &, table_entry*, char);

//Read Fasta File
int read_fasta(char*, tables&);

#endif
