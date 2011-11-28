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
/*******************************/
/* Class that models the entry */
/* of the hash table           */
/*******************************/
#ifndef TABLE_ENTRY_H
#define TABLE_ENTRY_H

#include "RNA_seq.h"

class table_entry{
    private:
    RNA_seq* short_read;
    //Table List
    table_entry* r_next;
    table_entry* r_prev;
    table_entry* l_next;
    table_entry* l_prev;
    //Fingerprints
    unsigned long long left_fingerprint;
    unsigned long long right_fingerprint;
    //Fragment Junction Offset
    //int junction_offset;
    //Chain List
    table_entry* chain_next;
    table_entry* chain_prev;
    //Sequence frequency
    long frequency;

    public:
    //Chains Linked
    vector<unsigned long long> D_delta;
    vector<unsigned long long> A_delta;
    //Constructors
    table_entry(String<Dna5>, unsigned long long, unsigned long long);
    //table_entry(String<Dna5>, string, string, unsigned long long, unsigned long long);
    //table_entry(String<Dna5>, string, string, int, unsigned long long, unsigned long long);
    //table_entry(String<Dna5>, string, string, int, int, long, int, unsigned long long, unsigned long long);
    ~table_entry();

    //Copy constructor
    table_entry(const table_entry&);
    table_entry& operator=(const table_entry&);

    //Set Methods
    void set_l_next(table_entry*);
    void set_l_prev(table_entry*);
    void set_r_next(table_entry*);
    void set_r_prev(table_entry*);
    void set_left_fingerprint(unsigned long long);
    void set_right_fingerprint(unsigned long long);
    //void set_junction_offset(int);
    void set_chain_next(table_entry*);
    void set_chain_prev(table_entry*);

    //Get Methods
    RNA_seq* get_short_read() const;
    table_entry* get_l_next() const;
    table_entry* get_l_prev() const;
    table_entry* get_r_next() const;
    table_entry* get_r_prev() const;
    unsigned long long get_left_fingerprint() const;
    unsigned long long get_right_fingerprint() const;
    //int get_junction_offset() const;
    table_entry* get_chain_next() const;
    table_entry* get_chain_prev() const;
    long get_frequency() const;

    //Increase / Decrease Sequence Frequency
    void increase_freq();
    void decrease_freq();

    //Linking chains methods
    void push_A_link(unsigned long long);
    int size_A_link()const;
    unsigned long long & at_A_link(int);

    void push_D_link(unsigned long long);
    int size_D_link()const;
    unsigned long long & at_D_link(int);

};

#endif
