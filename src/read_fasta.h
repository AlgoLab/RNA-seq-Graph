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
string binary_conversion(string);

//Encode binary sequence into a number
unsigned long long fingerprint(string);

//Parse Fasta Information
table_entry* parse_fasta(String<Dna5>, string);

//Add entries in the Hash Table
void add_entry(hash_map &, table_entry*, char);

//Read Fasta File
int read_fasta(char*, tables&);

#endif
