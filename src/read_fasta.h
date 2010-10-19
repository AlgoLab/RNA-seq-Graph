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
};

typedef std::map<unsigned long long, element_table> Map;


//Hash tables
struct tables{
    Map left_map;
    Map right_map;
};

//Convert DNA sequence into binary
string binary_conversion(string);

//Encode binary sequence into a number
unsigned long long fingerprint(string);

//Parse Fasta Information
table_entry* parse_fasta(String<Dna5>, string);

//Add entries in the Hash Table
void add_entry(Map &, table_entry*, char);

//Read Fasta File
void read_fasta(char*, tables&);

#endif
