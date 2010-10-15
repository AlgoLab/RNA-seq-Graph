#ifndef BUILD_CHAINS_H
#define BUILD_CHAINS_H

#include "read_fasta.h"

//RNA-seqs Print Functions
void print_spliced(tables);
void print_unspliced(tables);
void print_hash_table(tables, char);

//Chains Building (half or specified overlap)
void build_unspliced_chains(tables);
void build_unspliced_chains(tables, int);

//Chains Printing
void print_unspliced_chains(tables);
void print_unspliced_chains(tables, int);

//Chains Merging
string merge_chains(string, int, string);
map<unsigned long long, string> merge_unspliced_chains(tables);

#endif
