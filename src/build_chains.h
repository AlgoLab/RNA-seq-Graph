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