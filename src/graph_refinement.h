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
