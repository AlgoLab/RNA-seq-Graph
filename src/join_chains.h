#ifndef JOIN_CHAINS_H
#define JOIN_CHAINS_H

#include <queue>
#include "read_fasta.h"
#include "graph_refinement.h"

struct links_pair{
    unsigned long long D_chain;
    unsigned long long A_chain;
};

struct small_frag{
    links_pair frag_links;
    CharString frag;
    ::std::vector<links_pair> other_links;
};

int get_left_linked_read(string, tables&, int);
int get_right_linked_read(string, tables&, int);

::std::map<unsigned long long, unsigned long long> chain_back_merging(map<unsigned long long, string>&, int);
::std::map<unsigned long long, unsigned long long> chains_unify(map<unsigned long long, string>&, unsigned int);

void link_fragment_chains(tables&, ::std::map<unsigned long long, string> &, int, char*);

void print_graph(::std::vector<table_entry*>, ::std::map<unsigned long long, string>, 
                 ::std::map<unsigned long long, unsigned long long>, char*);

void check_cutted_frags(CharString, ::std::vector<table_entry*> &, 
                        map<unsigned long long, string> &, unsigned int);
#endif
