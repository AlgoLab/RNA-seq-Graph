#include "read_fasta.h"

struct delta_link{
    String<Dna5> seq;
    unsigned long long l_fingerprint;
    unsigned long long r_fingerprint;
    vector<unsigned long long> D_delta;
    vector<unsigned long long> A_delta;
};

void get_left_linked_read(string, ::std::vector<delta_link> &, int);
void get_right_linked_read(string, ::std::vector<delta_link> &, int);

::std::vector<delta_link> link_fragment_chains(const tables&, ::std::map<unsigned long long, string>);

void print_graph(::std::vector<delta_link>);
