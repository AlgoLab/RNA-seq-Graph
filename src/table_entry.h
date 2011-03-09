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
    int junction_offset;
    //Chain List
    table_entry* chain_next;
    table_entry* chain_prev;
    //Sequence frequency
    long frequency;
    
    //Chains Linked
    vector<unsigned long long> D_delta;
    vector<unsigned long long> A_delta;

    public:
    //Constructors
    table_entry(String<Dna5>, unsigned long long, unsigned long long);
    table_entry(String<Dna5>, string, string, unsigned long long, unsigned long long);
    table_entry(String<Dna5>, string, string, int, unsigned long long, unsigned long long);
    table_entry(String<Dna5>, string, string, int, int, long, int, unsigned long long, unsigned long long);
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
    void set_junction_offset(int);
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
    int get_junction_offset() const;
    table_entry* get_chain_next() const;
    table_entry* get_chain_prev() const;
    long get_frequency() const;

    //Increase / Decrease Sequence Frequency
    void increase_freq();
    void decrease_freq();

    //Linking chains methods
    void push_A_link(unsigned long long);
    int size_A_link()const;
    unsigned long long at_A_link(int)const;

    void push_D_link(unsigned long long);
    int size_D_link()const;
    unsigned long long at_D_link(int)const;
};

#endif
