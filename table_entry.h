#ifndef TABLE_ENTRY_H
#define TABLE_ENTRY_H

#include "RNA_seq.h"

class table_entry{
    private:
    RNA_seq* short_read;
    //Table List
    table_entry* next;
    table_entry* prev;
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
    void set_next(table_entry*);
    void set_prev(table_entry*);
    void set_left_fingerprint(unsigned long long);
    void set_right_fingerprint(unsigned long long);
    void set_junction_offset(int);
    void set_chain_next(table_entry*);
    void set_chain_prev(table_entry*);
    //Get Methods
    RNA_seq* get_short_read();
    table_entry* get_next();
    table_entry* get_prev();
    unsigned long long get_left_fingerprint();
    unsigned long long get_right_fingerprint();
    int get_junction_offset();
    table_entry* get_chain_next();
    table_entry* get_chain_prev();
long get_frequency();
    //Increase / Decrease Sequence Frequency
    void increase_freq();
    void decrease_freq();
};

#endif
