#include "table_entry.h"

//Constructors
table_entry::table_entry(String<Dna5> seq,
                         unsigned long long left_f, unsigned long long right_f){
    short_read = new RNA_seq(seq);
    this->next = NULL;
    this->prev = NULL;
    this->left_fingerprint = left_f;
    this->right_fingerprint = right_f;
    this->junction_offset = 0;
    this->chain_next = NULL;
    this->chain_prev = NULL;
    this->frequency = 1;
}

table_entry::table_entry(String<Dna5> seq, string source, string gene_id, 
                         unsigned long long left_f, unsigned long long right_f){
    short_read = new RNA_seq(seq, source, gene_id);
    this->next = NULL;
    this->prev = NULL;
    this->left_fingerprint = left_f;
    this->right_fingerprint = right_f;
    this->junction_offset = 0;
    this->chain_next = NULL;
    this->chain_prev = NULL;
    this->frequency = 1;
}

table_entry::table_entry(String<Dna5> seq, string source, string gene_id, 
                         int gene_strand,
                         unsigned long long left_f, unsigned long long right_f){
    short_read = new RNA_seq(seq, source, gene_id, gene_strand);
    this->next = NULL;
    this->prev = NULL;
    this->left_fingerprint = left_f;
    this->right_fingerprint = right_f;
    this->junction_offset = 0;
    this->chain_next = NULL;
    this->chain_prev = NULL;
    this->frequency = 1;
}

table_entry::table_entry(String<Dna5> seq, string source, string gene_id, 
                         int gene_strand, int transcript_id, long offset, int clone_end,
                         unsigned long long left_f, unsigned long long right_f){
    short_read = new RNA_seq(seq, source, gene_id, gene_strand, transcript_id, offset, clone_end);
    this->next = NULL;
    this->prev = NULL;
    this->left_fingerprint = left_f;
    this->right_fingerprint = right_f;
    this->junction_offset = 0;
    this->chain_next = NULL;
    this->chain_prev = NULL;
    this->frequency = 1;
}

table_entry::~table_entry(){
    if (short_read != NULL)
        delete short_read;
}

table_entry::table_entry(const table_entry& rhs){
    short_read = new RNA_seq(*(rhs.short_read));
    prev = rhs.prev;
    next = rhs.next;
    left_fingerprint = rhs.left_fingerprint;
    right_fingerprint = rhs.right_fingerprint;
    junction_offset = rhs.junction_offset;
    chain_next = rhs.chain_next;
    chain_prev = rhs.chain_prev;
    frequency = rhs.frequency;
}

table_entry& table_entry::operator=(const table_entry& rhs){
    if(&rhs == this)
        return *this;
    if (short_read != NULL)
        delete short_read;
    short_read = new RNA_seq(*(rhs.short_read));
    prev = rhs.prev;
    next = rhs.next;
    left_fingerprint = rhs.left_fingerprint;
    right_fingerprint = rhs.right_fingerprint;
    junction_offset = rhs.junction_offset;
    chain_next = rhs.chain_next;
    chain_prev = rhs.chain_prev;
    frequency = rhs.frequency;
    return *this;
}

//Set Methods
void table_entry::set_next(table_entry* next){
    this->next = next;
}

void table_entry::set_prev(table_entry* prev){
    this->prev = prev;
}

void table_entry::set_left_fingerprint(unsigned long long left_f){
    this->left_fingerprint = left_f;
}

void table_entry::set_right_fingerprint(unsigned long long right_f){
    this->right_fingerprint = right_f;
}

void table_entry::set_junction_offset(int junction_offset){
    this->junction_offset = junction_offset;
}

void table_entry::set_chain_next(table_entry* chain_next){
    this->chain_next = chain_next;
}

void table_entry::set_chain_prev(table_entry* chain_prev){
    this->chain_prev = chain_prev;
}

//Get Methods
RNA_seq* table_entry::get_short_read(){
    return short_read;
}

table_entry* table_entry::get_next(){
    return next;
}

table_entry* table_entry::get_prev(){
    return prev;
}

unsigned long long table_entry::get_left_fingerprint(){
    return left_fingerprint;
}

unsigned long long table_entry::get_right_fingerprint(){
    return right_fingerprint;
}
int table_entry::get_junction_offset(){
    return junction_offset;
}

table_entry* table_entry::get_chain_next(){
    return chain_next;
}

table_entry* table_entry::get_chain_prev(){
    return chain_prev;
}

long table_entry::get_frequency(){
    return frequency;
}

///Increase / Decrease Sequence Frequency

void table_entry::increase_freq(){
    this->frequency++;
}

void table_entry::decrease_freq(){
    this->frequency--;
}

