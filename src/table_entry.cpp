#include "table_entry.h"

//Constructors
table_entry::table_entry(String<Dna5> seq,
                         unsigned long long left_f, unsigned long long right_f){
    short_read = new RNA_seq(seq);

    this->l_next = NULL;
    this->l_prev = NULL;
    this->r_next = NULL;
    this->r_prev = NULL;

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
 
    this->l_next = NULL;
    this->l_prev = NULL;
    this->r_next = NULL;
    this->r_prev = NULL;

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
 
    this->l_next = NULL;
    this->l_prev = NULL;
    this->r_next = NULL;
    this->r_prev = NULL;

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

    this->l_next = NULL;
    this->l_prev = NULL;
    this->r_next = NULL;
    this->r_prev = NULL;

    this->left_fingerprint = left_f;
    this->right_fingerprint = right_f;
    this->junction_offset = 0;
    this->chain_next = NULL;
    this->chain_prev = NULL;
    this->frequency = 1;
}

//Distructor
table_entry::~table_entry(){
    if (short_read != NULL){
        delete short_read;
    }
    //if(next != NULL){
    //    delete next;
    //}
}

//Copy Constructor
//in which the pointers next and prev are
//setted at NULL
table_entry::table_entry(const table_entry& rhs){
    short_read = new RNA_seq(*(rhs.short_read));

    l_next = rhs.l_next;
    l_prev = rhs.l_prev;
    r_next = rhs.r_next;
    r_prev = rhs.r_prev;

    left_fingerprint = rhs.left_fingerprint;
    right_fingerprint = rhs.right_fingerprint;
    junction_offset = rhs.junction_offset;
    chain_next = rhs.chain_next;
    chain_prev = rhs.chain_prev;
    frequency = rhs.frequency;
}

//Overloading =
table_entry& table_entry::operator=(const table_entry& rhs){
    if(&rhs == this)
        return *this;
    if (short_read != NULL){
        delete short_read;
    }
    short_read = new RNA_seq(*(rhs.short_read));

    l_next = rhs.l_next;
    l_prev = rhs.l_prev;
    r_next = rhs.r_next;
    r_prev = rhs.r_prev;

    left_fingerprint = rhs.left_fingerprint;
    right_fingerprint = rhs.right_fingerprint;
    junction_offset = rhs.junction_offset;
    chain_next = rhs.chain_next;
    chain_prev = rhs.chain_prev;
    frequency = rhs.frequency;
    return *this;
}

//Set Methods
void table_entry::set_l_next(table_entry* next){
    this->l_next = next;
}

void table_entry::set_l_prev(table_entry* prev){
    this->l_prev = prev;
}

void table_entry::set_r_next(table_entry* next){
    this->r_next = next;
}

void table_entry::set_r_prev(table_entry* prev){
    this->r_prev = prev;
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
RNA_seq* table_entry::get_short_read() const{
    return short_read;
}

table_entry* table_entry::get_l_next() const{
    return l_next;
}

table_entry* table_entry::get_l_prev() const{
    return l_prev;
}

table_entry* table_entry::get_r_next() const{
    return r_next;
}

table_entry* table_entry::get_r_prev() const{
    return r_prev;
}

unsigned long long table_entry::get_left_fingerprint() const{
    return left_fingerprint;
}

unsigned long long table_entry::get_right_fingerprint() const{
    return right_fingerprint;
}
int table_entry::get_junction_offset() const{
    return junction_offset;
}

table_entry* table_entry::get_chain_next() const{
    return chain_next;
}

table_entry* table_entry::get_chain_prev() const{
    return chain_prev;
}

long table_entry::get_frequency() const{
    return frequency;
}

///Increase / Decrease Sequence Frequency

void table_entry::increase_freq(){
    this->frequency++;
}

void table_entry::decrease_freq(){
    this->frequency--;
}

//Linking Chains Methods
void table_entry::push_A_link(unsigned long long el){
    A_delta.push_back(el);
}

int table_entry::size_A_link()const{
    return A_delta.size();
}

unsigned long long table_entry::at_A_link(int pos)const{
    return A_delta[pos];
}

void table_entry::push_D_link(unsigned long long el){
    D_delta.push_back(el);
}

int table_entry::size_D_link()const{
    return D_delta.size();
}

unsigned long long table_entry::at_D_link(int pos)const{
    return D_delta[pos];
}

