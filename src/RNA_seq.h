#ifndef RNA_SEQ_H
#define RNA_SEQ_H

#include <string>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;
using namespace std;

class RNA_seq{
    private:
    String<Dna5> sequence;  
    string source;
    string gene_id;
    int gene_strand;
    int transcript_id;
    //Offset starting from the first digit of the transcript
    long offset;
    int clone_end;
    
    public:
    //Constructors
    RNA_seq(String<Dna5>);
    RNA_seq(String<Dna5>, string, string);
    RNA_seq(String<Dna5>, string, string, int);
    RNA_seq(String<Dna5>, string, string, int, int, long, int);
    //Copy Constructor    
    RNA_seq(const RNA_seq&);
    //Overloading =
    RNA_seq& operator=(const RNA_seq&);
    //Set Methods
    void set_RNA_seq_sequence(String<Dna5>);
    void set_RNA_seq_source(string);
    void set_RNA_seq_gene_id(string);
    void set_RNA_seq_gene_strand(int);
    void set_RNA_seq_transcript_id(int);
    void set_RNA_seq_offset(long);
    void set_RNA_seq_clone_end(int);
    
    //Get Mehtods
    String<Dna5> get_RNA_seq_sequence() const;
    string get_RNA_seq_source() const;
    string get_RNA_seq_gene_id() const;
    int get_RNA_seq_gene_strand() const;
    int get_RNA_seq_transcript_id() const;
    long get_RNA_seq_offset() const;
    int get_RNA_seq_clone_end() const;
};

#endif
