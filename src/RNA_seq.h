/**
* RNA-seq-Graph
* Method for reconstructing the Isoform Graph of a gene from RNA-seq data
*
* Copyright (C) 2011 Stefano Beretta <ste.beretta(-at-)gmail.com>
*
* Distributed under the terms of the GNU Affero General Public License (AGPL)
*
*
* This file is part of RNA-seq-Graph.
*
* RNA-seq-Graph is free software: you can redistribute it and/or modify
* it under the terms of the GNU Affero General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* RNA-seq-Graph is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU Affero General Public License for more details.
*
* You should have received a copy of the GNU Affero General Public License
* along with RNA-seq-Graph.  If not, see <http://www.gnu.org/licenses/>.
**/
/********************************/
/* Class the models the RNA-seq */
/* sequence                     */
/********************************/
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
    //string source;
    //string gene_id;
    //int gene_strand;
    //int transcript_id;
    //Offset starting from the first digit of the transcript
    //long offset;
    //int clone_end;
    
    public:
    //Constructors
    RNA_seq(const String<Dna5>&);
    //RNA_seq(String<Dna5>, string, string);
    //RNA_seq(String<Dna5>, string, string, int);
    //RNA_seq(String<Dna5>, string, string, int, int, long, int);
    //Copy Constructor    
    RNA_seq(const RNA_seq&);
    //Overloading =
    RNA_seq& operator=(const RNA_seq&);
    //Set Methods
    void set_RNA_seq_sequence(const String<Dna5>& seq) {
        sequence = seq;
    }
    //void set_RNA_seq_source(string);
    //void set_RNA_seq_gene_id(string);
    //void set_RNA_seq_gene_strand(int);
    //void set_RNA_seq_transcript_id(int);
    //void set_RNA_seq_offset(long);
    //void set_RNA_seq_clone_end(int);
    
    //Get Mehtods
    const String<Dna5>& get_RNA_seq_sequence() const {
        return sequence;
    }
    //string get_RNA_seq_source() const;
    //string get_RNA_seq_gene_id() const;
    //int get_RNA_seq_gene_strand() const;
    //int get_RNA_seq_transcript_id() const;
    //long get_RNA_seq_offset() const;
    //int get_RNA_seq_clone_end() const;
};

#endif
