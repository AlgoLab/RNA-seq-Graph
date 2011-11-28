/**
* RNA-seq-Graph
* Method for reconstructing the Isoform Graph of a gene from RNA-seq data, without the genome information
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
* A full copy of the GNU Affero General Public License is reported into
* the file COPYING. However more informations ca be found at
* <http://www.gnu.org/licenses/>.
**/
/***********************/
/* Missing Refseq      */
/***********************/
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <set>

#include <seqan/sequence.h>
#include <seqan/file.h>

#define READ_LEN 64

using namespace seqan;
using namespace std;

int main(int argc, char* argv[]){
    if(argc < 2){
        std::cerr << "Usage: extract-reads <Refseq_file> [--genes <Gene1> ... <GeneN>]" << std::endl;
        std::cerr << std::endl;
        return 1;
    }
    vector<string> genes;
    if(argc > 2){
        int j=0;
        for(int i=3; i<argc; ++i){
            genes.push_back(string(argv[i]));
            //std::cerr << genes[j] << std::endl;
            ++j;
        }
    }
    bool exp_genes = false;
    if(genes.size() > 0) exp_genes = true;

    long file_dimension = 0;

    //Files
    ifstream f(argv[1],::std::ios_base::in|::std::ios_base::binary|ios::ate);

    if(f.is_open()){
        file_dimension = f.tellg();
        f.close();
    }
    else{
        std::cerr << "Unable to open file " << argv[1] << std::endl;
        return 1;	
    }
    long read_size = 0;
    int perc = 0;
    std::fstream ref;
    ref.open(argv[1], std::ios_base::in | std::ios_base::binary);
    std::ofstream out_file;
    out_file.open("Refseq_reads.fa");

    if(ref.is_open() && out_file.is_open()){        
        String<char> fasta_tag;
        String<Dna5> read_seq;

        //Ref sequences
        long num_reads = 0;
        std::cerr << "Processing Reference Sequences file..." << std::endl;
        clock_t tStart = clock();
        
        while(!ref.eof()){
            readMeta(ref, fasta_tag, Fasta());
            read(ref, read_seq, Fasta());
            string read_tag;
            assign(read_tag,fasta_tag);
            read_size+=((length(fasta_tag)*sizeof(char))+(length(read_seq)*sizeof(char)));
            if ((read_size/(double)file_dimension)*100 - 1 >= perc) {
                perc++;
                if(perc%10 == 0){
                    std::cerr << "Processed: " << perc << "%" << std::endl;
                }
            }
            int sp_g = -1;
            if(!exp_genes){
                //Find gene name
                size_t pos_g = read_tag.find(" ");

                for(unsigned int i = 0;i<=length(read_seq)-READ_LEN;i++){
                    num_reads++;
                    out_file << ">" << read_tag.substr(0,pos_g) << "_" << i+1 << std::endl;
                    out_file << infix(read_seq,i,i+READ_LEN) << std::endl;
                }
            }else{ 
                for(int i=0; i<argc-3; ++i){
                    if(read_tag.find(genes[i]) != string::npos){
                        sp_g = i;
                    }
                }
                if(sp_g != -1){
                    for(unsigned int i = 0;i<=length(read_seq)-READ_LEN;i++){
                        num_reads++;
                        out_file << ">" << genes[sp_g] << "_" << i+1 << std::endl;
                        out_file << infix(read_seq,i,i+READ_LEN) << std::endl;
                    }
                }
            }
        }
        std::cerr << "Processing Reference Sequences file...done!" << std::endl;
        std::cerr << "Processing took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
        std::cerr << " seconds." << std::endl;
        std::cerr << "Total Refseq Extracted Reads: " << num_reads << std::endl;

        ref.close();
        out_file.close();
    }
}
