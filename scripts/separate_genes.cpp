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
/***********************************/
/* Pre-processing of the reads for */
/* error correction based on a     */
/* consensus procedure with a      */
/* frequence threshold for SNPs    */
/***********************************/

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
    if(argc < 3){
        std::cerr << "Usage: separate_genes <RNA-seq_file> --gene <Gene>" << std::endl;
        std::cerr << std::endl;
        return 1;
    }
  
    long file_dimension = 0;
    ifstream f(argv[1],::std::ios_base::in|::std::ios_base::binary|ios::ate);

    if(f.is_open() ){
        file_dimension = f.tellg();
        f.close();
    }
    else{
        std::cerr << "Unable to open files " << std::endl;
        return 1;	
    }
    long read_size = 0;
    int perc = 0;
    std::fstream data_file;
     data_file.open(argv[1], std::ios_base::in | std::ios_base::binary);
     if(data_file.is_open()){
         string gene(argv[3]);
         std::cerr << "Filtering " << gene << " gene" << std::endl;
         String<char> fasta_tag;
         String<Dna5> read_seq;
         
         //RNA-seq Reads
         clock_t tStart = clock();
         std::cerr << "Processing RNA-seq file..." << std::endl;
         while(!data_file.eof()){
             readMeta(data_file, fasta_tag, Fasta());
             read(data_file, read_seq, Fasta());
             read_size+=((length(fasta_tag)*sizeof(char))+(length(read_seq)*sizeof(char)));
             if ((read_size/(double)file_dimension)*100 - 1 >= perc) {
                 perc++;
                 if(perc%10 == 0){
                     std::cerr << "Processed: " << perc << "%" << std::endl;
                 }
             }
             
             string read_tag;
             assign(read_tag,fasta_tag);
             if(read_tag.find(gene) != string::npos){
                 std::cout << ">" << read_tag << std::endl;
                 std::cout << read_seq << std::endl;
             }
         }
         
         std::cerr << "Processing RNA-seq file...done!" << std::endl;
         std::cerr << "Processing took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
         std::cerr << " seconds." << std::endl << std::endl;
         data_file.close();
         
     }else{
         std::cerr << "Unable to open RNA-seq file" << std::endl;
         std::cerr << std::endl;
         return 1;
     }
     return 0;
}
