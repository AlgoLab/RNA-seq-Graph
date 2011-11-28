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
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;
using namespace std;

int main(int argc, char* argv[]){
    clock_t tStart = clock();

    MultiSeqFile multiSeqFile;
    if (argc < 4 || !open(multiSeqFile.concat, argv[3], OPEN_RDONLY)){
        std::cout << "Usage: Filter_RNA-seq <quality-threshold> <read-length> <RNA-seq-file>" << std::endl;
        std::cout << std::endl;
        return 1;
    }
    //Output File  
    string name(argv[3]);
    size_t k = name.rfind(".");
    if(k != string::npos && k>0){
        name = name.substr(0,k);
    }
    
    name += ".fa";
    std::cout << name << std::endl;
    std::ofstream out_file;
    out_file.open(name.c_str());

    //Threshold
    int threshold;
    std::istringstream ss(argv[1]);
    ss >> threshold;
    //Read Length
    int read_length;
    std::istringstream ss_length(argv[2]);
    ss_length >> read_length;

    //RNASeq Processing
    AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);
    unsigned seqCount = length(multiSeqFile); 

    String<Dna5Q> seq;
    CharString qual;
    CharString id;
    for (unsigned i = 0; i < seqCount; ++i){
        assignSeq(seq, multiSeqFile[i], format);    // read sequence
        assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
        assignSeqId(id, multiSeqFile[i], format);   // read sequence id

        //std::cout << id << std::endl;
        //std::cout << seq << std::endl;
        // convert ascii to values from 0..62
        int pos = 0;
        string filt_seq = "";
        for (unsigned j = 0; j < length(qual) && j < length(seq); ++j){
            assignQualityValue(seq[j], (int)(ordValue(qual[j]) - 33));
            //std::cout << getQualityValue(seq[j]) << " ";
            if(getQualityValue(seq[j]) >= threshold){
                filt_seq += seq[j];
                pos++;
            }else{
                if(pos >= read_length){
                    //std::cout << filt_seq << std::endl;
                    write(out_file, filt_seq, id, Fasta());
                }
                pos=0;
                filt_seq = "";
            }
        }
        if(pos >= read_length){
            //std::cout << filt_seq << std::endl;
            write(out_file, filt_seq, id, Fasta());
        }
        //std::cout << std::endl;
    }
    
    std::cerr << "Loading " << seqCount << " sequences took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    std::cerr << " seconds." << std::endl << std::endl;
    out_file.close();
    return 0;
}
