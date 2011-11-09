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

/************************************************/
/* Convert a DNA sequnece on alphabet {A,C,G,T} */
/* into a binary sequence                       */
/************************************************/
unsigned long long fingerprint(const string& seq){
    unsigned long long number = 0;
    for(unsigned int i=0; i<seq.length(); i++){
        number = number<<2;
        if(seq.at(i) == 'N' || seq.at(i) == 'n'){
            number |= 0;
        }
        if(seq.at(i) == 'A' || seq.at(i) == 'a'){
            number |= 0;
        }
        if(seq.at(i) == 'C' || seq.at(i) == 'c'){
            number |= 1;
        }
        if(seq.at(i) == 'G' || seq.at(i) == 'g'){
            number |= 2;
        }
        if(seq.at(i) == 'T' || seq.at(i) == 't'){
            number |= 3;
        }
    }//End_For
    return number;
}//End_Method

int main(int argc, char* argv[]){
    if(argc < 2){
        std::cerr << "Usage: missing-ref <Refseq_file> [--genes <Gene1> ... <GeneN>]" << std::endl;
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
    long file_dimension_2 = 0;
    //Files
    ifstream f(argv[1],::std::ios_base::in|::std::ios_base::binary|ios::ate);
    ifstream f_2("RNASEQ_NEW.fa",::std::ios_base::in|::std::ios_base::binary|ios::ate);

    if(f.is_open() && f_2.is_open()){
        file_dimension = f.tellg();
        file_dimension_2 = f_2.tellg();
        f.close();
        f_2.close();
    }
    else{
        std::cerr << "Unable to open file " << argv[1] << " or RNASEQ_NEW.fa" << std::endl;
        return 1;	
    }
    long read_size = 0;
    int perc = 0;
    std::fstream data_file;
    data_file.open("RNASEQ_NEW.fa", std::ios::in | std::ios_base::binary);
    std::fstream ref;
    ref.open(argv[1], std::ios_base::in | std::ios_base::binary);
    std::ofstream out_file;
    out_file.open("MISSING_REFSEQ.fa");
    //Keep also 2 strings: the one preceding and the one following the extrected one
    std::map<unsigned long long, pair<string, bool> > ref_map;

    if(data_file.is_open() && ref.is_open()){        
        
        String<char> fasta_tag;
        String<Dna5> read_seq;
        std::vector<string> v;
        std::vector<String<char> > tags;
        //Posizione del read nel vettore v
        long num_read = 0;

        //Ref sequences
        long num_fing = 0;
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
                genes.push_back(read_tag.substr(0,pos_g));
                sp_g = genes.size() -1;
                //std::cout << read_tag.substr(pos_s+1,pos_f-pos_s-1) << std::endl;
                for(unsigned int i = 0;i<=length(read_seq)-READ_LEN;i++){
                    num_fing++;
                    string read;
                    assign(read,infix(read_seq,i,i+READ_LEN));
                    unsigned long long fing = fingerprint(read.substr(0,READ_LEN/2));
                    if(ref_map.find(fing) == ref_map.end()){
                        ref_map[fing] = make_pair(read,false);
                    }else{
                        std::cerr << "DUPLICATO" << std::endl;
                    }
                }
            }else{
                for(int i=0; i<argc-3; ++i){
                    if(read_tag.find(genes[i]) != string::npos){
                        sp_g = i;
                    }
                }
                if(sp_g != -1){
                    for(unsigned int i = 0;i<=length(read_seq)-READ_LEN;i++){
                        num_fing++;
                        string read;
                        assign(read,infix(read_seq,i,i+READ_LEN));
                        unsigned long long fing = fingerprint(read.substr(0,READ_LEN/2));
                        if(ref_map.find(fing) == ref_map.end()){
                            ref_map[fing] = make_pair(read,false);
                        }else{
                            std::cerr << "DUPLICATO" << std::endl;
                        }
                    }
                }
            }
        }
        ref.close();
        std::cerr << "Processing Reference Sequences file...done!" << std::endl;
        std::cerr << "Processing took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
        std::cerr << " seconds." << std::endl;
        std::cerr << "Total Refseq Extracted Fingerprints: " << num_fing << std::endl;
        std::cerr << "Total Refseq Extracted Fingerprints Unique: " << ref_map.size() << std::endl << std::endl;

        //RNA-seq
        tStart = clock();
        std::cerr << "Processing RNA-seq file..." << std::endl;
        read_size = 0;
        perc = 0;
        while(!data_file.eof()){
            while(num_read < 5000000 && !data_file.eof()){
                readMeta(data_file, fasta_tag, Fasta());
                read(data_file, read_seq, Fasta());
                string read_tag;
                assign(read_tag,fasta_tag);
                read_size+=((length(fasta_tag)*sizeof(char))+(length(read_seq)*sizeof(char)));
                if ((read_size/(double)file_dimension_2)*100 - 1 >= perc) {
                    perc++;
                    if(perc%10 == 0){
                        std::cerr << "Processed: " << perc << "%" << std::endl;
                    }
                }
                if(read_tag.find("l:") != string::npos){
                    string read;
                    assign(read,prefix(read_seq,READ_LEN/2));
                    if(ref_map.find(fingerprint(read)) != ref_map.end()){
                        ref_map[fingerprint(read)].second = true;
                    }
                }
                if(read_tag.find("r:") != string::npos){
                    string read;
                    assign(read,suffix(read_seq,READ_LEN/2));
                    if(ref_map.find(fingerprint(read)) != ref_map.end()){
                        ref_map[fingerprint(read)].second = true;
                    }
                }
            }
            std::cerr << "Processing RNA-seq file...done!" << std::endl;
            std::cerr << "Processing took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
            std::cerr << " seconds." << std::endl;
        }
        std::map<unsigned long long, pair<string, bool> >::iterator it;
        for(it = ref_map.begin(); it != ref_map.end(); ++it){
            if(!it->second.second){
                out_file << ">REFSEQ" << std::endl;
                out_file << it->second.first << std::endl;
            }
        }
        data_file.close();
        ref.close();
        out_file.close();
    }
}