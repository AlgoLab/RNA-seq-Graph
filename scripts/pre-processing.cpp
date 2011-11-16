/***********************************/
/* Pre-processing of the reads for */
/* error correction based on a     */
/* consensus procedure with Refseq */
/* sequences                       */
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
//Number of RNASEQ processed simultaneously
#define RNASEQ_BLOCK 10000000

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
    if(argc < 3){
        std::cerr << "Usage: preproc <RNA-seq_file> <Refseq_file> [--genes <Gene1> ... <GeneN>]" << std::endl;
        std::cerr << std::endl;
        return 1;
    }
    vector<string> genes;
    if(argc > 3){
        int j=0;
        for(int i=4; i<argc; ++i){
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
    ifstream f_2(argv[2],::std::ios_base::in|::std::ios_base::binary|ios::ate);
        
    if(f.is_open() && f_2.is_open()){
        file_dimension = f.tellg();
        file_dimension_2 = f_2.tellg();
        f_2.close();
        f.close();
    }
    else{
        std::cerr << "Unable to open files " << std::endl;
        return 1;	
    }
    long read_size = 0;
    int perc = 0;
    std::fstream data_file;
    std::fstream ref;
    std::ofstream out_file;
    std::ofstream not_mapped;
    data_file.open(argv[1], std::ios_base::in | std::ios_base::binary);
    ref.open(argv[2], std::ios_base::in | std::ios_base::binary);
    out_file.open("RNASEQ_NEW.fa");
    not_mapped.open("RNASEQ_NOT_MAPPED");

    //Refseq map: gene_id, sequence(64bp), if_mapped_flag
    //std::map<unsigned long long, pair<pair<int,string>,bool> > ref_map;
    //Refseq map: gene_id, sequence(32bp),next_refseq,if_mapped_flag
    std::map<unsigned long long, pair<pair<pair<int,string>,unsigned long long>,bool> > ref_map;
    //Data Map
    std::map<unsigned long long, unsigned long long > data_map;
    if(data_file.is_open() && ref.is_open() && out_file.is_open() && not_mapped.is_open()){        
        
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
            if ((read_size/(double)file_dimension_2)*100 - 1 >= perc) {
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
                //std::cerr << read_tag.substr(0,pos_f) << " GENE" << std::endl;
                for(unsigned int i = 0;i<=length(read_seq)-READ_LEN;i++){
                    num_fing++;
                    string read;
                    assign(read,infix(read_seq,i,i+READ_LEN));
                    unsigned long long fing = fingerprint(read.substr(0,READ_LEN/2));
                    if(ref_map.find(fing) == ref_map.end()){
                        ref_map[fing] = make_pair(make_pair(make_pair(sp_g,read.substr(0,READ_LEN/2)),
                                                            fingerprint(read.substr(READ_LEN/2))),false);
                    }
                }
            }else{
                for(int i=0; i<argc-4; ++i){
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
                            ref_map[fing] = make_pair(make_pair(make_pair(sp_g,read.substr(0,READ_LEN/2)),
                                                                fingerprint(read.substr(READ_LEN/2))),false);
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
        
        //Mapping
        int found_neigh = 0;

        //RNA-seq Reads
        tStart = clock();
        std::cerr << "Processing RNA-seq file..." << std::endl;
        read_size = 0;
        perc = 0;
        long tot_read = 0;
        while(!data_file.eof()){
            while(num_read < RNASEQ_BLOCK && !data_file.eof()){
                readMeta(data_file, fasta_tag, Fasta());
                read(data_file, read_seq, Fasta());
                read_size+=((length(fasta_tag)*sizeof(char))+(length(read_seq)*sizeof(char)));
                if ((read_size/(double)file_dimension)*100 - 1 >= perc) {
                    perc++;
                    if(perc%10 == 0){
                        std::cerr << "Processed: " << perc << "%" << std::endl;
                    }
                }
                for(unsigned int i = 0;i<=length(read_seq)-READ_LEN;i++){
                    string read;
                    assign(read,infix(read_seq,i,i+READ_LEN));
                    string read_tag;
                    assign(read_tag,fasta_tag);
                    v.push_back(read);
                    tags.push_back(read_tag);
                    num_read = v.size()-1;
                    tot_read ++;
                    unsigned long long left_f = fingerprint(read.substr(0,READ_LEN/2));
                    unsigned long long right_f = fingerprint(read.substr(READ_LEN/2));
                    
                    if(data_map.find(left_f) == data_map.end()){
                        data_map[left_f] = 0;
                    }
                    if(data_map.find(right_f) == data_map.end()){
                        data_map[right_f] = 0;
                    }
                }
            }
            std::cerr << "RNA-seq Extracted Reads: " << tot_read << std::endl;
            int num_neigh = 0;
            //Generazione vicini a distanza 1
            char sub[] = {'A','C','G','T'};
            //std::map<unsigned long long, pair<pair<int, string>,bool> >::iterator it;
            std::map<unsigned long long, pair<pair<pair<int,string>,unsigned long long>,bool> >::iterator it;
            for(it = ref_map.begin(); it != ref_map.end(); ++it){
                string s = it->second.first.first.second;
                //Controllo il read stesso
                if(data_map.find(fingerprint(s))!=data_map.end()){
                    found_neigh++;
                    data_map[fingerprint(s)] = fingerprint(s);
                    it->second.second = true;
                    //std::cerr << s << " " << fingerprint(s) << std::endl;
                }
                //... ed i vicini a distanza 1
                for(unsigned int i=0; i<s.length(); ++i){
                    //std::cout << s << std::endl;
                    for(int j=0; j<4; ++j){
                        if(s.at(i)!=sub[j]){
                            //Vicino a distanza 1
                            string rep = s;
                            rep.replace(i,1,1,sub[j]);
                            num_neigh++;
                            //std::cerr << rep << std::endl;
                            if(data_map.find(fingerprint(rep))!=data_map.end()){
                                found_neigh++;
                                data_map[fingerprint(rep)] = fingerprint(s);
                                it->second.second = true;
                                //std::cerr << rep << std::endl;
                                //std::cerr << data_map[fingerprint(rep)].second << std::endl;
                            }
                        }
                    }
                }
            }
            std::cerr << "Total Refseq Extracted Fingerprints Unique: " << ref_map.size() << std::endl;
            std::cerr << "Total Refseq Fingerprint Neighbor (DH<1) Generated: " << num_neigh << std::endl;
            std::cerr << "Refseq (and Neighbor) - RNAseq Fingerprint Mappings: " << found_neigh << std::endl << std::endl;
            //Re-construction of RNA-seq reads
            for(unsigned int i=0; i<v.size(); ++i){
                //std::cerr << v[i] << "\n";
                //std::cerr << v[i].substr(0,READ_LEN/2) << "\n";
                //std::cerr << v[i].substr(READ_LEN/2) << "\n";
                unsigned long long left_f = fingerprint(v[i].substr(0,READ_LEN/2));
                unsigned long long right_f = fingerprint(v[i].substr(READ_LEN/2));
                
                if(data_map[left_f] == 0 && data_map[right_f] == 0){
                    //std::cout << ">" << tags[i] << "\n";
                    not_mapped << ">" << tags[i] << "\n";
                    //std::cout << v[i] << "\n";
                    not_mapped << v[i] << "\n";
                }else{
                    //std::cout << ">";
                    out_file << ">";
                    string new_read = "";
                    if(data_map[left_f] != 0){
                        //std::cout << "l:" << genes[ref_map[data_map[left_f]].first.first] << "|";
                        out_file << "l:" << genes[ref_map[data_map[left_f]].first.first.first] << "|";
                        new_read += ref_map[data_map[left_f]].first.first.second;
                    }else{
                        //std::cerr << "l: " << v[i].substr(0,READ_LEN/2) << std::endl;
                        new_read += v[i].substr(0,READ_LEN/2);
                    }
                    if(data_map[right_f] != 0){
                        //std::cout << "r:" << genes[ref_map[data_map[right_f]].first.first] << "|";
                        out_file << "r:" << genes[ref_map[data_map[right_f]].first.first.first] << "|";
                        new_read += ref_map[data_map[right_f]].first.first.second;
                    }else{
                        //std::cerr << "r: " << v[i].substr(READ_LEN/2) << std::endl;
                        new_read += v[i].substr(READ_LEN/2);
                    }
                    //std::cout << "mRNA\n";
                    out_file << "mRNA\n";
                    //std::cout << new_read << "\n";
                    out_file << new_read << "\n";
                }
            }
            //Reset tables...
            data_map.clear();
            v.clear();
            tags.clear();
            num_read = 0;
        }
        std::cerr << "Processing RNA-seq file...done!" << std::endl;
        std::cerr << "Processing took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
        std::cerr << " seconds." << std::endl;
        std::cerr << "Total RNA-seq Extracted Fingerprints: " << tot_read*2 << std::endl << std::endl;
        //std::cerr << "RNA-seq Fingerprints (Unique): " << data_map.size() << std::endl;
        data_file.close();
        data_file.clear();
        not_mapped.close();
        //Verify how many Refseq reads are mapped into RNA-seq reads
        long ref_mapped = 0;
        long ref_not_mapped = 0;
        //std::map<unsigned long long, pair<pair<int, string>,bool> >::iterator it;
        std::map<unsigned long long, pair<pair<pair<int,string>,unsigned long long>,bool> >::iterator it;
        std::map<unsigned long long, bool> refseq_not_mapped;
        for(it = ref_map.begin(); it != ref_map.end(); ++it){
            if(it->second.second){
                ref_mapped++;
            }else{ 
                ref_not_mapped++;
                refseq_not_mapped[it->first] = false;
            }
        }
        std::cerr << "Total Refseq Extracted Fingerprints Unique: " << ref_map.size() << std::endl;
        std::cerr << "Total Refseq Fingerprints Mapped into RNA-seq Fingerprints: " << ref_mapped << std::endl;
        std::cerr << "Total Refseq Fingerprints NOT Mapped into RNA-seq Fingerprints: " << ref_not_mapped << std::endl;
        std::cerr << std::endl;
        out_file.flush();

        //REVERSE AND COMLPEMENT
        num_fing = 0;
        read_size = 0;
        perc = 0;
        v.clear();
        tags.clear();
        num_read = 0;
        ref_map.clear();
        data_map.clear();
        file_dimension = 0;
        //Files
        ifstream f_tmp("RNASEQ_NOT_MAPPED",::std::ios_base::in|::std::ios_base::binary|ios::ate);
        if(f_tmp.is_open()){
            file_dimension = f_tmp.tellg();
            f_tmp.close();
        }else{
            std::cerr << "Unable to open files " << std::endl;
            return 1;	
        }
        data_file.open("RNASEQ_NOT_MAPPED", std::ios_base::in | std::ios_base::binary);
        ref.open(argv[2], std::ios_base::in | std::ios_base::binary);
        //Check if there are unmapped reads
        if(!data_file.is_open() || !ref.is_open() || file_dimension == 0){
            std::cerr << "No unmapped RNA-seq" << std::endl;
            return 1;
        }
        std::cerr << "Processing Reversed and Complemented Sequences ..." << std::endl;
        tStart = clock();
        while(!ref.eof()){
            readMeta(ref, fasta_tag, Fasta());
            read(ref, read_seq, Fasta());
            //Reversed and complemented Sequence
            String<Dna5> rev_comp = Dna5StringReverseComplement(read_seq);
            read_seq = rev_comp;
            string read_tag;
            assign(read_tag,fasta_tag);
            read_size+=((length(fasta_tag)*sizeof(char))+(length(read_seq)*sizeof(char)));
            if ((read_size/(double)file_dimension_2)*100 - 1 >= perc) {
                perc++;
                if(perc%10 == 0){
                    std::cerr << "Processed: " << perc << "%" << std::endl;
                }
            }
            int sp_g = -1;
            if(!exp_genes){
                //Find gene name
                size_t pos_g = read_tag.find_last_of(" ");
                genes.push_back(read_tag.substr(0,pos_g));
                sp_g = genes.size() - 1;
                //std::cout << read_tag.substr(0,pos_g) << std::endl;
                for(unsigned int i = 0;i<=length(read_seq)-READ_LEN;i++){
                    num_fing++;
                    string read;
                    assign(read,infix(read_seq,i,i+READ_LEN));
                    unsigned long long fing = fingerprint(read.substr(0,READ_LEN/2));
                    if(ref_map.find(fing) == ref_map.end()){
                        ref_map[fing] = make_pair(make_pair(make_pair(sp_g,read.substr(0,READ_LEN/2)),
                                                            fingerprint(read.substr(READ_LEN/2))),false);
                    }
                }
            }else{
                for(int i=0; i<argc-4; ++i){
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
                            ref_map[fing] = make_pair(make_pair(make_pair(sp_g,read.substr(0,READ_LEN/2)),
                                                                fingerprint(read.substr(READ_LEN/2))),false);
                        }
                    }
                }
            }
        }
        ref.close();
        std::cerr << "Processing Reversed and Complemented Sequences ...done!" << std::endl;
        std::cerr << "Processing took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
        std::cerr << " seconds." << std::endl;
        std::cerr << "Total Refseq Extracted Fingerprints: " << num_fing << std::endl;
        std::cerr << "Total Refseq Extracted Fingerprints Unique: " << ref_map.size() << std::endl << std::endl;
        
        tStart = clock();
        std::cerr << "Processing RNA-seq file..." << std::endl;
        read_size = 0;
        perc = 0;
        tot_read = 0;
        while(!data_file.eof()){
            while(num_read < RNASEQ_BLOCK && !data_file.eof()){
                readMeta(data_file, fasta_tag, Fasta());
                read(data_file, read_seq, Fasta());
                read_size+=((length(fasta_tag)*sizeof(char))+(length(read_seq)*sizeof(char)));
                if ((read_size/(double)file_dimension)*100 - 1 >= perc) {
                    perc++;
                    if(perc%10 == 0){
                        std::cerr << "Processed: " << perc << "%" << std::endl;
                    }
                }
                for(unsigned int i = 0;i<=length(read_seq)-READ_LEN;i++){
                    string read;
                    assign(read,infix(read_seq,i,i+READ_LEN));
                    string read_tag;
                    assign(read_tag,fasta_tag);
                    v.push_back(read);
                    tags.push_back(read_tag);
                    num_read = v.size()-1;
                    tot_read ++;
                    unsigned long long left_f = fingerprint(read.substr(0,READ_LEN/2));
                    unsigned long long right_f = fingerprint(read.substr(READ_LEN/2));
                    
                    if(data_map.find(left_f) == data_map.end()){
                        data_map[left_f] = 0;
                    }
                    if(data_map.find(right_f) == data_map.end()){
                        data_map[right_f] = 0;
                    }
                }
            }
            std::cerr << "RNA-seq Extracted Reads: " << tot_read << std::endl;
            int num_neigh = 0;
            //Generazione vicini a distanza 1
            char sub[] = {'A','C','G','T'};
            //std::map<unsigned long long, pair<pair<int, string>,bool> >::iterator it;
            std::map<unsigned long long, pair<pair<pair<int,string>,unsigned long long>,bool> >::iterator it;
            for(it = ref_map.begin(); it != ref_map.end(); ++it){
                string s = it->second.first.first.second;
                //Controllo il read stesso
                if(data_map.find(fingerprint(s))!=data_map.end()){
                    found_neigh++;
                    data_map[fingerprint(s)] = fingerprint(s);
                    it->second.second = true;
                    //std::cerr << s << " " << fingerprint(s) << std::endl;
                }
                //... ed i vicini a distanza 1
                for(unsigned int i=0; i<s.length(); ++i){
                    //std::cout << s << std::endl;
                    for(int j=0; j<4; ++j){
                        if(s.at(i)!=sub[j]){
                            //Vicino a distanza 1
                            string rep = s;
                            rep.replace(i,1,1,sub[j]);
                            num_neigh++;
                            //std::cerr << rep << std::endl;
                            if(data_map.find(fingerprint(rep))!=data_map.end()){
                                found_neigh++;
                                data_map[fingerprint(rep)] = fingerprint(s);
                                it->second.second = true;
                                //std::cerr << rep << std::endl;
                                //std::cerr << data_map[fingerprint(rep)].second << std::endl;
                            }
                        }
                    }
                }
            }
            std::cerr << "Total Refseq Extracted Fingerprints Unique: " << ref_map.size() << std::endl;
            std::cerr << "Total Refseq Fingerprint Neighbor (DH<1) Generated: " << num_neigh << std::endl;
            std::cerr << "Refseq (and Neighbor) - RNAseq Fingerprint Mappings: " << found_neigh << std::endl << std::endl;
            //Re-construction of RNA-seq reads
            for(unsigned int i=0; i<v.size(); ++i){
                //std::cerr << v[i] << "\n";
                //std::cerr << v[i].substr(0,READ_LEN/2) << "\n";
                //std::cerr << v[i].substr(READ_LEN/2) << "\n";
                unsigned long long left_f = fingerprint(v[i].substr(0,READ_LEN/2));
                unsigned long long right_f = fingerprint(v[i].substr(READ_LEN/2));
                
                if(data_map[left_f] != 0 || data_map[right_f] != 0){
                    //std::cout << ">";
                    out_file << ">";
                    string new_read = "";
                    if(data_map[left_f] != 0){
                        //std::cout << "l:" << genes[ref_map[data_map[left_f]].first.first] << "|";
                        out_file << "l:" << genes[ref_map[data_map[left_f]].first.first.first] << "|";
                        new_read += ref_map[data_map[left_f]].first.first.second;
                    }else{
                        //std::cerr << "l: " << v[i].substr(0,READ_LEN/2) << std::endl;
                        new_read += v[i].substr(0,READ_LEN/2);
                    }
                    if(data_map[right_f] != 0){
                        //std::cout << "r:" << genes[ref_map[data_map[right_f]].first.first] << "|";
                        out_file << "r:" << genes[ref_map[data_map[right_f]].first.first.first] << "|";
                        new_read += ref_map[data_map[right_f]].first.first.second;
                    }else{
                        //std::cerr << "r: " << v[i].substr(READ_LEN/2) << std::endl;
                        new_read += v[i].substr(READ_LEN/2);
                    }
                    //std::cout << "mRNA\n";
                    out_file << "mRNA\n";
                    //std::cout << new_read << "\n";
                    out_file << Dna5StringReverseComplement(new_read) << "\n";
                }
            }
            //Reset tables...
            data_map.clear();
            v.clear();
            tags.clear();
            num_read = 0;
        }
        std::cerr << "Processing RNA-seq file...done!" << std::endl;
        std::cerr << "Processing took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
        std::cerr << " seconds." << std::endl;
        std::cerr << "Total RNA-seq Extracted Fingerprints: " << tot_read*2 << std::endl << std::endl;
        //std::cerr << "RNA-seq Fingerprints (Unique): " << data_map.size() << std::endl;
        data_file.close();
        out_file.close();
    }else{
        std::cerr << "Unable to open files" << std::endl;
        std::cerr << std::endl;
        return 1;
    }
    return 0;
}
