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

//Compute the edit distance between two strings
unsigned int levenshtein_distance(const string &s1, const string & s2) {
        const size_t len1 = s1.size(), len2 = s2.size();
        vector<unsigned int> col(len2+1), prevCol(len2+1);
 
        for (unsigned int i = 0; i < prevCol.size(); i++)
                prevCol[i] = i;
        for (unsigned int i = 0; i < len1; i++) {
                col[0] = i+1;
                for (unsigned int j = 0; j < len2; j++)
                        col[j+1] = min( min( 1 + col[j], 1 + prevCol[1 + j]), prevCol[j] + (s1[i]==s2[j] ? 0 : 1) );
                col.swap(prevCol);
        }
        return prevCol[len2];
}

int main(int argc, char* argv[]){
    if(argc < 3){
        std::cerr << "Usage: preproc <RNA-seq_file> <Refseq_file>" << std::endl;
        std::cerr << std::endl;
        return 1;
    }

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
    out_file.open("RNASEQ_annotated.fa");
    not_mapped.open("RNASEQ_not_annotated.fa");

    //Refseq seuqences
    std::vector<String<Dna5> > ref_sequences;
    //Refseq gene noames
    std::vector<string> genes;
    //Refseq map: id refseq sequence, position
    std::map<unsigned long long, vector<pair<int,int > > > ref_map;
 
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
            //Tags
            readMeta(ref, fasta_tag, Fasta());
            string read_tag;
            assign(read_tag,fasta_tag);
            //Sequence
            read(ref, read_seq, Fasta());
            read_size+=((length(fasta_tag)*sizeof(char))+(length(read_seq)*sizeof(char)));
            if ((read_size/(double)file_dimension_2)*100 - 1 >= perc) {
                perc++;
                if(perc%10 == 0){
                    std::cerr << "Processed: " << perc << "%" << std::endl;
                }
            }
            if(length(read_seq) >= READ_LEN/2){
                //Find gene name
                size_t pos_g = read_tag.find(" ");
                genes.push_back(read_tag.substr(0,pos_g));
                ref_sequences.push_back(read_seq);
                int sp_g = ref_sequences.size() -1;
                //std::cout << "Refseq: " << sp_g  << " " << read_tag.substr(0,pos_g) << " " << length(read_seq)-READ_LEN/2 << " " << length(read_seq) << " " << READ_LEN/2 << std::endl;
                //std::cerr << read_tag.substr(0,pos_f) << " GENE" << std::endl;
                for(unsigned int i = 0;i<=length(read_seq)-READ_LEN/2;i++){
                    num_fing++;
                    string read;
                    assign(read,infix(read_seq,i,i+READ_LEN/2));
                    unsigned long long fing = fingerprint(read);
                    //if(ref_map.find(fing) == ref_map.end()){
                    ref_map[fing].push_back(make_pair(sp_g,i)); 
                    //}   
                    //ref_map[fing] = make_pair(make_pair(make_pair(sp_g,read.substr(0,READ_LEN/2)),
                    //                                    fingerprint(read.substr(READ_LEN/2))),false);
                    //}else{
                    
                    //}
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
        //String<Dna5> rev_comp = Dna5StringReverseComplement(read_seq);
        //RNA-seq Reads
        tStart = clock();
        std::cerr << "Processing RNA-seq file..." << std::endl;
        read_size = 0;
        perc = 0;
        long tot_read = 0;
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
            //for(unsigned int i = 0;i<=length(read_seq)-READ_LEN;i++){ //All 64bp Reads
            for(unsigned int i = 0;i<1;i++){ //Only 1Read: 1x75bp = 1x64bp
                string read;
                assign(read,infix(read_seq,i,i+READ_LEN));
                string read_tag;
                assign(read_tag,fasta_tag);
                tot_read++;
                unsigned long long left_f = fingerprint(read.substr(0,READ_LEN/2));
                unsigned long long right_f = fingerprint(read.substr(READ_LEN/2));
                bool left_mapped = false;
                bool right_mapped = false;
                bool mapped = false;
                if(ref_map.find(left_f) != ref_map.end() && ref_map.find(right_f) != ref_map.end()){
                    std::cout << "Caso 1" << std::endl;
                    out_file << ">";
                    out_file << "l:" << genes[ref_map[left_f][0].first] << "|";
                    out_file << "r:" << genes[ref_map[right_f][0].first] << std::endl;
                    out_file << read_seq << std::endl;
                    mapped = true;
                }else{
                    if(ref_map.find(left_f) != ref_map.end()){
                        std::cout << "Caso 2" << std::endl;
                        for(int j=0; j<ref_map[left_f].size();++j){
                            string ref_right;
                            int pos = ref_map[left_f].at(j).second + (READ_LEN/2);
                            assign(ref_right,infix(ref_sequences[ref_map[left_f].at(j).first],
                                                   pos,pos+READ_LEN/2));
                                                      
                            if(levenshtein_distance(read.substr(READ_LEN/2),ref_right)<=1){
                                out_file << ">";
                                out_file << "l:" << genes[ref_map[left_f].at(j).first] << "|";
                                out_file << "r:" << genes[ref_map[left_f].at(j).first] << std::endl;
                                out_file << read.substr(0,READ_LEN/2) + ref_right << std::endl;
                                mapped = true;
                            }
                        }
                    }
                    if(ref_map.find(right_f) != ref_map.end()){
                        std::cout << "Caso 3" << std::endl;
                        for(int j=0; j<ref_map[right_f].size();++j){
                            string ref_left;
                            assign(ref_left,infix(ref_sequences[ref_map[right_f].at(j).first],
                                                   ref_map[right_f].at(j).second-READ_LEN/2,READ_LEN/2));
                           
                            if(levenshtein_distance(read.substr(0,READ_LEN/2),ref_left)<=1){
                                out_file << ">";
                                out_file << "l:" << genes[ref_map[right_f].at(j).first] << "|";
                                out_file << "r:" << genes[ref_map[right_f].at(j).first] << std::endl;
                                out_file << ref_left + read.substr(READ_LEN/2) << std::endl;
                                mapped = true;
                            }
                        }
                    }
                }
                if(!mapped){
                    not_mapped << ">" << read_tag << std::endl;
                    not_mapped << read << std::endl;
                }
            }
        }
        std::cerr << "RNA-seq Extracted Reads: " << tot_read << std::endl;
            /*
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
            */
        ref.close();
        data_file.close();
        out_file.close();
    }else{
        std::cerr << "Unable to open files" << std::endl;
        std::cerr << std::endl;
        return 1;
    }
    return 0;
}
