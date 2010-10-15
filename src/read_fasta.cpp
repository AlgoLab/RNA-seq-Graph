#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <sstream>

#include <seqan/sequence.h>
#include <seqan/file.h>

#include "read_fasta.h"

using namespace seqan;
using namespace std;

/************************************************/
/* Convert a DNA sequnece on alphabet {A,C,G,T} */
/* into a binary sequence                       */
/************************************************/
string binary_conversion(string seq){
    string bin = "";
    for(unsigned int i=0; i<seq.length(); i++){
        if(seq.at(i) == 'A' || seq.at(i) == 'a'){
            bin += "00";
        }
        if(seq.at(i) == 'C' || seq.at(i) == 'c'){
            bin += "01";
        }
        if(seq.at(i) == 'G' || seq.at(i) == 'g'){
            bin += "10";
        }
        if(seq.at(i) == 'T' || seq.at(i) == 't'){
            bin += "11";
        }
    }//End_For
    return bin;
}//End_Method

/**********************************************/
/* Encode a binary sequence (max 64 bit) into */
/* a decimal (long long) number               */
/**********************************************/
unsigned long long fingerprint(string s){
    string bin = binary_conversion(s);
    unsigned long long number = 0;
    for(unsigned int i=0; i<bin.length(); i++){
        number = number<<1;
        if(bin.at(i)=='1'){
            number |= 1;
        }//End_If
    }//End_For
    return number;
}//End_Method

////Parse Fasta Information
table_entry* parse_fasta(String<Dna5> seq, string meta){
    table_entry* el = NULL;
    string left_seq;
    assign(left_seq,prefix(seq,length(seq)/2));
    string right_seq;
    assign(right_seq,suffix(seq,length(seq)/2));
    int source_pos = meta.find("/source=");
    int gene_id_pos = meta.find("/gene_id=");
    int gene_strand_pos = meta.find("/gene_strand=");

    if(source_pos != -1 && gene_id_pos != -1 && gene_strand_pos != -1){
        source_pos += 8;
        gene_id_pos += 9;
        gene_strand_pos += 13;
        int next_stop = meta.find_first_of(" ",source_pos);
        string source = meta.substr(source_pos,next_stop - source_pos);

        next_stop = meta.find_first_of(" ",gene_id_pos);
        string gene_id = meta.substr(gene_id_pos,next_stop - gene_id_pos);

        next_stop = meta.find_first_of(" ",gene_strand_pos);
        int gene_strand;
        ::std::istringstream(meta.substr(gene_strand_pos,next_stop - gene_strand_pos)) >> gene_strand;

        int gb_pos = meta.find("/gb=", 0);
        int clone_end_pos = meta.find("/clone_end=", 0);

        if(gb_pos != -1 && clone_end_pos != -1){
            gb_pos += 4;
            clone_end_pos += 11;
            next_stop = meta.find_first_of("_",gb_pos);
            int transcript_id;
            ::std::istringstream(meta.substr(gb_pos,next_stop - gb_pos)) >> transcript_id;

            gb_pos = next_stop+1;
            next_stop = meta.find_first_of(" ",gb_pos);
            long offset;
            ::std::istringstream(meta.substr(gb_pos,next_stop - gb_pos)) >> offset;


            next_stop = meta.find_first_of(" ",clone_end_pos);
            int clone_end;
            ::std::istringstream(meta.substr(clone_end_pos,next_stop - 1 - clone_end_pos)) >> clone_end;

            el = new table_entry(seq, source, gene_id, gene_strand, transcript_id, offset, clone_end,
                                 fingerprint(left_seq), fingerprint(right_seq));
        }else{
            el = new table_entry(seq, source, gene_id, gene_strand, fingerprint(left_seq), fingerprint(right_seq));
        }//End_If
    }else{
        ::std::cout << "Invalid fasta entry" << ::std::endl;
    }//End_If
    return el;
}//End_Method

///Add entries in the "Hash Table"
void add_entry(Map &m, table_entry* entry, char l_r){
    if(l_r == 'l'){
        Map::iterator it;
        it = m.find(entry->get_left_fingerprint());
        if(it == m.end()){
            element_table el;
            el.unspliced = 1;
            el.p = entry;
            m[entry->get_left_fingerprint()] = el;
        }else{
            table_entry* temp = m[entry->get_left_fingerprint()].p;
            table_entry* prev;
            bool found = 0;
            do{
                if(temp->get_right_fingerprint() == entry->get_right_fingerprint()){
                    found = 1;
                }
                prev = temp;
                temp = temp->get_next();
            }while(temp != NULL && !found);
            if(found){
                prev->increase_freq();
					 delete entry;
            }else{
		m[entry->get_left_fingerprint()].unspliced = 0;
                prev->set_next(entry);
                entry->set_prev(prev);
            }//End_If
        }//End_If
    }else{
        Map::iterator it;
        it = m.find(entry->get_right_fingerprint());
        if(it == m.end()){
            element_table el;
            el.unspliced = 1;
            el.p = entry;
            m[entry->get_right_fingerprint()] = el;
        }else{
            table_entry* temp = m[entry->get_right_fingerprint()].p;
            table_entry* prev;
            bool found = 0;
            do{
                if(temp->get_left_fingerprint() == entry->get_left_fingerprint()){
                    found = 1;
                    //break;
                }
                prev = temp;
                temp = temp->get_next();
            }while(temp != NULL && !found);
            if(found){
                prev->increase_freq();
					 delete entry;
            }else{
		m[entry->get_right_fingerprint()].unspliced = 0;
                prev->set_next(entry);
                entry->set_prev(prev);
            }//End_If
        }//End_I
    }//End_If
}//End_Method

///Read a fasta file
tables read_fasta(char* file_name){
    //Struct with hash tables
    tables t;

    string name = file_name;
    string extension = name.substr(name.find_last_of(".") + 1);
    if(extension == "fa" || extension == "fas" || extension == "fasta"){
        ::std::fstream fstrm;
        fstrm.open(file_name, ::std::ios_base::in | ::std::ios_base::binary);
        if(fstrm.is_open()){
            String<char> fasta_tag;
            String<Dna5> fasta_seq;
            int c = 0;
            while(!fstrm.eof()){
                c++;
                //::std::cout << c << ::std::endl;
                readMeta(fstrm, fasta_tag, Fasta());
                read(fstrm, fasta_seq, Fasta());

                table_entry* tab = parse_fasta(fasta_seq,toCString(fasta_tag));
                table_entry* r_tab = new table_entry(*tab);
                add_entry(t.left_map, tab, 'l');
                add_entry(t.right_map, r_tab, 'r');
                //::std::cout << tab->get_short_read()->get_RNA_seq_sequence() << ::std::endl;
            }//End_while
            fstrm.close();
        }else{
            ::std::cout << "Unable to open file " << file_name << ::std::endl;
        }//End_if
    }else{
        ::std::cout << "Not a fasta file (.fa, .fas or .fasta)" << ::std::endl;
    }//End_if
    return t;
}//End_Method

