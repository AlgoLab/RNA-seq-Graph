#include <iostream>

#include "read_fasta.h"

/********************************/
/* Print unspliced RNA-seq      */
/* These reads have exactly one */
/* entry in each hash table     */
/* (left and right)             */
/********************************/
void print_unspliced(tables table){
    int c = 1;
    Map::iterator it;
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        if((*it).second.unspliced){
            table_entry* t = (*it).second.p;
            if(table.right_map[t->get_right_fingerprint()].unspliced){
                String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2)  << ::std::endl;
                c++;
                ::std::cout << ::std::endl;
            }//End_If
        }//End_If
    }//End_For
}//End_Method

/*******************************/
/* Print left or right table   */
/*******************************/
void print_hash_table(tables table,char l_r){
    if(l_r == 'l'){
        int c = 1;
        long sum = 0;
        Map::iterator it;
        for(it=table.left_map.begin(); it != table.left_map.end(); it++){
            table_entry* t = (*it).second.p;
            while(t != NULL){
                String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2) << " " << t->get_frequency() << ::std::endl;
                sum += t->get_frequency();
                t = t->get_next();
                c++;
            }//End_While
            ::std::cout << ::std::endl;
            //c++;
        }//End_For
        ::std::cout << "Tot " << sum << ::std::endl;
    }else{
        int c = 1;
        long sum = 0;
        Map::iterator it;
        for(it=table.right_map.begin(); it != table.right_map.end(); it++){
            table_entry* t = (*it).second.p;
            while(t != NULL){
                String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2) << " " << t->get_frequency() << ::std::endl;
                sum += t->get_frequency();
                t = t->get_next();
                c++;
            }//End_While
            ::std::cout << ::std::endl;
            //c++;
        }//End_For
        ::std::cout << "Tot " << sum << ::std::endl;
    }//End_If
}//End_Method

/*************************************/
/* Build chains of unspliced RNA-seq */
/* considering half sequence overlap */
/*************************************/
void build_unspliced_chains(tables table){
    Map::iterator it = table.left_map.begin();
    while(it != table.left_map.end()){
        Map::iterator it_temp = it;
        if((*it).second.unspliced && table.right_map[(*it).second.p->get_right_fingerprint()].unspliced){
            table_entry* t = (*it).second.p;
            while(t->get_chain_next() != NULL){
                t = t->get_chain_next();
            }//End_While
            if(table.left_map.find(t->get_right_fingerprint()) != table.left_map.end()){
                if(table.left_map[t->get_right_fingerprint()].unspliced && 
                   table.right_map[table.left_map[t->get_right_fingerprint()].p->get_right_fingerprint()].unspliced){
                    table_entry* t_temp = table.left_map[t->get_right_fingerprint()].p;
                    t->set_chain_next(t_temp);
                    t_temp->set_chain_prev(t);
                }else{
                    it = ++it_temp;
                }//End_If
            }else{
                it = ++it_temp;
            }//End_If
        }else{
            it = ++it_temp;
        }//End_If
    }//End_While
}//End_Method

/*************************************/
/* Build chains of unspliced RNA-seq */
/* wit parametrized overlapping      */
/*************************************/
void build_unspliced_chains(tables table, int overlap){
    Map::iterator it = table.left_map.begin();
    while(it != table.left_map.end()){
        Map::iterator it_temp = it;
        if((*it).second.unspliced && table.right_map[(*it).second.p->get_right_fingerprint()].unspliced){
            table_entry* t = (*it).second.p;
            while(t->get_chain_next() != NULL){
                t = t->get_chain_next();
            }//End_While
            String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
            string segment;
            assign(segment,infix(seq,overlap,length(seq)/2+overlap));
            unsigned long long fprint = fingerprint(segment);
    
            if(table.left_map.find(fprint) != table.left_map.end()){
                if(table.left_map[fprint].unspliced &&
                   table.right_map[table.left_map[fprint].p->get_right_fingerprint()].unspliced){
                    table_entry* t_temp = table.left_map[fprint].p;
                    t->set_chain_next(t_temp);
                    t_temp->set_chain_prev(t);
                }else{
                    it = ++it_temp;
                }//End_If
            }else{
                it = ++it_temp;
            }//End_If
        }else{
            it = ++it_temp;
        }//End_If
    }//End_While
}//End_Method

/*************************************/
/* Print chains of unspliced RNA_seq */
/*************************************/
void print_unspliced_chains(tables table){
    int c = 1;
    Map::iterator it;
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        table_entry* t = (*it).second.p;
        if((*it).second.unspliced  && 
           table.right_map[(*it).second.p->get_right_fingerprint()].unspliced && 
           t->get_chain_prev() == NULL){
            String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
            ::std::cout << c << " " << prefix(seq,length(seq)/2) << suffix(seq,length(seq)/2);
            while(t->get_chain_next() != NULL){
                t = t->get_chain_next();
                seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << suffix(seq,length(seq)/2);
            }//End_While
            c++;
            ::std::cout << ::std::endl << ::std::endl;
        }//End_If
    }//End_For
}//End_Method

/*************************************/
/* Print chains of unspliced RNA_seq */
/* built with overlap                */
/*************************************/
void print_unspliced_chains(tables table, int overlap){
    int c = 1;
    Map::iterator it;
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        table_entry* t = (*it).second.p;
        if((*it).second.unspliced &&
           table.right_map[(*it).second.p->get_right_fingerprint()].unspliced && 
           t->get_chain_prev() == NULL){
            String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
            ::std::cout << c << " " << prefix(seq,length(seq)/2) << suffix(seq,length(seq)/2);
            while(t->get_chain_next() != NULL){
                t = t->get_chain_next();
                seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << suffix(seq,length(seq)-overlap);
            }//End_While
            c++;
            ::std::cout << ::std::endl << ::std::endl;
        }//End_If
    }//End_For
}//End_Method

/******************************/
/* Merge two chains based on  */
/* their overlap              */
/******************************/
string merge_chains(string head, int offset, string to_be_chained){
    if(to_be_chained.length() > head.length() - offset){
        long l = to_be_chained.length() - (head.length() - offset);
        //::std::cout << head.substr(offset) << ::std::endl;
        //::std::cout << to_be_chained.substr(0, to_be_chained.length() - l) << ::std::endl;
        if(head.substr(offset) == to_be_chained.substr(0, to_be_chained.length() - l)){ 
            //::std::cout << to_be_chained.substr(to_be_chained.length() - l) << ::std::endl;
            head.append(to_be_chained.substr(to_be_chained.length() - l));
        }else{
            return "";
        }
    }
    return head;
}

/*************************************/
/* Merge previous built unspliced    */
/* chains with half sequence overlap */
/*************************************/
void merge_unspliced_chains(tables table){
    int c = 1;
    map<unsigned long long, string> chain_map;
    Map::iterator it;
    // 1 - Build a new table with the chains
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        table_entry* t = (*it).second.p;
        //Uncomment and substitute with the following to consider only chain longer than simple RNA-seqs
        //if((*it).second.unspliced && t->get_chain_prev() == NULL && t->get_chain_next() != NULL){
        if((*it).second.unspliced && t->get_chain_prev() == NULL){
            //Chains longer than simple RNA-seq
            table_entry* t_temp = t;
            String<Dna5> seq = t_temp->get_short_read()->get_RNA_seq_sequence();
            string chain;
            assign(chain,seq);
            t_temp = t_temp->get_chain_next();
            while(t_temp != NULL){
                seq = t_temp->get_short_read()->get_RNA_seq_sequence();
                append(chain,suffix(seq,length(seq)/2));
                //::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2);
                //::std::cout << " " << t_temp->get_short_read()->get_RNA_seq_transcript_id();
                //::std::cout << " " << t_temp->get_short_read()->get_RNA_seq_offset();
                //::std::cout << ::std::endl;
                t_temp = t_temp->get_chain_next();
            }
            //::std::cout << c << " " << chain << ::std::endl;
            chain_map[t->get_left_fingerprint()] = chain;
            c++;
        }//End_If
    }//End_For
    // 2 - Merge the chains built
    unsigned int rna_seq_length = length((*table.left_map.begin()).second.p->get_short_read()->get_RNA_seq_sequence());
    map<unsigned long long, string>::iterator iter;
    c = 0;
    for(iter=chain_map.begin(); iter != chain_map.end(); iter++){
        //c++;
        //::std::cout << c << " " << (*iter).second << ::std::endl << ::std::endl;
        for(unsigned int j=1; j<rna_seq_length/2; j++){
            string segment;
            assign(segment,infix((*iter).second,j,rna_seq_length/2+j));
            if(chain_map.find(fingerprint(segment)) != chain_map.end()){
                string merged = merge_chains((*iter).second, j, chain_map[fingerprint(segment)]);
                if(merged != ""){
                    (*iter).second = merged;
                    chain_map.erase(fingerprint(segment));
                }
            }
        }
    }
    for(iter=chain_map.begin(); iter != chain_map.end(); iter++){
        c++;
        ::std::cout << c << " " << (*iter).second << " - " << (*iter).second.length() << ::std::endl << ::std::endl;
    }
}//End_Method

int main(int argc, char* argv[]){
    if(argc < 3){
        ::std::cout << ::std::endl;
        ::std::cout << "Usage: read_input <fasta_file> <option>" << ::std::endl;
        ::std::cout << "options:" << ::std::endl;
        ::std::cout << "\t 1 - Print hash table (left)" << ::std::endl;
        ::std::cout << "\t 2 - Print unspliced RNA-seq sequences" << ::std::endl;
        ::std::cout << "\t 3 - Build chains of unspliced reads with half sequence overlap" << ::std::endl;
        ::std::cout << "\t 4 - Build chains of unspliced reads with specific overlap" << ::std::endl;
        ::std::cout << "\t 5 - Merge chains built of unspliced reads with half sequence overlap" << ::std::endl;
        ::std::cout << ::std::endl;
        return 1;
    }//End_If

    tables table = read_fasta(argv[1]);
    switch(::std::atoi(argv[2])){
    case 1:
        ::std::cout << "Select left or right table (l/r): ";
        char l_r;
        ::std::cin >> l_r;
        print_hash_table(table,l_r);
        break;
    case 2:
        print_unspliced(table);
        break;
    case 3:
        build_unspliced_chains(table);
        print_unspliced_chains(table);
        break;
    case 4:
        ::std::cout << "Insert overlap: ";
        int overlap;
        ::std::cin >> overlap;
        build_unspliced_chains(table,overlap);
        print_unspliced_chains(table,overlap);
        break;
    case 5:
        build_unspliced_chains(table);
        merge_unspliced_chains(table);
        break;
    default:
        ::std::cout << "Wrong Option..." << ::std::endl;
    }//End_Switch
}//End_Main
