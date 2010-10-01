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
                String<Dna> seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2)  << ::std::endl;
                c++;
                ::std::cout << ::std::endl;
            }//End_If
        }//End_If
    }//End_For
}

/*******************************/
/* Print left or right table   */
/*******************************/
void print_hash_table(tables table,char l_r){
    if(l_r == 'l'){
        int c = 1;
        Map::iterator it;
        for(it=table.left_map.begin(); it != table.left_map.end(); it++){
            table_entry* t = (*it).second.p;
            while(t != NULL){
                String<Dna> seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2)  << ::std::endl;
                t = t->get_next();
            }//End_While
            ::std::cout << ::std::endl;
            c++;
        }//End_For
    }else{
        int c = 1;
        Map::iterator it;
        for(it=table.right_map.begin(); it != table.right_map.end(); it++){
            table_entry* t = (*it).second.p;
            while(t != NULL){
                String<Dna> seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2)  << ::std::endl;
                t = t->get_next();
            }//End_While
            ::std::cout << ::std::endl;
            c++;
        }//End_For
    }//End_If
}//End_Method

/*************************************/
/* Build chains of unspliced RNA-seq */
/*************************************/
void build_unspliced_chains(tables table){
    Map::iterator it = table.left_map.begin();
    while(it != table.left_map.end()){
        Map::iterator it_temp = it;
        if((*it).second.unspliced){
            table_entry* t = (*it).second.p;
            while(t->get_chain_next() != NULL){
                t = t->get_chain_next();
            }//End_While
            if(table.left_map.find(t->get_right_fingerprint()) != table.left_map.end()){
                // table.left_map[t->get_right_fingerprint()];
                if(table.left_map[t->get_right_fingerprint()].unspliced){
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
    }//End_For
}

/*************************************/
/* Print chains of unspliced RNA_seq */
/*************************************/
void print_unspliced_chains(tables table){
    int c = 1;
    Map::iterator it;
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        table_entry* t = (*it).second.p;
        if((*it).second.unspliced && t->get_chain_prev() == NULL){
            String<Dna> seq = t->get_short_read()->get_RNA_seq_sequence();
            ::std::cout << c << " " << prefix(seq,length(seq)/2) << suffix(seq,length(seq)/2);
            while(t->get_chain_next() != NULL){
                t = t->get_chain_next();
                seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << suffix(seq,length(seq)/2) << ::std::endl;
            }//End_While
            c++;
            ::std::cout << ::std::endl;
        }//End_If
    }//End_For
}

int main(int argc, char* argv[]){
    tables table = read_fasta(argv[1]);
    //print_hash_table(table,'l');
    //print_unspliced(table);
    build_unspliced_chains(table);
    print_unspliced_chains(table);
}
