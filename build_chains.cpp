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

/***********************************/
/* Merge previous built unspliced  */
/* chains                          */
/***********************************/
void merge_unspliced_chains(tables table){
    int c = 1;
    Map::iterator it;
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        table_entry* t = (*it).second.p;
        if((*it).second.unspliced && t->get_chain_prev() == NULL && t->get_chain_next() != NULL){
            //Chains longer than simple RNA-seq
            table_entry* t_temp = t;
            while(t_temp != NULL){
                String<Dna5> seq = t_temp->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2);
                ::std::cout << " " << t_temp->get_short_read()->get_RNA_seq_transcript_id();
                ::std::cout << " " << t_temp->get_short_read()->get_RNA_seq_offset();
                ::std::cout << ::std::endl;
                t_temp = t_temp->get_chain_next();
            }
            ::std::cout << ::std::endl;
            c++;
        }//End_If
    }//End_For
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
