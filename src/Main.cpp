#include "build_chains.h"
#include "join_chains.h"

int main(int argc, char* argv[]){
    if(argc < 3){
        ::std::cout << ::std::endl;
        ::std::cout << "Usage: build_RNA_seq_graph <fasta_file> <option>" << ::std::endl;
        ::std::cout << "options:" << ::std::endl;
        ::std::cout << "\t 1 - Print left hash table" << ::std::endl;
        ::std::cout << "\t 2 - Print right hash table" << ::std::endl;
        ::std::cout << "\t 3 - Print unspliced RNA-seq sequences" << ::std::endl;
        ::std::cout << "\t 4 - Print spliced RNA-seq sequences" << ::std::endl;
        ::std::cout << "\t 5 - Print half spliced RNA-seq sequences" << ::std::endl;
        ::std::cout << "\t 6 - Build chains of unspliced reads with half sequence overlap" << ::std::endl;
        ::std::cout << "\t 7 - Build chains of unspliced reads with specific overlap" << ::std::endl;
        ::std::cout << "\t 8 - Merge chains built of unspliced reads with half sequence overlap" << ::std::endl;
        ::std::cout << "\t 9 - Join merged chains built of unspliced reads with half sequence overlap" << ::std::endl;
        ::std::cout << ::std::endl;
        return 1;
    }//End_If

    tables table;
    if(read_fasta(argv[1], table) == 0){
        map<unsigned long long, string> chains;
        switch(::std::atoi(argv[2])){
        case 1:
            print_hash_table(table,'l');
            break;
        case 2:
            print_hash_table(table,'r');
            break;
        case 3:
            print_unspliced(table);
            break;
        case 4:
            print_spliced(table);
            break;
        case 5:
            print_half_spliced(table);
            break;
        case 6:
            build_unspliced_chains(table);
            print_unspliced_chains(table);
            break;
        case 7:
            ::std::cout << "Insert overlap: ";
            int overlap;
            ::std::cin >> overlap;
            build_unspliced_chains(table,overlap);
            print_unspliced_chains(table,overlap);
            break;
        case 8:
            build_unspliced_chains(table);
            chains = merge_unspliced_chains(table);
            print_merged_chains(chains);
            break;
        case 9:
            build_unspliced_chains(table);
            chains = merge_unspliced_chains(table);
            link_fragment_chains(table,chains);
            break;
        default:
            ::std::cout << "Wrong Option..." << ::std::endl;
        }//End_Switch
        
        //Maps deallocation
        for(hash_map::iterator it = table.left_map.begin(); it != table.left_map.end(); ++it){
            table_entry* t = it->second.p;
            table_entry* p = t;
            while(t != NULL){
                p = p->get_l_next();
                delete t;
                t = p;
            }
        }
    }
    /*for(Map::iterator it = table.right_map.begin(); it != table.right_map.end(); ++it){
        table_entry* t = it->second.p;
        table_entry* p = t;
        while(t != NULL){
            p = p->get_r_next();
            delete t;
            t = p;
        }
        //delete it->second.p;
        }*/
}//End_Main
