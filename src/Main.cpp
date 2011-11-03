#include "build_chains.h"
#include "join_chains.h"

/******************************/
/* Print the advanced options */
/******************************/
void show_advanced(){
    ::std::cout << "\t --debug {1-8}" << ::std::endl;
    ::std::cout << "\t\t 1 - Print left hash table" << ::std::endl;
    ::std::cout << "\t\t 2 - Print right hash table" << ::std::endl;
    ::std::cout << "\t\t 3 - Print unspliced RNA-seq sequences" << ::std::endl;
    ::std::cout << "\t\t 4 - Print spliced RNA-seq sequences" << ::std::endl;
    ::std::cout << "\t\t 5 - Print half spliced RNA-seq sequences" << ::std::endl;
    ::std::cout << "\t\t 6 - Build chains of unspliced reads with half sequence overlap" << ::std::endl;
    ::std::cout << "\t\t 7 - Build chains of unspliced reads with specific overlap" << ::std::endl;
    ::std::cout << "\t\t 8 - Merge chains built of unspliced reads with half sequence overlap" << ::std::endl;
    ::std::cout << ::std::endl;
}

/***************************/
/* Print the usage options */
/***************************/
void show_usage(bool adv){
    ::std::cout << ::std::endl;
    ::std::cout << "Usage: build_RNA_seq_graph [options] --reads <RNA-seq_file>" << ::std::endl;
    ::std::cout << "options:" << ::std::endl;
    ::std::cout << "\t --ref_level {1-5}" << ::std::endl;
    ::std::cout << "\t\t 1 - Standard Algorithm (Default option)" << ::std::endl;
    ::std::cout << "\t\t 2 - Add tiny blocks" << ::std::endl;
    ::std::cout << "\t\t 3 - Add linking edges" << ::std::endl;
    ::std::cout << "\t\t 4 - Add small blocks" << ::std::endl;
    ::std::cout << "\t\t 5 - Refine overlapping nodes" << ::std::endl;
    ::std::cout << ::std::endl;
    ::std::cout << "\t -o <graphML_out_file> (Default: std output)" << ::std::endl;
    ::std::cout << ::std::endl;
    if(!adv){
        ::std::cout << "\t --advanced (Sperimental)" << ::std::endl;
        ::std::cout << ::std::endl;
        ::std::cout << "\t --help (Print this screen)" << ::std::endl;
        ::std::cout << ::std::endl;
    }else{
        show_advanced();
    }
}


/****************/
/* Program Main */
/****************/
int main(int argc, char* argv[]){
    if(argc < 2){
        show_usage(false);
        return 1;
    }//End_If

    if(std::string(argv[1]) == "--advanced"){
        show_usage(true);
        return 0;
    }

    if(std::string(argv[1]) == "--help"){
        show_usage(false);
        return 0;
    }
    
    if(argc < 3){
        show_usage(false);
        return 1;
    }//End_If

    //Read the options
    int debug = 0;
    int ref_level = 1;
    char* graphML_out_file = NULL;
    char* reads = NULL;
    for (int i = 1; i < argc; i++) {
        if (i + 1 != argc){
            if (string(argv[i]) == "--ref_level") {
                ref_level = atoi(argv[++i]);
            } else if (string(argv[i]) == "-o") {
                graphML_out_file = argv[++i];
            } else if (string(argv[i]) == "--reads") {
                reads = argv[++i];
            } else if (string(argv[i]) == "--debug") {
                debug = atoi(argv[++i]);
            } else {
                std::cout << "Invalid arguments, please try again.\n";
                return 1;
            }
        }
    }
    if(reads == NULL){
        std::cout << "Invalid arguments, please try again.\n";
        show_usage(false);
        return 1;
    }

    tables table;
    if(read_fasta(reads, table) == 0){
        map<unsigned long long, string> chains;
        switch(debug){
        case 0:
            build_unspliced_chains(table);
            merge_unspliced_chains(table, chains);
            link_fragment_chains(table, chains, ref_level, graphML_out_file);
            break;
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
            merge_unspliced_chains(table, chains);
            print_merged_chains(chains);
            break;
        default:
            ::std::cout << "Wrong Debug Option..." << ::std::endl;
        }//End_Switch
    
        //Map deallocation
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
}//End_Main
