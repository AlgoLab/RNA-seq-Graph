#include "build_chains.h"

/*
void test(tables table){
    Map::iterator it;
    int c = 1;
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        table_entry* t = (*it).second.p;
        if((*it).second.unspliced && 
           table.right_map[t->get_right_fingerprint()].unspliced &&
           ((t->get_short_read()->get_RNA_seq_transcript_id() == 3 && 
            t->get_short_read()->get_RNA_seq_offset() <= 2764 &&
            t->get_short_read()->get_RNA_seq_offset() >= 2700) ||
           (t->get_short_read()->get_RNA_seq_transcript_id() == 2 && 
            t->get_short_read()->get_RNA_seq_offset() <= 2764 &&
            t->get_short_read()->get_RNA_seq_offset() >= 2700) ||
            (t->get_short_read()->get_RNA_seq_transcript_id() == 2 && 
             t->get_short_read()->get_RNA_seq_offset() <= 3071 &&
             t->get_short_read()->get_RNA_seq_offset() >= 3007))
           ){
            String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
            ::std::cout << t->get_short_read()->get_RNA_seq_transcript_id() << " " << t->get_short_read()->get_RNA_seq_offset() << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2) << " - " << c << ::std::endl;
            //::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2) << " - " << t->get_short_read()->get_RNA_seq_transcript_id() << " " << t->get_short_read()->get_RNA_seq_offset() << ::std::endl;
            c++;
        }
    }
}
*/

int main(int argc, char* argv[]){
    if(argc < 3){
        ::std::cout << ::std::endl;
        ::std::cout << "Usage: read_input <fasta_file> <option>" << ::std::endl;
        ::std::cout << "options:" << ::std::endl;
        ::std::cout << "\t 1 - Print hash table (left)" << ::std::endl;
        ::std::cout << "\t 2 - Print unspliced RNA-seq sequences" << ::std::endl;
        ::std::cout << "\t 3 - Print spliced RNA-seq sequences" << ::std::endl;
        ::std::cout << "\t 4 - Build chains of unspliced reads with half sequence overlap" << ::std::endl;
        ::std::cout << "\t 5 - Build chains of unspliced reads with specific overlap" << ::std::endl;
        ::std::cout << "\t 6 - Merge chains built of unspliced reads with half sequence overlap" << ::std::endl;
        ::std::cout << ::std::endl;
        return 1;
    }//End_If

    tables table = read_fasta(argv[1]);
    map<unsigned long long, string> chains;
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
        print_spliced(table);
        break;
    case 4:
        build_unspliced_chains(table);
        print_unspliced_chains(table);
        break;
    case 5:
        ::std::cout << "Insert overlap: ";
        int overlap;
        ::std::cin >> overlap;
        build_unspliced_chains(table,overlap);
        print_unspliced_chains(table,overlap);
        break;
    case 6:
        build_unspliced_chains(table);
        chains = merge_unspliced_chains(table);
        break;
        /* case 7:
        test(table);
        break;*/
    default:
        ::std::cout << "Wrong Option..." << ::std::endl;
    }//End_Switch
}//End_Main
