#include <iostream>

#include <cassert>

#include "build_chains.h"

/**********************************/
/* Print spliced RNA-seq          */
/* These reads have more than one */
/* entry in each hash table       */
/* (left and right)               */
/**********************************/
void print_spliced(const tables &table){
    int c = 1;
    hash_map::const_iterator it;
    std::cerr << "Printing Spliced Reads...";
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        if(!(*it).second.unspliced){
            table_entry* t = (*it).second.p;
            while(t != NULL){
                String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2) << ::std::endl;
                c++;
                t = t->get_l_next();
            }//End_While
            ::std::cout << ::std::endl;
        }//End_If
    }//End_For

    for(it=table.right_map.begin(); it != table.right_map.end(); it++){
        if(!(*it).second.unspliced){
            table_entry* t = (*it).second.p;
            while(t != NULL){
                String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2) << ::std::endl;
                c++;
                t = t->get_r_next();
            }//End_While
            ::std::cout << ::std::endl;
        }//End_If
    }//End_For
    std::cerr << "done!" << std::endl;
}//End_Method

/**********************************/
/* Print half spliced RNA-seq     */
/* These reads have more than one */
/* entry in each hash table       */
/* (left and right)               */
/**********************************/
void print_half_spliced(const tables &table){
    hash_map::const_iterator seq_it;
    int len = length(table.left_map.begin()->second.p->get_short_read()->get_RNA_seq_sequence());
    int c = 0;
    std::cerr << "Printing Half Spliced Reads...";
    //Look for half spliced RNA-seqs
    for(seq_it=table.left_map.begin(); seq_it != table.left_map.end(); seq_it++){
        if(!(*seq_it).second.unspliced){
            int a[4] = {0, 0, 0, 0};
            int sum = 0;
            table_entry* t = (*seq_it).second.p;
            while(sum <= 1 && t != NULL){
                string str;
                assign(str,t->get_short_read()->get_RNA_seq_sequence());
                if(str[len/2] == 'A' || str[len/2] == 'a'){
                    if(a[0] == 0){
                        a[0] = 1;
                        ++sum;
                    }//End_If
                }//End_If
                if(str[len/2] == 'C' || str[len/2] == 'c'){
                    if(a[1] == 0){
                        a[1] = 1;
                        ++sum;
                    }//End_If
                }//End_If
                if(str[len/2] == 'G' || str[len/2] == 'g'){
                    if(a[2] == 0){
                        a[2] = 1;
                        ++sum;
                    }//End_If
                }//End_If
                if(str[len/2] == 'T' || str[len/2] == 't'){
                    if(a[3] == 0){
                        a[3] = 1;
                        ++sum;
                    }//End_If
                }//End_If
                t = t->get_l_next();
            }//End_While
            if(sum > 1){
                table_entry* t = (*seq_it).second.p;
                while(t != NULL){
                    c++;
                    ::std::cout << c << " " << t->get_short_read()->get_RNA_seq_sequence() << ::std::endl;
                    t = t->get_l_next();
                }
            }//End_If
        }//End_If
    }//End_For
    for(seq_it=table.right_map.begin(); seq_it != table.right_map.end(); seq_it++){
        if(!(*seq_it).second.unspliced){
            int a[4] = {0, 0, 0, 0};
            int sum = 0;
            table_entry* t = (*seq_it).second.p;
            while(sum <= 1 && t != NULL){
                string str;
                assign(str,t->get_short_read()->get_RNA_seq_sequence());
                if(str[len/2 - 1] == 'A' || str[len/2 - 1] == 'a'){
                    if(a[0] == 0){
                        a[0] = 1;
                        ++sum;
                    }//End_If
                }//End_If
                if(str[len/2 - 1] == 'C' || str[len/2 - 1] == 'c'){
                    if(a[1] == 0){
                        a[1] = 1;
                        ++sum;
                    }//End_If
                }//End_If
                if(str[len/2 - 1] == 'G' || str[len/2 - 1] == 'g'){
                    if(a[2] == 0){
                        a[2] = 1;
                        ++sum;
                    }//End_If
                }//End_If
                if(str[len/2 - 1] == 'T' || str[len/2 - 1] == 't'){
                    if(a[3] == 0){
                        a[3] = 1;
                        ++sum;
                    }//End_If
                }//End_If
                t = t->get_r_next();
            }//End_While
            if(sum > 1){
                table_entry* t = (*seq_it).second.p;
                while(t != NULL){
                    c++;
                    ::std::cout << c << " " << t->get_short_read()->get_RNA_seq_sequence() << ::std::endl;
                    t = t->get_r_next();
                }
            }//End_If
        }//End_If
    }//End_For
    std::cerr << "done!" << std::endl;
}

/********************************/
/* Print unspliced RNA-seq      */
/* These reads have exactly one */
/* entry in each hash table     */
/* (left and right)             */
/********************************/
void print_unspliced(const tables &table){
    int c = 1;
    hash_map::const_iterator it;
    std::cerr << "Printing Unspliced Reads...";
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        if((*it).second.unspliced){
            table_entry* t = (*it).second.p;
            if(table.right_map.find(t->get_right_fingerprint())->second.unspliced){
                String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2)  << ::std::endl;
                c++;
                ::std::cout << ::std::endl;
            }//End_If
        }//End_If
    }//End_For
    std::cerr << "done!" << std::endl;
}//End_Method

/*******************************/
/* Print left or right table   */
/*******************************/
void print_hash_table(const tables &table, char l_r){
    if(l_r == 'l'){
	std::cerr << "Printing Left Hash Table...";
        int c = 1;
        long sum = 0;
        hash_map::const_iterator it;
        for(it=table.left_map.begin(); it != table.left_map.end(); it++){
            table_entry* t = (*it).second.p;
            while(t != NULL){
                String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2) << ::std::endl;
                      //<< " " << t->get_frequency() << ::std::endl;
                sum += t->get_frequency();
                t = t->get_l_next();
                c++;
            }//End_While
            ::std::cout << ::std::endl;
            //c++;
        }//End_For
 	std::cerr << "done!" << std::endl;
        std::cerr << "Tot " << sum << ::std::endl;
    }else{
	std::cerr << "Printing Right Hash Table...";
        int c = 1;
        long sum = 0;
        hash_map::const_iterator it;
        for(it=table.right_map.begin(); it != table.right_map.end(); it++){
            table_entry* t = (*it).second.p;
            while(t != NULL){
                String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << c << " " << prefix(seq,length(seq)/2) << " " << suffix(seq,length(seq)/2) << ::std::endl;
                //<< " " << t->get_frequency() << ::std::endl;
                sum += t->get_frequency();
                t = t->get_r_next();
                c++;
            }//End_While
            ::std::cout << ::std::endl;
            //c++;
        }//End_For
	std::cerr << "done!" << std::endl;
        std::cerr << "Tot " << sum << ::std::endl;
    }//End_If
}//End_Method

/*************************************/
/* Build chains of unspliced RNA-seq */
/* considering half sequence overlap */
/*************************************/
void build_unspliced_chains(tables &table){
    hash_map::iterator it;
    clock_t tStart = clock();
    std::cerr << "Composing Unspliced Chains...";
    for(it = table.left_map.begin(); it != table.left_map.end(); ++it){
        if((*it).second.unspliced && table.right_map[(*it).second.p->get_right_fingerprint()].unspliced){
            table_entry* t = (*it).second.p;
            if(table.left_map.find(t->get_right_fingerprint()) != table.left_map.end()){
                if(table.left_map[t->get_right_fingerprint()].unspliced && 
                   table.right_map[table.left_map[t->get_right_fingerprint()].p->get_right_fingerprint()].unspliced &&
                   table.left_map[t->get_right_fingerprint()].p->get_chain_prev() == NULL){
                    table_entry* t_temp = table.left_map[t->get_right_fingerprint()].p;
                    t->set_chain_next(t_temp);
                    t_temp->set_chain_prev(t);
                }
            }//End_If
        }
    }//End_For
    std::cerr << "done!" << std::endl;
    std::cerr << "Chains composition took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    std::cerr << " seconds." << std::endl << std::endl;
}//End_Method

/*************************************/
/* Build chains of unspliced RNA-seq */
/* wit parametrized overlapping      */
/*************************************/
void build_unspliced_chains(tables &table, int overlap){
    hash_map::iterator it;
    clock_t tStart = clock();
    std::cerr << "Composing Unspliced Chains...";
    for(it = table.left_map.begin(); it != table.left_map.end(); ++it){
        if((*it).second.unspliced && table.right_map[(*it).second.p->get_right_fingerprint()].unspliced){
            table_entry* t = (*it).second.p;
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
                }
            }
        }
    }//End_For
    std::cerr << "done!" << std::endl;
    std::cerr << "Chains composition took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    std::cerr << " seconds." << std::endl << std::endl;
}//End_Method

/*************************************/
/* Print chains of unspliced RNA_seq */
/*************************************/
void print_unspliced_chains(const tables &table){
    int c = 1;
    hash_map::const_iterator it;
    std::cerr << "Printing Chains with Half Overlap...";
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        table_entry* t = (*it).second.p;
        if((*it).second.unspliced  && 
           table.right_map.find((*it).second.p->get_right_fingerprint())->second.unspliced && 
           t->get_chain_prev() == NULL){
            String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
            ::std::cout << c << " " << prefix(seq,length(seq)/2) << suffix(seq,length(seq)/2);
            while(t->get_chain_next() != NULL){
                t = t->get_chain_next();
                seq = t->get_short_read()->get_RNA_seq_sequence();
                std::cout << suffix(seq,length(seq)/2);
            }//End_While
            c++;
            std::cout << ::std::endl << ::std::endl;
        }//End_If
    }//End_For
    std::cerr << "done!" << std::endl;
}//End_Method

/*************************************/
/* Print chains of unspliced RNA_seq */
/* built with overlap                */
/*************************************/
void print_unspliced_chains(const tables &table, int overlap){
    int c = 1;
    hash_map::const_iterator it;
    std::cerr << "Printing Chains with Overlap " << overlap << "...";
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        table_entry* t = (*it).second.p;
        if((*it).second.unspliced &&
           table.right_map.find((*it).second.p->get_right_fingerprint())->second.unspliced && 
           t->get_chain_prev() == NULL){
            String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
            std::cout << c << " " << prefix(seq,length(seq)/2) << suffix(seq,length(seq)/2);
            while(t->get_chain_next() != NULL){
                t = t->get_chain_next();
                seq = t->get_short_read()->get_RNA_seq_sequence();
                ::std::cout << suffix(seq,length(seq)-overlap);
            }//End_While
            c++;
            std::cout << ::std::endl << ::std::endl;
        }//End_If
    }//End_For
    std::cerr << "done!" << std::endl;
}//End_Method

/******************************/
/* Merge two chains based on  */
/* their overlap              */
/******************************/
static string
merge_chains(const string& head, int offset, const string& to_be_chained){
	//std::cerr << head << "  " << offset << "  " << to_be_chained << std::endl;
	string result= "";
    if(to_be_chained.length() > head.length() - offset){
        long l = to_be_chained.length() - (head.length() - offset);
        //::std::cout << head.substr(offset) << ::std::endl;
        //::std::cout << to_be_chained.substr(0, to_be_chained.length() - l) << ::std::endl;
        if(head.substr(offset) == to_be_chained.substr(0, to_be_chained.length() - l)){ 
            //::std::cout << to_be_chained.substr(to_be_chained.length() - l) << ::std::endl;
		result= head + to_be_chained.substr(to_be_chained.length() - l);
        }//End_If
    }else{
        if(head.substr(offset, to_be_chained.length()) == to_be_chained){
            result= head;
        }
    }//End_If
    return result;
}//End_Method

/*************************************/
/* Merge previous built unspliced    */
/* chains with half sequence overlap */
/*************************************/
void merge_unspliced_chains(tables &table, map<unsigned long long, string>& chain_map){
    int c = 0;
    unsigned int rna_seq_length = length((*table.left_map.begin()).second.p->get_short_read()->get_RNA_seq_sequence());
    hash_map::const_iterator it;
    clock_t tStart = clock();
    std::cerr << "Building Unspliced Chains...";
    // 1 - Build a new table with the chains
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        table_entry* t = (*it).second.p;
        //Uncomment and substitute with the following to consider only chain longer than simple RNA-seqs
        //if((*it).second.unspliced && t->get_chain_prev() == NULL && t->get_chain_next() != NULL){
        if((*it).second.unspliced && t != NULL &&
           table.right_map.find((*it).second.p->get_right_fingerprint())->second.unspliced && 
           t->get_chain_prev() == NULL){
            
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
            }//End_While
            //::std::cout << c << " " << chain << ::std::endl;
            chain_map[t->get_left_fingerprint()] = chain;
	    //Eliminazione read utilizzati
		
	    table_entry* p = t;
            while(t != NULL){
		//table.left_map.erase(p->get_left_fingerprint());
		//table.right_map.erase(p->get_right_fingerprint());
		table.left_map[p->get_left_fingerprint()].p = NULL;
		table.right_map[p->get_right_fingerprint()].p = NULL;
                p = p->get_chain_next();
                delete t;
                t = p;
            }
            c++;
        }//End_If
    }//End_For
    hash_map::iterator iterator = table.left_map.begin();
    while(iterator != table.left_map.end()){
	if(iterator->second.p == NULL){
		table.left_map.erase(iterator++);
	}else{
		++iterator;
	}
    }
    iterator = table.right_map.begin();
    while(iterator != table.right_map.end()){
	if(iterator->second.p == NULL){
		table.right_map.erase(iterator++);
	}else{
		++iterator;
	}
    }
    std::cerr << "done!" << std::endl;
    std::cerr << "Built " << c << " chains took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    std::cerr << " seconds." << std::endl << std::endl;
    // 2 - Merge the chains built  
    
    map<unsigned long long, string>::iterator iter;
    c = 0;
    /****************************************************/
    /* Merging depth -> length of the considered prefix */
    /****************************************************/
    //unsigned int merg_depth = rna_seq_length*2;
    std::cerr << "Merging Chains...";
    tStart = clock();
    for(iter=chain_map.begin(); iter != chain_map.end(); iter++){
        //c++;
        //std::cout << c << " " << (*iter).second << ::std::endl << ::std::endl;
	/****************************************************/
	/* Merging depth -> length of the considered prefix */
	/****************************************************/
        for(unsigned int j=1; j<iter->second.length()-(rna_seq_length/2)/* && (j<=merg_depth)*/; j++){
            string segment;
            assign(segment,infix((*iter).second,j,rna_seq_length/2+j));
	    map<unsigned long long, string>::const_iterator it_segm= chain_map.find(fingerprint(segment));
            if((it_segm != chain_map.end()) && (it_segm != iter)){
                string merged = merge_chains((*iter).second, j, chain_map[fingerprint(segment)]);
                if(merged != ""){
                    assert(iter->first != fingerprint(segment));
                    (*iter).second = merged;
		    //if(merged.length() < rna_seq_length/2)
			//std::cerr << merged << " " << merged.length() << std::endl;
                    chain_map.erase(fingerprint(segment));
		    c++;
                }//End_If
            }//End_If
        }//End_For
    }//End_For
    std::cerr << "done! " << std::endl;
    std::cerr << "Merging to " << chain_map.size() << " chains took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    std::cerr << " seconds." << std::endl << std::endl;
    //std::cerr << "Fine merging\n";
}//End_Method

/**************************/
/* Print chains after the */
/* merging phase          */ 
/**************************/
void print_merged_chains(::std::map<unsigned long long, string> & chains){
    int c = 0;
    std::map<unsigned long long, string>::iterator ch_it;
    std::cerr << "Printing Chains...";
    for(ch_it = chains.begin(); ch_it != chains.end(); ++ch_it){
        c++;
        std::cout << c << " " << ch_it->second << ::std::endl << ::std::endl;
    }
    std::cerr << "done!" << std::endl;
}

