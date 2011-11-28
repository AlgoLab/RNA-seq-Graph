/**
* RNA-seq-Graph
* Method for reconstructing the Isoform Graph of a gene from RNA-seq data
*
* Copyright (C) 2011 Stefano Beretta <ste.beretta(-at-)gmail.com>
*
* Distributed under the terms of the GNU Affero General Public License (AGPL)
*
*
* This file is part of RNA-seq-Graph.
*
* RNA-seq-Graph is free software: you can redistribute it and/or modify
* it under the terms of the GNU Affero General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* RNA-seq-Graph is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU Affero General Public License for more details.
*
* You should have received a copy of the GNU Affero General Public License
* along with RNA-seq-Graph.  If not, see <http://www.gnu.org/licenses/>.
**/
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <ctime>
#include <set>

#include <seqan/sequence.h>
#include <seqan/file.h>

#include "read_fasta.h"

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

/***********************************/
/* Parse Fasta Information:        */
/* Input-> DNA sequence            */
/* Output-> table entry element    */
/*                                 */
/* Notes-> called by parse_fasta   */
/***********************************/
table_entry* build_entry(String<Dna5> seq){
    table_entry* el = NULL;
    string left_seq;
    assign(left_seq,prefix(seq,length(seq)/2));
    string right_seq;
    assign(right_seq,suffix(seq,length(seq)/2));
    el = new table_entry(seq,fingerprint(left_seq),fingerprint(right_seq));
    return el;
}
/*
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
        el = new table_entry(seq, "region", "id", 0, fingerprint(left_seq), fingerprint(right_seq));
    }//End_If
    return el;
}//End_Method
*/

/***************************************/
/* Add entries in the "Hash Table":    */
/* Input-> hash table                  */
/*      -> entry to be added           */
/* Output-> entry added to the table   */
/*                                     */
/* Notes-> called by parse_fasta       */
/***************************************/
void add_entry(tables &t, table_entry* entry){
    hash_map::iterator it;
    it = t.left_map.find(entry->get_left_fingerprint());
    if(it == t.left_map.end()){
        element_table el;
        el.unspliced = 1;
        el.half_spliced = 0;
        el.p = entry;
        t.left_map[entry->get_left_fingerprint()] = el;
    }else{
        table_entry* temp = t.left_map[entry->get_left_fingerprint()].p;
        table_entry* prev;
        bool found = 0;
        do{
            if(temp->get_right_fingerprint() == entry->get_right_fingerprint()){
                found = 1;
            }
            prev = temp;
            temp = temp->get_l_next();
        }while(temp != NULL && !found);
        if(found){
            prev->increase_freq();
            delete entry;
            return;
        }else{
            t.left_map[entry->get_left_fingerprint()].unspliced = 0;
            prev->set_l_next(entry);
            entry->set_l_prev(prev);
        }//End_If
    }//End_If
    
    it = t.right_map.find(entry->get_right_fingerprint());
    if(it == t.right_map.end()){
        element_table el;
        el.unspliced = 1;
        el.half_spliced = 0;
        el.p = entry;
        t.right_map[entry->get_right_fingerprint()] = el;
    }else{
        table_entry* temp = t.right_map[entry->get_right_fingerprint()].p;
        table_entry* prev;
        do{
            prev = temp;
            temp = temp->get_r_next();
        }while(temp != NULL);
        t.right_map[entry->get_right_fingerprint()].unspliced = 0;
        prev->set_r_next(entry);
        entry->set_r_prev(prev);
    }//End_If
}//End_Method


/*****************************/
/* Try to correct reads that */
/* that are spliced due to   */
/* errors                    */
/*****************************/
void pruning(tables &table, int max_mismatch, int diff_threshold){
    //PROBLEMA:
    //quando elimino un elemento da una lista (e quindi da una tabella)
    //lo devo fare anche per l'altra lista (e quindi per l'altra tabella)
    //!!!!!
    hash_map::iterator it;
    clock_t tStart = clock();
    std::cerr << "Pruning Reads...";
    std::set< table_entry* >s;
    for(it=table.left_map.begin(); it != table.left_map.end(); it++){
        if(!(*it).second.unspliced){
            table_entry* t = (*it).second.p;
            s.insert(t);
            t = t->get_l_next();	
            while(t != NULL){
                String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
		std::set< table_entry* >::iterator spliced_it;
		bool new_s = true;
		for ( spliced_it=s.begin() ; spliced_it != s.end();){
			int diff = 0;
			for(unsigned int l=0; l<length(seq); ++l){
				if((*spliced_it)->get_short_read()->get_RNA_seq_sequence()[l] != seq[l]) diff++;
			}
			if(diff <= max_mismatch){
				if((*spliced_it)->get_frequency() < t->get_frequency() && 
				   ((*spliced_it)->get_frequency()/(double)t->get_frequency())*100 < diff_threshold){
					(*spliced_it)->get_l_prev()->set_l_next((*spliced_it)->get_l_next());
					(*spliced_it)->get_l_next()->set_l_prev((*spliced_it)->get_l_prev());
					delete(*spliced_it);
					s.erase(spliced_it++);
				}else{
					if((t->get_frequency()/(double)(*spliced_it)->get_frequency())*100 <diff_threshold){
						new_s = false;
						if(t->get_l_prev() == NULL){
							(*it).second.p = t->get_l_next();
						}else{
							t->get_l_prev()->set_l_next(t->get_l_next());
						}
						if(t->get_l_next() == NULL){
							t->get_l_prev()->set_l_next(NULL);
						}else{
							t->get_l_next()->set_l_prev(t->get_l_prev());
						}
					}
					++spliced_it;
				}
			}else{
				++spliced_it;
			}
		}
		if(new_s){
                	s.insert(t);
			t = t->get_l_next();
		}else{
			table_entry* tmp = t;
			delete(t);
			t = tmp->get_l_next();
		}
                
            }//End_While*/
	    if(s.size() == 1){
		(*it).second.unspliced = true;
	    }
        }//End_If
	s.clear();
    }//End_For
    std::cerr << "Meta' " << std::endl;
    for(it=table.right_map.begin(); it != table.right_map.end(); it++){
        if(!(*it).second.unspliced){
            table_entry* t = (*it).second.p;
            s.insert(t);
            t = t->get_r_next();	
            while(t != NULL){
                String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
		std::set< table_entry* >::iterator spliced_it;
		bool new_s = true;
		for ( spliced_it=s.begin() ; spliced_it != s.end();){
			int diff = 0;
			for(unsigned int l=0; l<length(seq); ++l){
				if((*spliced_it)->get_short_read()->get_RNA_seq_sequence()[l] != seq[l]) diff++;
			}
			if(diff <= max_mismatch){
				if((*spliced_it)->get_frequency() < t->get_frequency() && 
				   ((*spliced_it)->get_frequency()/(double)t->get_frequency())*100 < diff_threshold){
					(*spliced_it)->get_r_prev()->set_r_next((*spliced_it)->get_r_next());
					(*spliced_it)->get_r_next()->set_r_prev((*spliced_it)->get_r_prev());
					delete(*spliced_it);
					s.erase(spliced_it++);
				}else{
					if((t->get_frequency()/(double)(*spliced_it)->get_frequency())*100 <diff_threshold){
						new_s = false;
						if(t->get_r_prev() == NULL){
							(*it).second.p = t->get_r_next();
						}else{
							t->get_r_prev()->set_r_next(t->get_r_next());
						}
						if(t->get_r_next() == NULL){
							t->get_r_prev()->set_r_next(NULL);
						}else{
							t->get_r_next()->set_r_prev(t->get_r_prev());
						}
					}
					++spliced_it;
				}
			}else{
				++spliced_it;
			}
		}
		if(new_s){
                	s.insert(t);
			t = t->get_r_next();
		}else{
			table_entry* tmp = t;
			delete(t);
			t = tmp->get_r_next();
		}
            }//End_While*/
	    if(s.size() == 1){
		(*it).second.unspliced = true;
	    }
        }//End_If
	s.clear();
    }//End_For
    std::cerr << "done!" << std::endl;
    std::cerr << "Pruning took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    std::cerr << " seconds." << std::endl << std::endl;
}

/******************************/
/* Read a fasta file:         */
/* Input-> RNA-seq file name  */
/* Output-> Hash table        */
/*                            */
/* Notes-> called by the main */
/******************************/
int read_fasta(char* file_name, tables &t){
    string name = file_name;
    string extension = name.substr(name.find_last_of(".") + 1);
    if(extension == "fa" || extension == "fas" || extension == "fasta"){
        unsigned long long file_dimension = 1;
	ifstream f(file_name,ios::in|ios::binary|ios::ate);
	if(f.is_open()){
            file_dimension = f.tellg();
            f.close();
	}
	else{
            ::std::cerr << "Unable to open file " << file_name << ::std::endl;
            return 1;	
	}
        unsigned long long read_size = 0;
	int perc = 0;
        ::std::fstream fstrm;
        fstrm.open(file_name, ::std::ios_base::in | ::std::ios_base::binary);
        if(fstrm.is_open()){
	    clock_t tStart = clock();
            ::std::cerr << "Processing RNA-seq file..." << ::std::endl;
            String<char> fasta_tag;
            String<Dna5> fasta_seq;
            int c = 0;
            while(!fstrm.eof()){
                
                //::std::cerr << c << ::std::endl;
                readMeta(fstrm, fasta_tag, Fasta());
                read(fstrm, fasta_seq, Fasta());
                read_size+=(length(fasta_tag)*sizeof(char)+length(fasta_seq)*sizeof(char));
                if ((read_size/(double)file_dimension)*100 - 1 >= perc) {
                    perc++;
                    if(perc%10 == 0){
                        std::cerr << "Processed: " << perc << "%" << endl;
                    }
                }
                //Parse RNA-seq Sequence
                if(length(fasta_seq) >= READ_LEN){
		    //All the substrings of length READ_LEN
		    /*
                    for(unsigned int i = 0;i<=length(fasta_seq)-READ_LEN;i++){
                        //table_entry* tab = parse_fasta(infix(fasta_seq,i,i+READ_LEN),toCString(fasta_tag));
			table_entry* tab = build_entry(infix(fasta_seq,i,i+READ_LEN));
                        add_entry(t, tab);
			c++;
                    }*/
		    //Only the first and the last substrings of length READ_LEN
		    table_entry* tab = build_entry(infix(fasta_seq,0,READ_LEN));
                    add_entry(t, tab);
		    tab = build_entry(infix(fasta_seq,length(fasta_seq)-READ_LEN,length(fasta_seq)));
                    add_entry(t, tab);
		    c+=2;
                }else{
                    ::std::cerr << "Invalid fasta entry" << ::std::endl;
                }
                //::std::cerr << tab->get_short_read()->get_RNA_seq_sequence() << ::std::endl;
            }//End_while
            ::std::cerr << "Processing RNA-seq file...done!" << ::std::endl;
	    std::cerr << "Loaded " << c << " sequences took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
            std::cerr << " seconds." << std::endl << std::endl;
            fstrm.close();

	    //Reads pruning
	    //int max_mismatch = 1, diff_threshold = 20;
	    //pruning(t,max_mismatch, diff_threshold);
        
        }else{
            ::std::cerr << "Unable to open file " << file_name << ::std::endl;
            return 1;
        }//End_if
    }else{
	::std::cerr << "Error opening RNA-seq file " << file_name << ::std::endl;
        ::std::cerr << "Not a fasta file (.fa, .fas or .fasta)" << ::std::endl;
        return 2;
    }//End_if
    return 0;
}//End_Method

