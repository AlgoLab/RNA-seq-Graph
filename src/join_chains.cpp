#include <vector>

#include "join_chains.h"

void print_graph(::std::vector<delta_link> links, const map<unsigned long long, string> chains){
    map<unsigned long long, string>::const_iterator ch_iter;

    map<unsigned long long, int> graph_nodes;
    int node_id = 0;
    for(ch_iter = chains.begin(); ch_iter != chains.end(); ++ch_iter){
        ++node_id;
        graph_nodes[ch_iter->first] = node_id;
    }
    //Graph Initialization
    bool graph[node_id][node_id];
    for(int i=0; i<node_id; ++i){
        for(int j=0; j<node_id; ++j){
            graph[i][j] = 0;
        }
    }

    //Adding edges
    for(unsigned int i=0; i<links.size(); ++i){
        //::std::cout << i << " " << links[i].seq << ::std::endl;
        for(unsigned int j=0; j<links[i].D_delta.size(); ++j){
            //::std::cout << links[i].D_delta[j] << ::std::endl;
            for(unsigned int k=0; k<links[i].A_delta.size(); ++k){
                //::std::cout << links[i].A_delta[k] << ::std::endl;
                graph[ graph_nodes[links[i].A_delta[k]]-1 ][ graph_nodes[links[i].D_delta[j]]-1 ] = 1;
                graph[ graph_nodes[links[i].D_delta[j]]-1 ][ graph_nodes[links[i].A_delta[k]]-1 ] = 1;
            }//End_For
        }//End_For
        //::std::cout << ::std::endl;
    }//End_For

    //Graph Visualization
    ::std::cout << "CHAINS GRAPH" << ::std::endl << ::std::endl;
    ::std::cout << "\t";
    for(int i=0; i<node_id; ++i){
        ::std::cout << i+1 << "\t";
    }
    ::std::cout << ::std::endl << ::std::endl;
    for(int i=0; i<node_id; ++i){
        ::std::cout << i+1 << "\t";
        for(int j=0; j<node_id; ++j){
            ::std::cout << graph[i][j] << "\t";
        }
        ::std::cout << ::std::endl;
    }
    ::std::cout << ::std::endl;
}//End_Method

void get_left_linked_read(string chain, ::std::vector<delta_link> &links, int delta){
    string r_t;
    for(unsigned int i=0; i<links.size(); ++i){
        int q = 1;
        bool stop = 0;
        //Try to join or cut the chain on the right
        while(q <= delta && !stop){
            ::seqan::assign(r_t,::seqan::infix(chain,::seqan::length(chain)-2*delta+q,::seqan::length(chain)-delta+q));
            //::std::cout << r_t << ::std::endl;
            if(fingerprint(r_t) == links[i].l_fingerprint){
                string head;
                //::std::cout << delta-q << ::std::endl;
                assign(head,::seqan::prefix(chain,delta));
                links[i].D_delta.push_back(fingerprint(head));
                stop = 1;
            }//End_If
            ++q;
        }//End_While
    }//End_For
}//End_Method

void get_right_linked_read(string chain, ::std::vector<delta_link> &links, int delta){
    string l_t;
    for(unsigned int i=0; i<links.size(); ++i){
        int q = delta-1;
        bool stop = 0;
        //Try to join or cut the chain on the left
        while(q >= 0 && !stop){
            ::seqan::assign(l_t,::seqan::infix(chain,q,delta+q));
            //::std::cout << r_t << ::std::endl;
            if(fingerprint(l_t) == links[i].r_fingerprint){
                string head;
                //::std::cout << q << ::std::endl;
                assign(head,::seqan::prefix(chain,delta));
                links[i].A_delta.push_back(fingerprint(head));
                stop = 1;
            }//End_If
            --q;
        }//End_While
    }//End_For
}//End_Mehtod

::std::vector<delta_link> link_fragment_chains(const tables& table, const map<unsigned long long, string> chains){
    ::std::vector<delta_link> links;
    Map::const_iterator seq_it;
    int len = length(table.left_map.begin()->second.p->get_short_read()->get_RNA_seq_sequence());
    //Look for half spliced RNA-seqs
    for(seq_it=table.left_map.begin(); seq_it != table.left_map.end(); seq_it++){
        if(!(*seq_it).second.unspliced){
            int a[4] = {0, 0, 0, 0};
            bool dupl = 0;
            table_entry* t = (*seq_it).second.p;
            while(!dupl && t != NULL){
                string str;
                assign(str,t->get_short_read()->get_RNA_seq_sequence());
                if(str[len/2] == 'A' || str[len/2] == 'a'){
                    if(a[0] == 0){
                        a[0] = 1;
                    }else{
                        dupl = 1;
                    }//End_If
                }//End_If
                if(str[len/2] == 'C' || str[len/2] == 'c'){
                    if(a[1] == 0){
                        a[1] = 1;
                    }else{
                        dupl = 1;
                    }//End_If
                }//End_If
                if(str[len/2] == 'G' || str[len/2] == 'g'){
                    if(a[2] == 0){
                        a[2] = 1;
                    }else{
                        dupl = 1;
                    }//End_If
                }//End_If
                if(str[len/2] == 'T' || str[len/2] == 't'){
                    if(a[3] == 0){
                        a[3] = 1;
                    }else{
                        dupl = 1;
                    }//End_If
                }//End_If
                t = t->get_next();
            }//End_While
            if(!dupl){
                t = (*seq_it).second.p;
                while(t != NULL){
                    String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
                    delta_link dl;
                    dl.seq = seq;
                    dl.l_fingerprint = t->get_left_fingerprint();
                    dl.r_fingerprint = t->get_right_fingerprint();
                    links.push_back(dl);
                    t = t->get_next();
                }//End_While
            }//End_If
        }//End_If
    }//End_For
    for(seq_it=table.right_map.begin(); seq_it != table.right_map.end(); seq_it++){
        if(!(*seq_it).second.unspliced){
            int a[4] = {0, 0, 0, 0};
            bool dupl = 0;
            table_entry* t = (*seq_it).second.p;
            while(!dupl && t != NULL){
                string str;
                assign(str,t->get_short_read()->get_RNA_seq_sequence());
                if(str[len/2 - 1] == 'A' || str[len/2 - 1] == 'a'){
                    if(a[0] == 0){
                        a[0] = 1;
                    }else{
                        dupl = 1;
                    }//End_If
                }//End_If
                if(str[len/2 - 1] == 'C' || str[len/2 - 1] == 'c'){
                    if(a[1] == 0){
                        a[1] = 1;
                    }else{
                        dupl = 1;
                    }//End_If
                }//End_If
                if(str[len/2 - 1] == 'G' || str[len/2 - 1] == 'g'){
                    if(a[2] == 0){
                        a[2] = 1;
                    }else{
                        dupl = 1;
                    }//End_If
                }//End_If
                if(str[len/2 - 1] == 'T' || str[len/2 - 1] == 't'){
                    if(a[3] == 0){
                        a[3] = 1;
                    }else{
                        dupl = 1;
                    }//End_If
                }//End_If
                t = t->get_next();
            }//End_While
            if(!dupl){
                t = (*seq_it).second.p;
                while(t != NULL){
                    String<Dna5> seq = t->get_short_read()->get_RNA_seq_sequence();
                    delta_link dl;
                    dl.seq = seq;
                    dl.l_fingerprint = t->get_left_fingerprint();
                    dl.r_fingerprint = t->get_right_fingerprint();
                    links.push_back(dl);
                    t = t->get_next();
                }//End_While
            }//End_If
        }//End_If
    }//End_For

    //Look for linking RNA-seqs
    map<unsigned long long, string>::const_iterator chain_it;
    //Come prova impostiamo delta a 1/2 della lunghezza dei read
    int delta = len/2;
    for(chain_it = chains.begin(); chain_it != chains.end(); ++chain_it){
        get_left_linked_read(chain_it->second, links, delta);
        get_right_linked_read(chain_it->second, links, delta);
    }//End_For

    print_graph(links,chains);
    return links;
}//End_Method


