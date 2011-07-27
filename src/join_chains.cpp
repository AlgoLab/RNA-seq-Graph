#include <vector>
#include <fstream>
#include <iostream>
#include <seqan/find.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>

#include "join_chains.h"
#include "build_chains.h"

/******************************/
/* Build and print the graph  */
/******************************/
void print_graph(::std::vector<table_entry*> links, const map<unsigned long long, string> chains,
                 map<unsigned long long, unsigned long long> mapping, char* graphML_out_file){
    map<unsigned long long, string>::const_iterator ch_iter;
    map<unsigned long long, int> graph_nodes;
    ofstream out_file;
    out_file.open("RNA-seq-graph.txt");
    long node_id = 0;
    //::std::cout << "graph: {" << ::std::endl;
    //::std::cout << "\tnode.shape\t: circle" << ::std::endl;
    //::std::cout << "\tnode.color\t: blue" << ::std::endl;
    //::std::cout << "\tnode.height\t: 80" << ::std::endl;
    //::std::cout << "\tnode.width\t: 80" << ::std::endl;

    //GraphML

    typedef boost::property<boost::vertex_name_t, int, 
        boost::property<boost::vertex_color_t, std::string> > VertexProperty;

    typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS,
         VertexProperty> Graph;

    std::vector<std::string> vertex_names;
    for(ch_iter = chains.begin(); ch_iter != chains.end(); ++ch_iter){
        vertex_names.push_back(ch_iter->second);
        ++node_id;
        graph_nodes[ch_iter->first] = node_id;
        
        //GDL output
        //::std::cout << "\t node: {" << ::std::endl;
        //::std::cout << "\t\t title: \"" << node_id << "\"" << ::std::endl;
        //::std::cout << "\t\t label: \"" << node_id << " - " << ch_iter->second.length() << "\"" << ::std::endl;
        //::std::cout << "\t\t //" << ch_iter->second << ::std::endl;
        //::std::cout << "\t}" << ::std::endl;
        //File output
        out_file << "node#" << node_id << " " << ch_iter->second << "\n"; 
    }

    //Graph Initialization
    //bool graph[node_id][node_id];
    //for(long i=0; i<node_id; ++i){
    //    for(long j=0; j<node_id; ++j){
    //        graph[i][j] = 0;
    //    }
    //}

    Graph ug(node_id);
    //Adding edges
    int num_edges = 0;
    for(unsigned int i=0; i<links.size(); ++i){
        //::std::cout << i << " " << links[i].seq << ::std::endl;
        for(int j=0; j<links[i]->size_D_link(); ++j){
            //::std::cout << "j " << j << ::std::endl;
            //::std::cout << links[i].D_delta[j] << ::std::endl;
            for(int k=0; k<links[i]->size_A_link(); ++k){
                if(graph_nodes.find(mapping[links[i]->at_D_link(j)]) != graph_nodes.end() &&
                   graph_nodes.find(mapping[links[i]->at_A_link(k)]) != graph_nodes.end()){
                    if(graph_nodes[mapping[links[i]->at_D_link(j)]] != graph_nodes[mapping[links[i]->at_A_link(k)]]){
// && 
//                       graph[graph_nodes[mapping[links[i]->at_D_link(j)]]-1][graph_nodes[mapping[links[i]->at_A_link(k)]]-1] == 0){
                        //::std::cout << "k " << k << ::std::endl;
                        //::std::cout << links[i].A_delta[k] << ::std::endl;
                        //graph[graph_nodes[mapping[links[i]->at_D_link(j)]]-1][graph_nodes[mapping[links[i]->at_A_link(k)]]-1] = 1;
                        //GDL output
                        //::std::cout << "\t edge: {" << ::std::endl;
                        //::std::cout << "\t\t source: \"" << graph_nodes[mapping[links[i]->at_D_link(j)]] << "\"" << ::std::endl;
                        //::std::cout << "\t\t target: \"" << graph_nodes[mapping[links[i]->at_A_link(k)]] << "\"" << ::std::endl;
                        //::std::cout << "\t}" << ::std::endl;
                        //File output
                        num_edges++;
                        int source = graph_nodes[mapping[links[i]->at_D_link(j)]];
                        int target = graph_nodes[mapping[links[i]->at_A_link(k)]];
                        out_file << "edge#" << num_edges << " ";
                        out_file << graph_nodes[mapping[links[i]->at_D_link(j)]] << ";";
                        out_file << graph_nodes[mapping[links[i]->at_A_link(k)]] << "\n";
                        ::boost::add_edge(source-1,target-1,ug);
                    }
                }
            }//End_For
        }//End_For
    }//End_For

    ::boost::dynamic_properties dp;
    ::boost::graph_traits<Graph>::vertex_iterator v, v_end;
    
    for(tie(v,v_end) = vertices(ug); v!=v_end;++v){
        put(::boost::vertex_name_t(), ug, *v, vertex_names[(*v)].length());
        put(::boost::vertex_color_t(), ug, *v, vertex_names[(*v)]);
    }
    
    dp.property("length", get(::boost::vertex_name_t(), ug));
    dp.property("seq", get(::boost::vertex_color_t(), ug));
    if(graphML_out_file == NULL){
        ::boost::write_graphml(::std::cout, ug, dp, true);
    }else{
        if(::std::strstr(graphML_out_file, ".graphml") == NULL){
            strcat(graphML_out_file, ".graphml");
        }
     
        ::std::ofstream out_stream;
        out_stream.open(graphML_out_file);
        ::boost::write_graphml(out_stream, ug, dp, true);
    }
    //::std::cout << "}" << ::std::endl;
    out_file.close();
}//End_Method

/**********************************/
/* Look for half spliced reads to */
/* join chains on the right       */
/**********************************/
int get_left_linked_read(string chain, tables& table, int delta){
    string r_t;
    int q = 0;
    int right_cut = 0;
    //Try to join and cut the chain on the right
    while(q <= delta){
        ::seqan::assign(r_t,::seqan::infix(chain,::seqan::length(chain)-2*delta+q,::seqan::length(chain)-delta+q));
        //::std::cout << r_t << ::std::endl;
        hash_map::iterator l_it = table.left_map.find(fingerprint(r_t));
        if(l_it != table.left_map.end()){
            if((*l_it).second.half_spliced){
                string head;
                //::std::cout << delta-q << ::std::endl;
                assign(head,::seqan::prefix(chain,delta));
                //::std::cout << "CH R " << chain << ::std::endl;
                //::std::cout << "SE L " << ::seqan::prefix((*l_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << "  " << ::seqan::suffix((*l_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << ::std::endl;
                table_entry* t = (*l_it).second.p;
                while(t != NULL){
                    t->push_D_link(fingerprint(head));
                    t = t->get_l_next();
                }
                
                if(delta-q>right_cut){
                    //::std::cout << "UNO" << ::std::endl;
                    right_cut = delta-q;
                }
            }else{
                table_entry* t = (*l_it).second.p;
                while(t != NULL){
                    hash_map::iterator r_it = table.right_map.find(t->get_right_fingerprint());
                    if((*r_it).second.half_spliced){
                        string head;
                        //::std::cout << delta-q << ::std::endl;
                        //::std::cout << "CH R " << chain << ::std::endl;
                        //::std::cout << "SE R " << ::seqan::prefix((*r_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << "  " << ::seqan::suffix((*r_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << ::std::endl;
                        assign(head,::seqan::prefix(chain,delta));
                        //::std::cout << r_t << ::std::endl;
                        table_entry* right_t = (*r_it).second.p;
                        while(right_t != NULL){
                            if(fingerprint(r_t) == right_t->get_left_fingerprint()){
                                right_t->push_D_link(fingerprint(head));
                            }
                            right_t = right_t->get_r_next();
                        }
                        
                        if(delta-q>right_cut){
                            //::std::cout << "DUE" << ::std::endl;
                            right_cut = delta-q;
                        }
                    }
                    t = t-> get_l_next();
                }
            }
        }//End_If
        ++q;
    }//End_For
    return right_cut;
}//End_Method

/**********************************/
/* Look for half spliced reads to */
/* join chains on the left        */
/**********************************/
int get_right_linked_read(string chain, tables& table, int delta){
    string l_t;
    int q = delta;
    int left_cut = 0;
    //Try to join and cut the chain on the left
    while(q >= 0){
        ::seqan::assign(l_t,::seqan::infix(chain,q,delta+q));
        //::std::cout << r_t << ::std::endl;
        hash_map::iterator r_it = table.right_map.find(fingerprint(l_t));
        if(r_it != table.right_map.end()){
            if((*r_it).second.half_spliced){
                string head;
                //::std::cout << q << ::std::endl;
                //::std::cout << "CH L " << chain << ::std::endl;
                //::std::cout << "R " << ::seqan::prefix((*r_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << "  " << ::seqan::suffix((*r_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << ::std::endl;
                assign(head,::seqan::prefix(chain,delta));
                table_entry* t = (*r_it).second.p;
                while(t != NULL){
                    t->push_A_link(fingerprint(head));
                    t = t->get_r_next();
 
                }
                if(q>left_cut){
                    //::std::cout << "TRE " << q << ::std::endl;
                    left_cut = q;
                }
            }else{
                table_entry* t = (*r_it).second.p;
                while(t != NULL){
                    hash_map::iterator l_it = table.left_map.find(t->get_left_fingerprint());
                    if((*l_it).second.half_spliced){
                        string head;
                        //::std::cout << q << ::std::endl;
                        //::std::cout << "CH L " << chain << ::std::endl;
                        //::std::cout << "L " << ::seqan::prefix((*l_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << "  " << ::seqan::suffix((*l_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << ::std::endl;
                        assign(head,::seqan::prefix(chain,delta));
                        table_entry* left_t = (*l_it).second.p;
                        while(left_t != NULL){
                            if(fingerprint(l_t) == left_t->get_right_fingerprint()){
                                left_t->push_A_link(fingerprint(head));
                            }
                            left_t = left_t->get_l_next();
                        }
                        if(q>left_cut){
                            //::std::cout << "QUA " << q << ::std::endl;
                            left_cut = q;
                        }
                    }
                    t = t->get_r_next();
                }
            }
        }//End_If
        --q;
    }//End_While
    return left_cut;
}//End_Mehtod

/*************************************/
/* Merge chains looking from the end */
/*************************************/
::std::map<unsigned long long, unsigned long long> chain_back_merging(map<unsigned long long, string>& chains, int len){
    ::std::map<unsigned long long, unsigned long long> mapping;
    ::std::multimap<unsigned long long, unsigned long long> new_chains;
    ::std::map<unsigned long long, string>::iterator chain_it;
    //unsigned int c = 0;
    queue<unsigned long long> q;
    for(chain_it = chains.begin(); chain_it != chains.end(); ++chain_it){
        mapping[chain_it->first] = chain_it->first;
        //c++;
        string tail;
        ::seqan::assign(tail,::seqan::suffix(chain_it->second,chain_it->second.length()-len));
        //::std::cout << c  << " " << tail << " " << tail.length() << ::std::endl;
        unsigned long long f_print = fingerprint(tail);
        if(new_chains.find(f_print) != new_chains.end()){
            multimap<unsigned long long, unsigned long long>::iterator it;
            pair<multimap<unsigned long long,unsigned long long>::iterator,
                multimap<unsigned long long,unsigned long long>::iterator> ret;
            ret = new_chains.equal_range(f_print);
            bool erased = false;
            for(it=ret.first; it!=ret.second; ++it){
                string dub = chains[it->second];
                if(dub.length() > chain_it->second.length()){
                    long diff = dub.length() - chain_it->second.length();
                    if(::seqan::suffix(dub,diff) == chain_it->second){
                        //::std::cout << "Merging 1 " << diff << ::std::endl;
                        //::std::cout << ::seqan::suffix(dub,diff) << ::std::endl;
                        //::std::cout << chain_it->second << ::std::endl;
                        q.push(chain_it->first);
                        erased = true;
                        mapping[chain_it->first] = chains.find(it->second)->first;
                    }
                }else{
                    long diff = chain_it->second.length() - dub.length();
                    if(::seqan::suffix(chain_it->second,diff) == dub){
                        //::std::cout << "Merging 2 " << diff << ::std::endl;
                        //::std::cout << ::seqan::suffix(chain_it->second,diff) << ::std::endl;
                        //::std::cout << dub << ::std::endl;
                        q.push(chains.find(it->second)->first);
                        erased = true;
                        mapping[chains.find(it->second)->first] = chain_it->first;
                    }
                }
            }
            if(!erased){
                new_chains.insert(pair<unsigned long long, unsigned long long>(f_print,chain_it->first));
            }
        }else{
            new_chains.insert(pair<unsigned long long, unsigned long long>(f_print,chain_it->first));
        }
    }
    while(!q.empty()){
        chains.erase(q.front());
        q.pop();
    }
    /*
    c=0;
    ::std::cout << ::std::endl;
    for(chain_it = chains.begin(); chain_it != chains.end(); ++chain_it){
        c++;
        string tail;
        ::seqan::assign(tail,::seqan::suffix(chain_it->second,chain_it->second.length()-len));
        ::std::cout << c  << " " << tail << " " << tail.length() << ::std::endl;
    }
    */
    return mapping;
}

/***************************************/
/* Merge chains looking for substrings */
/***************************************/
::std::map<unsigned long long, unsigned long long> chains_unify(map<unsigned long long, string>& chains, unsigned int len){
    ::std::map<unsigned long long, string>::iterator chain_it;
    ::std::map<unsigned long long, string>::iterator chain_it2;
    ::std::map<unsigned long long, unsigned long long> mapping;
    queue<unsigned long long> q;
    for(chain_it = chains.begin(); chain_it != chains.end(); ++chain_it){
        int max_ch_len = 0;
        mapping[chain_it->first] = chain_it->first;
        if(chain_it->second.length() > len){
            CharString p = chain_it->second;
            Pattern<CharString, ShiftOr > pattern(p);
            for(chain_it2 = chains.begin(); chain_it2 != chains.end(); ++chain_it2){
                if(chain_it != chain_it2 && chain_it->second.length() < chain_it2->second.length()){
                    CharString text = chain_it2->second;
                    Finder<CharString> finder(text);
                    int ch_len = chain_it2->second.length();
                    if(find(finder,pattern)){
                        q.push(chain_it->first);
                        if(ch_len > max_ch_len){
                            mapping[chain_it->first] = chain_it2->first;
                            max_ch_len = ch_len;
                        }
                    }
                }
            }
        } 
    }
    while(!q.empty()){
        chains.erase(q.front());
        q.pop();
    }
    return mapping;
}

/***********************************/
/* Check if the fragments cutted   */
/* during the link phase are nodes */
/* too                             */
/***********************************/
void check_cutted_frags(CharString frag, ::std::vector<table_entry*> &links, 
                        map<unsigned long long, string> &chains, unsigned int min_length){
    if(length(frag)>min_length){
        ::std::queue<int> l_link;
        ::std::queue<int> r_link;
        Pattern<CharString, ShiftOr > pattern(frag);
        for(unsigned int i=0; i<links.size(); ++i){
            CharString text = links[i]->get_short_read()->get_RNA_seq_sequence();
            Finder<CharString> finder(text);
            find(finder,pattern);
            if(beginPosition(finder) < min_length){
                //::std::cout << "L link " << i << ::std::endl;
                l_link.push(i);
            }
            if(endPosition(finder) > length(text) - min_length){
                //::std::cout << "R link" << ::std::endl;
                r_link.push(i);
            }
        }
        
        if(l_link.size() != 0 && r_link.size() != 0){
            string head;
            assign(head,frag);
            for(unsigned int z=0; z<min_length*2 - length(frag);++z){
                head.append("A");
            }
            chains[fingerprint(head)] = toCString(frag);
            
            //::std::cout << toCString(frag) << ::std::endl;
            while(!l_link.empty()){
                links[l_link.front()]->push_D_link(fingerprint(head));
                l_link.pop();
            }
            while(!r_link.empty()){
                links[r_link.front()]->push_A_link(fingerprint(head));
                r_link.pop();
            }
        }        
    }
}

void link_fragment_chains(tables& table, map<unsigned long long, string> & chains, int ref_level, char* out_file){
    ::std::vector<table_entry*> linking_reads;
    hash_map::iterator seq_it;
    int len = length(table.left_map.begin()->second.p->get_short_read()->get_RNA_seq_sequence());
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
                        //}else{
                        //dupl = 1;
                    }//End_If
                }//End_If
                if(str[len/2] == 'C' || str[len/2] == 'c'){
                    if(a[1] == 0){
                        a[1] = 1;
                        ++sum;
                        //}else{
                        //dupl = 1;
                    }//End_If
                }//End_If
                if(str[len/2] == 'G' || str[len/2] == 'g'){
                    if(a[2] == 0){
                        a[2] = 1;
                        ++sum;
                        //}else{
                        //dupl = 1;
                    }//End_If
                }//End_If
                if(str[len/2] == 'T' || str[len/2] == 't'){
                    if(a[3] == 0){
                        a[3] = 1;
                        ++sum;
                        //}else{
                        //dupl = 1;
                    }//End_If
                }//End_If
                t = t->get_l_next();
            }//End_While
            if(sum > 1){
                (*seq_it).second.half_spliced = 1;
                table_entry* t = (*seq_it).second.p;
                while(t != NULL){
                    linking_reads.push_back(t);
                    //::std::cout << t->get_short_read()->get_RNA_seq_sequence() << ::std::endl;
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
                        //}else{
                        //dupl = 1;
                    }//End_If
                }//End_If
                if(str[len/2 - 1] == 'C' || str[len/2 - 1] == 'c'){
                    if(a[1] == 0){
                        a[1] = 1;
                        ++sum;
                        //}else{
                        //dupl = 1;
                    }//End_If
                }//End_If
                if(str[len/2 - 1] == 'G' || str[len/2 - 1] == 'g'){
                    if(a[2] == 0){
                        a[2] = 1;
                        ++sum;
                        //}else{
                        //dupl = 1;
                    }//End_If
                }//End_If
                if(str[len/2 - 1] == 'T' || str[len/2 - 1] == 't'){
                    if(a[3] == 0){
                        a[3] = 1;
                        ++sum;
                        //}else{
                        //dupl = 1;
                    }//End_If
                }//End_If
                t = t->get_r_next();
            }//End_While
            if(sum > 1){
                (*seq_it).second.half_spliced = 1;
                table_entry* t = (*seq_it).second.p;
                while(t != NULL){
                    linking_reads.push_back(t);
                    //::std::cout << t->get_short_read()->get_RNA_seq_sequence() << ::std::endl;
                    t = t->get_r_next();
                }
            }//End_If
        }//End_If
    }//End_For

    //Look for linking RNA-seqs
    map<unsigned long long, string>::iterator chain_it;
    //Come prova impostiamo delta a 1/2 della lunghezza dei read
    int delta = len/2;
    string new_chain;
    ::std::cerr << "Linking Chains...";
    for(chain_it = chains.begin(); chain_it != chains.end(); ++chain_it){
        //::std::cout << "Left" << ::std::endl;
        int right_cut = get_left_linked_read(chain_it->second, table, delta);
        //::std::cout << "Right" << ::std::endl;
        int left_cut = get_right_linked_read(chain_it->second, table, delta);
        //::std::cout << left_cut << " " << right_cut << ::std::endl;
        assign(new_chain,::seqan::infix(chain_it->second,left_cut,chain_it->second.length()-right_cut));
        //::std::cout << chain_it->second << ::std::endl;
        //::std::cout << new_chain << ::std::endl << ::std::endl;

        check_cutted_frags(::seqan::prefix(chain_it->second,left_cut),linking_reads,chains,delta/2);
        check_cutted_frags(::seqan::suffix(chain_it->second,chain_it->second.length()-right_cut),linking_reads,chains,delta/2);
        //::std::cout << "Pre " << ::seqan::prefix(chain_it->second,left_cut) << ::std::endl;
        //::std::cout << "Suf " << ::seqan::suffix(chain_it->second,chain_it->second.length()-right_cut) << ::std::endl;
        chain_it->second = new_chain;
        //::std::cout << "Fine catena" << ::std::endl;
    }//End_For
    ::std::cerr << "done!" << ::std::endl;
    //print_merged_chains(chains);
    ::std::cerr << "Merging Chains...";
    ::std::map<unsigned long long, unsigned long long> mapping = chains_unify(chains,delta/2);
    ::std::cerr << "done!" << ::std::endl;
    //print_merged_chains(chains);
    //::std::map<unsigned long long, unsigned long long> mapping = chain_back_merging(chains,delta);

    switch(ref_level){
    case 1:
        ::std::cerr << "Standard Algorithm...done!" << ::std::endl;
        break;
    case 2:
        ::std::cerr << "Graph Refinement..." << ::std::endl;
        ::std::cerr << "Step 1...";
        tiny_blocks(linking_reads,chains,delta,mapping);
        ::std::cerr << "done!" << ::std::endl;
        ::std::cerr << "Graph Refinement...done!" << ::std::endl;
        break;
    case 3:
        ::std::cerr << "Graph Refinement..." << ::std::endl;
        ::std::cerr << "Step 1...";
        tiny_blocks(linking_reads,chains,delta,mapping);
        ::std::cerr << "done!" << ::std::endl;
        ::std::cerr << "Step 2...";
        add_linking_reads(linking_reads,chains,delta/2);
        ::std::cerr << "done!" << ::std::endl;
        ::std::cerr << "Graph Refinement...done!" << ::std::endl;
        break;
    case 4:
        ::std::cerr << "Graph Refinement..." << ::std::endl;
        ::std::cerr << "Step 1...";
        tiny_blocks(linking_reads,chains,delta,mapping);
        ::std::cerr << "done!" << ::std::endl;
        ::std::cerr << "Step 2...";
        add_linking_reads(linking_reads,chains,delta/2);
        ::std::cerr << "done!" << ::std::endl;
        ::std::cerr << "Step 3...";
        small_blocks(linking_reads,chains,delta,mapping);
        ::std::cerr << "done!" << ::std::endl;
        ::std::cerr << "Graph Refinement...done!" << ::std::endl;
        break;
    case 5:
        ::std::cerr << "Graph Refinement..." << ::std::endl;
        ::std::cerr << "Step 1...";
        tiny_blocks(linking_reads,chains,delta,mapping);
        ::std::cerr << "done!" << ::std::endl;
        ::std::cerr << "Step 2...";
        add_linking_reads(linking_reads,chains,delta/2);
        ::std::cerr << "done!" << ::std::endl;
        ::std::cerr << "Step 3...";
        small_blocks(linking_reads,chains,delta,mapping);
        ::std::cerr << "done!" << ::std::endl;
        ::std::cerr << "Step 4...";
        check_overlapping_nodes(linking_reads,chains,delta,mapping,5,95);
        ::std::cerr << "done!" << ::std::endl;
        ::std::cerr << "Graph Refinement...done!" << ::std::endl;
        break;
    default:
        ::std::cerr << "Wrong Refinement Option..." << ::std::endl;
        ::std::cerr << "No refinement preformed!" << ::std::endl;
    }

    //std::cerr << "Exporting Output...";
    //linking_refinement(linking_reads,chains,delta,mapping);
    print_graph(linking_reads,chains, mapping, out_file);
    //std::cerr << "done!" << std::endl;
}//End_Method

