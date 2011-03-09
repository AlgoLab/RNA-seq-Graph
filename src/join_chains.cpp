#include <vector>
#include <fstream>
#include <iostream>
#include <queue>
#include <stack>
#include <seqan/find.h>

#include "join_chains.h"

/******************************/
/* Build and print the graph  */
/******************************/
void print_graph(::std::vector<table_entry*> links, const map<unsigned long long, string> chains){
    map<unsigned long long, string>::const_iterator ch_iter;
    map<unsigned long long, int> graph_nodes;
    ofstream out_file;
    out_file.open("RNA-seq-graph.txt");
    int node_id = 0;
    ::std::cout << "graph: {" << ::std::endl;
    ::std::cout << "\tnode.shape\t: circle" << ::std::endl;
    ::std::cout << "\tnode.color\t: blue" << ::std::endl;
    ::std::cout << "\tnode.height\t: 80" << ::std::endl;
    ::std::cout << "\tnode.width\t: 80" << ::std::endl;
    for(ch_iter = chains.begin(); ch_iter != chains.end(); ++ch_iter){
        ++node_id;
        graph_nodes[ch_iter->first] = node_id;
        //GDL output
        ::std::cout << "\t node: {" << ::std::endl;
        ::std::cout << "\t\t title: \"" << node_id << "\"" << ::std::endl;
        ::std::cout << "\t\t label: \"" << node_id << " - " << ch_iter->second.length() << "\"" << ::std::endl;
        ::std::cout << "\t\t //" << ch_iter->second << ::std::endl;
        ::std::cout << "\t}" << ::std::endl;
        //File output
        out_file << "node#" << node_id << " " << ch_iter->second << "\n"; 
    }
    //Graph Initialization
    bool graph[node_id][node_id];
    for(int i=0; i<node_id; ++i){
        for(int j=0; j<node_id; ++j){
            graph[i][j] = 0;
        }
    }

    //Adding edges
    int num_edges = 0;
    for(unsigned int i=0; i<links.size(); ++i){
        //::std::cout << i << " " << links[i].seq << ::std::endl;
        for(int j=0; j<links[i]->size_D_link(); ++j){
            //::std::cout << "j " << j << ::std::endl;
            //::std::cout << links[i].D_delta[j] << ::std::endl;
            for(int k=0; k<links[i]->size_A_link(); ++k){
                if(graph_nodes.find(links[i]->at_D_link(j)) != graph_nodes.end() &&
                   graph_nodes.find(links[i]->at_A_link(k)) != graph_nodes.end()){
                    if(graph_nodes[links[i]->at_D_link(j)] != graph_nodes[links[i]->at_A_link(k)] && 
                       graph[graph_nodes[links[i]->at_D_link(j)]-1][graph_nodes[links[i]->at_A_link(k)]-1] == 0){
                        //::std::cout << "k " << k << ::std::endl;
                        //::std::cout << links[i].A_delta[k] << ::std::endl;
                        graph[graph_nodes[links[i]->at_D_link(j)]-1][graph_nodes[links[i]->at_A_link(k)]-1] = 1;
                        //GDL output
                        ::std::cout << "\t edge: {" << ::std::endl;
                        ::std::cout << "\t\t source: \"" << graph_nodes[links[i]->at_D_link(j)] << "\"" << ::std::endl;
                        ::std::cout << "\t\t target: \"" << graph_nodes[links[i]->at_A_link(k)] << "\"" << ::std::endl;
                        ::std::cout << "\t}" << ::std::endl;
                        //File output
                        num_edges++;
                        out_file << "edge#" << num_edges << " ";
                        out_file << graph_nodes[links[i]->at_D_link(j)] << ";";
                        out_file << graph_nodes[links[i]->at_A_link(k)] << "\n";
                    }
                }
            }//End_For
        }//End_For
    }//End_For
    ::std::cout << "}" << ::std::endl;
    out_file.close();
}//End_Method

void add_linking_reads(::std::vector<table_entry*> & links, const map<unsigned long long, string> chains, int min_overlap){
    map<unsigned long long, string>::const_iterator ch_iter;
    map<unsigned long long, int> graph_nodes;
    //int c = 0;
    int node_id = 0;
    for(ch_iter = chains.begin(); ch_iter != chains.end(); ++ch_iter){
        ++node_id;
        graph_nodes[ch_iter->first] = node_id;
    }

    for(unsigned int i=0; i<links.size(); ++i){
        if(links[i]->size_D_link() != 0 && links[i]->size_A_link() == 0){
            string read_tail;
            assign(read_tail,::seqan::suffix(links[i]->get_short_read()->get_RNA_seq_sequence(),
                                             ::seqan::length(links[i]->get_short_read()->get_RNA_seq_sequence())-min_overlap));
            //::std::cout << read_tail << ::std::endl;
            for(ch_iter = chains.begin(); ch_iter != chains.end(); ++ch_iter){
                int q = 0;
                while(q < min_overlap*2){
                    string chain_head;
                    assign(chain_head,::seqan::infix(ch_iter->second,q,q+min_overlap));
                    if(read_tail == chain_head){
                        links[i]->push_A_link(ch_iter->first);
                        //::std::cout << "Aggiunto A " << graph_nodes[ch_iter->first] << ::std::endl;
                    }
                    ++q;
                }
            }
            /*
            c++;
            ::std::cout << c << " - " << links[i]->get_short_read()->get_RNA_seq_sequence() << ::std::endl;
            ::std::cout << links[i]->size_D_link() << " " << links[i]->size_A_link() << ::std::endl;
            for(int j=0; j<links[i]->size_D_link(); ++j){
                if(graph_nodes.find(links[i]->at_D_link(j)) != graph_nodes.end()){
                    ::std::cout << "D " << graph_nodes[links[i]->at_D_link(j)] << ::std::endl;
                }
            }
            */
        }
        if(links[i]->size_A_link() != 0 && links[i]->size_D_link() == 0){
            string read_head;
            assign(read_head,::seqan::prefix(links[i]->get_short_read()->get_RNA_seq_sequence(),min_overlap));
            //::std::cout << read_head << ::std::endl;
            for(ch_iter = chains.begin(); ch_iter != chains.end(); ++ch_iter){
                int q = 0;
                while(q < min_overlap*2){
                    string chain_tail;
                    assign(chain_tail,::seqan::infix(ch_iter->second,ch_iter->second.length()-q-min_overlap,
                                                     ch_iter->second.length()-q));
                    if(read_head == chain_tail){
                        links[i]->push_D_link(ch_iter->first);
                        //::std::cout << "Aggiunto D" << ::std::endl;
                    }
                    ++q;
                }
            }
            /*
            c++;
            ::std::cout << c << " - " << links[i]->get_short_read()->get_RNA_seq_sequence() << ::std::endl;
            ::std::cout << links[i]->size_D_link() << " " << links[i]->size_A_link() << ::std::endl;
            for(int k=0; k<links[i]->size_A_link(); ++k){
                if(graph_nodes.find(links[i]->at_A_link(k)) != graph_nodes.end()){
                    ::std::cout << "A " << graph_nodes[links[i]->at_A_link(k)] << ::std::endl;
                }
            }
            */
        }
    }
}

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
                //::std::cout << "SE L " << ::seqan::prefix((*l_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << "  " << ::seqan::suffix((*l_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << ::std::endl << ::std::endl;
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
                        //::std::cout << ::std::endl << delta-q << ::std::endl;
                        //::std::cout << "CH R " << chain << ::std::endl;
                        //::std::cout << "SE R " << ::seqan::prefix((*r_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << "  " << ::seqan::suffix((*r_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << ::std::endl << ::std::endl;
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
                //::std::cout << "SE R " << ::seqan::prefix((*r_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << "  " << ::seqan::suffix((*r_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << ::std::endl << ::std::endl;
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
                        //::std::cout << "SE L " << ::seqan::prefix((*l_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << "  " << ::seqan::suffix((*l_it).second.p->get_short_read()->get_RNA_seq_sequence(),delta) << ::std::endl << ::std::endl;
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
void chain_back_merging(map<unsigned long long, string>& chains, int len){
    multimap<unsigned long long, unsigned long long> new_chains;
    map<unsigned long long, string>::iterator chain_it;
    //unsigned int c = 0;
    queue<unsigned long long> q;
    for(chain_it = chains.begin(); chain_it != chains.end(); ++chain_it){
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
                    }
                }else{
                    long diff = chain_it->second.length() - dub.length();
                    if(::seqan::suffix(chain_it->second,diff) == dub){
                        //::std::cout << "Merging 2 " << diff << ::std::endl;
                        //::std::cout << ::seqan::suffix(chain_it->second,diff) << ::std::endl;
                        //::std::cout << dub << ::std::endl;
                        q.push(chains.find(it->second)->first);
                        erased = true;
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
}

/***************************************/
/* Merge chains looking for substrings */
/***************************************/
void chains_unify(map<unsigned long long, string>& chains){
    map<unsigned long long, string>::iterator chain_it;
    map<unsigned long long, string>::iterator chain_it2;
    queue<unsigned long long> q;
    for(chain_it = chains.begin(); chain_it != chains.end(); ++chain_it){
        for(chain_it2 = chains.begin(); chain_it2 != chains.end(); ++chain_it2){
            if(chain_it != chain_it2 && chain_it->second.length() < chain_it2->second.length()){
                CharString text = chain_it2->second;
                Finder<CharString> finder(text);
                Pattern<CharString, Horspool> pattern(chain_it->second);
                if(find(finder,pattern)){
                    q.push(chain_it->first);
                }
            }
        } 
    }
    while(!q.empty()){
        chains.erase(q.front());
        q.pop();
    }
}

int overlappedStringLength(string s1, string s2) {
    //Trim s1 so it isn't longer than s2
    if (s1.length() > s2.length()) s1 = s1.substr(s1.length() - s2.length());

    int *T = computeBackTrackTable(s2); //O(n)
    unsigned int m = 0;
    int i = 0;
    while (m + i < s1.length()) {
        if (s2[i] == s1[m + i]) {
            i += 1;
            //<-- removed the return case here, because |s1| <= |s2|
        } else {
            m += i - T[i];
            if (i > 0) i = T[i];
        }
    }
    delete[] T;
    
    return i; //<-- changed the return here to return characters matched
}

int* computeBackTrackTable(string s) {
    int *T = new int[s.length()];
    int cnd = 0;
    T[0] = -1;
    T[1] = 0;
    unsigned int pos = 2;
    while (pos < s.length()) {
        if (s[pos - 1] == s[cnd]) {
            T[pos] = cnd + 1;
            pos += 1;
            cnd += 1;
        } else if (cnd > 0) {
            cnd = T[cnd];
        } else {
            T[pos] = 0;
            pos += 1;
        }
    }
    return T;
}

/*******************************/
/* Find fragments with length  */
/* between l/2 and l and link  */
/* them to existing chains     */
/*******************************/
void small_blocks(::std::vector<table_entry*> & links, map<unsigned long long, string> &chains, int len){
    //map<unsigned long long, int> graph_nodes;
    map<unsigned long long, string>::iterator ch_iter;
    ::std::vector<small_frag> short_blocks;
    stack<unsigned int> s;
    for(unsigned int i=0; i<links.size(); ++i){
        for(unsigned int j=0; j<links.size(); ++j){
            if(i!=j && links[i]->size_D_link() != 0 && links[i]->size_A_link() == 0
               && links[j]->size_D_link() == 0 && links[j]->size_A_link() != 0 ){
                small_frag f;
                string s1,s2;
                ::seqan::assign(s1,links[i]->get_short_read()->get_RNA_seq_sequence());
                ::seqan::assign(s2,links[j]->get_short_read()->get_RNA_seq_sequence());
                
                unsigned int ov = s1.length()-overlappedStringLength(s1,s2);
                if(overlappedStringLength(s1,s2)!=0 && ov <= s1.length()/2){
                    assign(s2,::seqan::suffix(s2,overlappedStringLength(s1,s2)));
                    s1.append(s2);
                    f.frag_links.D_chain = i;
                    f.frag_links.A_chain = j;
                    f.frag = ::seqan::infix(s1,len,s1.length() - len);
                    short_blocks.push_back(f);
                }
            }
        }
    }
    
    for(unsigned int i=0; i<short_blocks.size(); ++i){
        bool sub_seq = false;
        for(unsigned int k=0; k<short_blocks.size(); ++k){
            if(i!=k && (::seqan::length(short_blocks[i].frag))>(::seqan::length(short_blocks[k].frag))){
                Finder<CharString> finder(short_blocks[i].frag);
                Pattern<CharString, ShiftAnd> pattern(short_blocks[k].frag);
                if(find(finder,pattern)){
                    links_pair erased_links;
                    erased_links.D_chain = short_blocks[i].frag_links.D_chain;
                    erased_links.A_chain = short_blocks[i].frag_links.A_chain;
                    //::std::cout << i << k << " - " << beginPosition(finder) << " " << endPosition(finder) << ::std::endl;
                    short_blocks[k].other_links.push_back(erased_links);
                    sub_seq = true;
                }
            }
        }
        if(sub_seq){
            s.push(i);
        }
    }
    
    while(!s.empty()){
        short_blocks.erase(short_blocks.begin()+s.top());
        s.pop();
    }
    
    for(unsigned int i=0; i<short_blocks.size(); ++i){
        string new_chain = "";
        string ch;
        assign(ch,short_blocks[i].frag);
        for (unsigned int k=0; k<len-length(short_blocks[i].frag);++k){
            new_chain.append("A");
        }
        new_chain.append(ch);
        if(chains.find(fingerprint(new_chain)) == chains.end()){
            chains[fingerprint(new_chain)] = ch;
            links[short_blocks[i].frag_links.D_chain]->push_A_link(fingerprint(new_chain));
            links[short_blocks[i].frag_links.A_chain]->push_D_link(fingerprint(new_chain));
        }
        for(unsigned int j=0; j<short_blocks[i].other_links.size(); ++j){
            links[short_blocks[i].other_links[j].D_chain]->push_A_link(fingerprint(new_chain));
            links[short_blocks[i].other_links[j].A_chain]->push_D_link(fingerprint(new_chain));
        }
    }
}

void link_fragment_chains(tables& table, map<unsigned long long, string> chains){
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
    for(chain_it = chains.begin(); chain_it != chains.end(); ++chain_it){
        //::std::cout << "Left" << ::std::endl;
        int right_cut = get_left_linked_read(chain_it->second, table, delta);
        //::std::cout << "Right" << ::std::endl;
        int left_cut = get_right_linked_read(chain_it->second, table, delta);
        //::std::cout << left_cut << " " << right_cut << ::std::endl;
        assign(new_chain,::seqan::infix(chain_it->second,left_cut,chain_it->second.length()-right_cut));
        //::std::cout << chain_it->second << ::std::endl;
        //::std::cout << ::seqan::prefix(chain_it->second,left_cut) << " " << new_chain << " ";
        //::std::cout << ::seqan::suffix(chain_it->second,chain_it->second.length()-right_cut) << ::std::endl << ::std::endl;
        //::std::cout << left_cut << " " << new_chain << " " << right_cut << ::std::endl << ::std::endl;
        chain_it->second = new_chain;
    }//End_For
    chains_unify(chains);
    //chain_back_merging(chains,delta);
    add_linking_reads(linking_reads,chains,delta/2);
    small_blocks(linking_reads,chains,delta);
    print_graph(linking_reads,chains);
}//End_Method


