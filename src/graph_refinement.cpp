#include <stack>
#include <seqan/find.h>
#include "graph_refinement.h"

/***************************/
/* These functions compute */
/* the ovelap between two  */
/* strings                 */
/***************************/
unsigned int overlappedStringLength(string s1, string s2) {
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
void small_blocks(::std::vector<table_entry*> & links, map<unsigned long long, string> &chains, unsigned int len,
                  map<unsigned long long, unsigned long long>& mapping){
    map<unsigned long long, string>::iterator ch_iter;
    ::std::vector<small_frag> short_blocks;
    stack<unsigned int> s;
    for(unsigned int i=0; i<links.size(); ++i){
        for(unsigned int j=0; j<links.size(); ++j){
            if(i!=j && links[i]->size_D_link() != 0 && links[i]->size_A_link() == 0
               && links[j]->size_D_link() == 0 && links[j]->size_A_link() != 0 ){
                
                string s1,s2;
                ::seqan::assign(s1,links[i]->get_short_read()->get_RNA_seq_sequence());
                ::seqan::assign(s2,links[j]->get_short_read()->get_RNA_seq_sequence());
                
		//Overlap between s1 and s2 grater or equal than s1/2
                unsigned int overlap = overlappedStringLength(s1,s2);
                
                if(overlap > 0 && overlap <= s1.length()/2){
                    assign(s2,::seqan::suffix(s2,overlap));
                    s1.append(s2);
                    //::std::cout << s1 << " " << s1.length() << ::std::endl;
                    small_frag f;
                    f.frag_links.D_chain = i;
                    f.frag_links.A_chain = j;
                    f.frag = ::seqan::infix(s1,len, s1.length() - len);
                    //::std::cout << overlap << " " << f.frag << " " << length(f.frag) << ::std::endl;
                    short_blocks.push_back(f);
                }
            }
        }
    }
    for(unsigned int i=0; i<short_blocks.size(); ++i){
        bool sub_seq = false;
        for(unsigned int k=0; k<short_blocks.size(); ++k){
            if(short_blocks[i].frag == short_blocks[k].frag && i>k){
                links_pair erased_links;
                erased_links.D_chain = short_blocks[i].frag_links.D_chain;
                erased_links.A_chain = short_blocks[i].frag_links.A_chain;
                short_blocks[k].other_links.push_back(erased_links);
                sub_seq = true;
            }
            if(i!=k && (::seqan::length(short_blocks[i].frag)) > (::seqan::length(short_blocks[k].frag))){
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
        //::std::cout << short_blocks[i].frag << " " << length(short_blocks[i].frag) << ::std::endl;
        string ch;
        assign(ch,::seqan::prefix(short_blocks[i].frag,len));
        //::std::cout << ch << " " << ch.length() << ::std::endl;
        if(chains.find(fingerprint(ch)) == chains.end()){
            chains[fingerprint(ch)] = ::seqan::toCString(short_blocks[i].frag);
            //::std::cout << ::seqan::toCString(short_blocks[i].frag) << " " << length(short_blocks[i].frag) << ::std::endl;
            mapping[fingerprint(ch)] = fingerprint(ch);
            links[short_blocks[i].frag_links.D_chain]->push_A_link(fingerprint(ch));
            links[short_blocks[i].frag_links.A_chain]->push_D_link(fingerprint(ch));
        }
        for(unsigned int j=0; j<short_blocks[i].other_links.size(); ++j){
            links[short_blocks[i].other_links[j].D_chain]->push_A_link(fingerprint(ch));
            links[short_blocks[i].other_links[j].A_chain]->push_D_link(fingerprint(ch));
        }
    }
}

/*******************************/
/* Find fragments with length  */
/* less than l/2 and link them */
/* to existing chains          */
/*******************************/
void tiny_blocks(::std::vector<table_entry*> & links, map<unsigned long long, string> &chains, int len,
                  map<unsigned long long, unsigned long long>& mapping){
    //map<unsigned long long, int> graph_nodes;
    map<unsigned long long, string>::iterator ch_iter;
    ::std::vector<small_frag> short_blocks;
    stack<unsigned int> s;

    ::std::cerr << "Link size: " << links.size() << ::std::endl;

    for(unsigned int i=0; i<links.size(); ++i){
        for(unsigned int j=0; j<links.size(); ++j){
            if(i!=j && links[i]->size_D_link() != 0 && links[j]->size_A_link() != 0){
               //&& links[j]->size_D_link() == 0 && links[j]->size_A_link() != 0 ){
                
                string s1,s2;
                ::seqan::assign(s1,links[i]->get_short_read()->get_RNA_seq_sequence());
                ::seqan::assign(s2,links[j]->get_short_read()->get_RNA_seq_sequence());
                
		//Overlap between s1 and s2 grater or equal than s1/2
                unsigned int overlap = overlappedStringLength(s1,s2);
                
                if(overlap > s1.length()/2 && overlap < s1.length()-6){
                    //::std::cout << s1 << ::std::endl;
                    //::std::cout << s2 << ::std::endl;
                    assign(s2,::seqan::suffix(s2,overlap));
                    s1.append(s2);
                    small_frag f;
                    f.frag_links.D_chain = i;
                    f.frag_links.A_chain = j;
                    f.frag = ::seqan::infix(s1,len, s1.length() - len);
                    //::std::cout << overlap << " " << f.frag << " " << length(f.frag) << ::std::endl;
                    short_blocks.push_back(f);
                }
            }
        }
        ::std::cerr << i << " ";
    }
    ::std::cerr << "Fine Primo Ciclo" << ::std::endl;
    ::std::cerr << "Short Blocks (initial) size: " << short_blocks.size() << ::std::endl;
    for(unsigned int i=0; i<short_blocks.size(); ++i){
        bool sub_seq = false;
        for(unsigned int k=0; k<short_blocks.size(); ++k){
            if(short_blocks[i].frag == short_blocks[k].frag && i<k){
                links_pair erased_links;
                erased_links.D_chain = short_blocks[i].frag_links.D_chain;
                erased_links.A_chain = short_blocks[i].frag_links.A_chain;
                short_blocks[k].other_links.push_back(erased_links);
                sub_seq = true;
            }
            if(i!=k && (::seqan::length(short_blocks[i].frag)) < (::seqan::length(short_blocks[k].frag))){
                Finder<CharString> finder(short_blocks[k].frag);
                Pattern<CharString, ShiftAnd> pattern(short_blocks[i].frag);
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
    ::std::cerr << "Fine Secondo Ciclo" << ::std::endl;
    ::std::cerr << "Short Blocks (final) size: " << short_blocks.size() << ::std::endl;
    for(unsigned int i=0; i<short_blocks.size(); ++i){//Start_For_1
        bool new_frag = true;
        Pattern<CharString, ShiftAnd> pattern(short_blocks[i].frag);
        for(ch_iter = chains.begin(); ch_iter != chains.end(); ++ch_iter){//Start_For_2
            if(ch_iter->second.length() > length(short_blocks[i].frag)){
                CharString ch_text = ch_iter->second;
                Finder<CharString> finder(ch_text);
                
                if(find(finder,pattern) ||
                   overlappedStringLength(ch_iter->second,toCString(short_blocks[i].frag))>length(short_blocks[i].frag)/2 ||
                   overlappedStringLength(toCString(short_blocks[i].frag),ch_iter->second)>length(short_blocks[i].frag)/2){//Start_If_3
                    new_frag = false;
                }//End_If_3
            }
        }//End_For_2
        
        
        if(new_frag){//Start_If_7
            //::std::cout << short_blocks[i].frag << " " << length(short_blocks[i].frag) << ::std::endl;
            string ch = "";
            for(unsigned int z = 0; z<len-length(short_blocks[i].frag); ++z){//Start_If_4
                ch.append("A");
            }//End_If_4
            ch.append(toCString(short_blocks[i].frag));
            //::std::cout << ch << " " << ch.length() << ::std::endl;
            if(chains.find(fingerprint(ch)) == chains.end()){//Start_If_5
                chains[fingerprint(ch)] = ::seqan::toCString(short_blocks[i].frag);
                //::std::cout << ::seqan::toCString(short_blocks[i].frag) <<" "<< length(short_blocks[i].frag)<<::std::endl;
                mapping[fingerprint(ch)] = fingerprint(ch);
                links[short_blocks[i].frag_links.D_chain]->push_A_link(fingerprint(ch));
                links[short_blocks[i].frag_links.A_chain]->push_D_link(fingerprint(ch));
                //::std::cout <<  links[short_blocks[i].frag_links.D_chain]->get_short_read()->get_RNA_seq_sequence() << ::std::endl;
                //::std::cout <<  links[short_blocks[i].frag_links.A_chain]->get_short_read()->get_RNA_seq_sequence() << ::std::endl;
                for(unsigned int j=0; j<short_blocks[i].other_links.size(); ++j){//Start_For_6
                    links[short_blocks[i].other_links[j].D_chain]->push_A_link(fingerprint(ch));
                    links[short_blocks[i].other_links[j].A_chain]->push_D_link(fingerprint(ch));
                }//End_For_6
            }//End_If_5
        }//End_If_7
    }//End_For_1
}//End_Method

/*******************************/
/* Adds links between existing */
/* nodes based on the match of */
/* a minimum overlap sequence  */
/*******************************/
void add_linking_reads(::std::vector<table_entry*> & links, const map<unsigned long long, string> chains, 
                       unsigned int min_overlap){
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
            assign(read_tail,
                   ::seqan::suffix(links[i]->get_short_read()->get_RNA_seq_sequence(),
                                   ::seqan::length(links[i]->get_short_read()->get_RNA_seq_sequence())-min_overlap));
            //::std::cout << read_tail << ::std::endl;
            for(ch_iter = chains.begin(); ch_iter != chains.end(); ++ch_iter){
                unsigned int q = 0;
                while(q < min_overlap*2){
                    if(ch_iter->second.length() > min_overlap*2){
                        string chain_head;
                        assign(chain_head,::seqan::infix(ch_iter->second,q,q+min_overlap));
                        if(read_tail == chain_head){
                            links[i]->push_A_link(ch_iter->first);
                            //::std::cout << "Aggiunto A " << graph_nodes[ch_iter->first] << ::std::endl;
                        }
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
                unsigned int q = 0;
                while(q < min_overlap*2){
                    if(ch_iter->second.length() > min_overlap*2){
                        string chain_tail;
                        assign(chain_tail,::seqan::infix(ch_iter->second,ch_iter->second.length()-q-min_overlap,
                                                         ch_iter->second.length()-q));
                        if(read_head == chain_tail){
                            links[i]->push_D_link(ch_iter->first);
                            //::std::cout << "Aggiunto D" << ::std::endl;
                        }
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

/*********************************/
/* Chack if some nodes could be  */
/* broken intoto two subnodes in */
/* order to link other existing  */
/* (and unlinked) nodes          */
/*********************************/
void linking_refinement(::std::vector<table_entry*> & links, map<unsigned long long, string> & chains, unsigned int len,
                        ::std::map<unsigned long long, unsigned long long> & mapping){
    for(unsigned int i=0; i<links.size(); ++i){
        //Linkato solo a dx
        if(links[i]->size_D_link() == 0 && links[i]->size_A_link() != 0){
            //::std::cout << "D link" << ::std::endl;
            CharString p = ::seqan::prefix(links[i]->get_short_read()->get_RNA_seq_sequence(),len);
            Pattern<CharString, ShiftOr > pattern(p);
            ::std::map<unsigned long long, string>::iterator chain_it;
            ::std::set<unsigned long long> modif_chains;
            for(chain_it = chains.begin(); chain_it != chains.end(); ++chain_it){ 
                
                CharString text = chain_it->second;
                Finder<CharString> finder(text);
                
                if(modif_chains.find(chain_it->first) == modif_chains.end() && find(finder,pattern)){
                    links[i]->push_D_link(chain_it->first);
                    if(chain_it->second.length()- endPosition(finder) > len){
                        //::std::cout << "D " << (i+1) << " " << beginPosition(finder) << ::std::endl;
                        CharString pre = ::seqan::prefix(chain_it->second, beginPosition(finder) + len);
                        string str_pre = ::seqan::toCString(pre);
                        CharString suf = ::seqan::suffix(chain_it->second, beginPosition(finder) + len);
                        string str_suf = ::seqan::toCString(suf);
                        //::std::cout << chain_it->second << " - " << chain_it->second.length() << ::std::endl;
                        //Sono sicuro che sia > len dato che la estraggo da un prefisso
                        //di lunghezza len...
                        chains[chain_it->first] = str_pre;
                        //::std::cout << str_pre << " - " << str_pre.length() << ::std::endl;
                        modif_chains.insert(chain_it->first);
                        //...ma il suffissopotrebbe essere piu' corto di len
                        string head;
                        if(str_suf.length() >= len){
                            head = ::seqan::toCString(::seqan::prefix(suf,len));
                            chains[fingerprint(head)] = str_suf;
                            mapping[fingerprint(head)] = fingerprint(head);
                        }else{
                            head = str_suf;
                            for(unsigned int z=0; z<len-str_suf.length();++z){
                                head.append("A");
                            }
                            chains[fingerprint(head)] = str_suf;
                            mapping[fingerprint(head)] = fingerprint(head);
                        }
                        //::std::cout << str_suf << " - " << str_suf.length() << ::std::endl << ::std::endl;
                        modif_chains.insert(fingerprint(head));
                        for(unsigned int z=0; z<links.size();++z){
                            for(int k=0; k<links[z]->size_D_link();++k){
                                if(links[z]->at_D_link(k) == chain_it->first){
                                    links[z]->at_D_link(k) = fingerprint(head);
                                }
                            }
                        }
                        //Aggiungere un link tra le due catene create
                        CharString l_part = chains[chain_it->first];
                        string new_link = ::seqan::toCString(::seqan::suffix(l_part,length(l_part) - len));
                        unsigned long long f_l = fingerprint(new_link);
                        new_link.append(head);
                        table_entry* t_new = new table_entry(new_link,f_l,fingerprint(head));
                        t_new->push_D_link(chain_it->first);
                        t_new->push_A_link(fingerprint(head));
                        links.push_back(t_new);
                    }
                }
            }
        }
        
        //Linkato solo a sx
        if(links[i]->size_A_link() == 0 && links[i]->size_D_link() != 0){
            //::std::cout << "A link" << ::std::endl;
            CharString p = ::seqan::suffix(links[i]->get_short_read()->get_RNA_seq_sequence(),len);
            Pattern<CharString, ShiftOr > pattern(p);
            ::std::map<unsigned long long, string>::iterator chain_it;
            ::std::set<unsigned long long> modif_chains;
            for(chain_it = chains.begin(); chain_it != chains.end(); ++chain_it){ 
                CharString text = chain_it->second;
                Finder<CharString> finder(text);
            
                if(modif_chains.find(chain_it->first) == modif_chains.end() && find(finder,pattern)){
                //if(find(finder,pattern)){
                    //::std::cout << "1 - if " << beginPosition(finder) << " " << endPosition(finder) << ::std::endl;
                    if(beginPosition(finder) == 0){
                        links[i]->push_A_link(chain_it->first);
                    }
                    if(endPosition(finder) > len){
                        //::std::cout << "A " << (i+1) << " " << beginPosition(finder) << ::std::endl;
                        CharString pre = ::seqan::prefix(chain_it->second, beginPosition(finder) + len);
                        string str_pre = ::seqan::toCString(pre);
                        CharString suf = ::seqan::suffix(chain_it->second, beginPosition(finder) + len);
                        string str_suf = ::seqan::toCString(suf);
                        chains[chain_it->first] = str_pre;
                        //::std::cout << str_pre << " - " << str_pre.length() << ::std::endl;
                        modif_chains.insert(chain_it->first);
                        string head;
                        if(str_suf.length() >= len){
                            head = ::seqan::toCString(::seqan::prefix(suf,len));
                            chains[fingerprint(head)] = str_suf;
                            mapping[fingerprint(head)] = fingerprint(head);
                        }else{
                            head = str_suf;
                            for(unsigned int z=0; z<len-str_suf.length();++z){
                                head.append("A");
                            }
                            chains[fingerprint(head)] = str_suf;
                            mapping[fingerprint(head)] = fingerprint(head);
                        }
                        //::std::cout << str_suf << " - " << str_suf.length() << ::std::endl << ::std::endl;
                        modif_chains.insert(fingerprint(head));
                        for(unsigned int z=0; z<links.size();++z){
                            for(int k=0; k<links[z]->size_D_link();++k){
                                if(links[z]->at_D_link(k) == chain_it->first){
                                    links[z]->at_D_link(k) = fingerprint(head);
                                }
                            }
                        }
                        //Aggiungere un link tra le due catene create
                        CharString l_part = chains[chain_it->first];
                        string new_link = ::seqan::toCString(::seqan::suffix(l_part,length(l_part) - len));
                        unsigned long long f_l = fingerprint(new_link);
                        new_link.append(head);
                        table_entry* t_new = new table_entry(new_link,f_l,fingerprint(head));
                        t_new->push_D_link(chain_it->first);
                        t_new->push_A_link(fingerprint(head));
                        links.push_back(t_new);

                        links[i]->push_A_link(fingerprint(head));
                    }
                }
            }
        }
    }
    //::std::cout << chains.size() << ::std::endl;
}

/********************************/
/* Chack if there is an overlap */
/* between existing nodes that  */
/* could be a node itself       */
/********************************/
void check_overlapping_nodes(::std::vector<table_entry*> & links, map<unsigned long long, string> & chains, int len,
                             ::std::map<unsigned long long, unsigned long long>& mapping, unsigned int min_overlap,
                             int ov_perc){
    ::std::map<unsigned long long, string>::iterator chain_it;
    ::std::map<unsigned long long, string>::iterator chain_it_2;
    ::std::vector<small_frag> short_blocks;
    stack<unsigned int> s;
    queue<unsigned long long> q;
    for(chain_it = chains.begin(); chain_it != chains.end(); ++chain_it){
        for(chain_it_2 = chains.begin(); chain_it_2 != chains.end(); ++chain_it_2){
            unsigned int ov = overlappedStringLength(chain_it->second,chain_it_2->second);
            if(chain_it != chain_it_2 && ov < (ov_perc*chain_it->second.length())/100 &&
               (ov_perc*ov < chain_it_2->second.length())/100 && ov > min_overlap){
                bool new_node = false;
                CharString pat_text=prefix(chain_it_2->second,ov);
                //::std::cout << chain_it->second << ::std::endl;
                //::std::cout << chain_it_2->second << ::std::endl;
                //::std::cout << ov << ::std::endl;
                Pattern<CharString, ShiftAnd> pattern(pat_text);
                for(unsigned int i=0; i<links.size();++i){
                    CharString link_read = links[i]->get_short_read()->get_RNA_seq_sequence();
                    Finder<CharString> finder(link_read);
                    if(find(finder,pattern) && (
                       prefix(link_read,beginPosition(finder)) == infix(chain_it->second,chain_it->second.length()-ov-beginPosition(finder),chain_it->second.length()-ov) ||
                       suffix(link_read,length(link_read) - endPosition(finder)) == infix(chain_it_2->second,ov,ov+endPosition(finder)))){
                        //::std::cout << link_read << ::std::endl;
                        //::std::cout << prefix(link_read,beginPosition(finder)) << ::std::endl;
                        //::std::cout << infix(chain_it->second,chain_it->second.length()-ov-beginPosition(finder),chain_it->second.length()-ov) << ::std::endl;
                        //::std::cout << suffix(link_read,length(link_read) - endPosition(finder)) << ::std::endl;
                        //::std::cout << infix(chain_it_2->second,ov,ov+endPosition(finder)) << ::std::endl;
                      
                        new_node = true;
                    }
                }
                if(new_node){
                    small_frag f;
                    f.frag_links.D_chain = chain_it->first;
                    f.frag_links.A_chain = chain_it_2->first;
                    f.frag = prefix(chain_it_2->second,ov);
                    short_blocks.push_back(f);
                }
            }else{
                if(chain_it != chain_it_2 && ov>=(ov_perc*chain_it->second.length())/100){
                    //::std::cout << "Chain_it sub-node of Chain_it_2" << ::std::endl;
                    //::std::cout << "Chain_it " << chain_it->second << ::std::endl;
                    //::std::cout << "Chain_it_2 " << chain_it_2->second << ::std::endl;
                    //::std::cout << ov << ::std::endl;
                    q.push(chain_it->first);
                }else{
                    if(chain_it != chain_it_2 && ov>=(ov_perc*chain_it_2->second.length())/100){
                        //::std::cout << "Chain_it_2 sub-node of Chain_it" << ::std::endl;
                        //::std::cout << "Chain_it " << chain_it->second << ::std::endl;
                        //::std::cout << "Chain_it_2 " <<chain_it_2->second << ::std::endl;
                        //::std::cout << ov << ::std::endl;
                        q.push(chain_it_2->first);
                    }
                }
            }
        }
    }

    for(unsigned int i=0; i<short_blocks.size(); ++i){
        bool sub_seq = false;
        for(unsigned int k=0; k<short_blocks.size(); ++k){
            if(short_blocks[i].frag == short_blocks[k].frag && i<k){
                links_pair erased_links;
                erased_links.D_chain = short_blocks[i].frag_links.D_chain;
                erased_links.A_chain = short_blocks[i].frag_links.A_chain;
                short_blocks[k].other_links.push_back(erased_links);
                sub_seq = true;
            }
            if(i!=k && (::seqan::length(short_blocks[i].frag)) < (::seqan::length(short_blocks[k].frag))){
                Finder<CharString> finder(short_blocks[k].frag);
                Pattern<CharString, ShiftAnd> pattern(short_blocks[i].frag);
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
    while(!q.empty()){
        chains.erase(q.front());
        q.pop();
    }

    for(unsigned int i=0; i<short_blocks.size(); ++i){
        //::std::cout << short_blocks[i].frag << " " << length(short_blocks[i].frag) << ::std::endl; 
        string ch = "";
        for(unsigned int z = 0; z<len-length(short_blocks[i].frag); ++z){
            ch.append("A");
        }
        ch.append(toCString(short_blocks[i].frag));
        //if(chains.find(fingerprint(ch)) == chains.end()){//Start_If_5
            //chains[fingerprint(ch)] = ::seqan::toCString(short_blocks[i].frag);
            //::std::cout << ::seqan::toCString(short_blocks[i].frag) <<" "<< length(short_blocks[i].frag)<<::std::endl;
            //mapping[fingerprint(ch)] = fingerprint(ch);
            //Add the first link
            string first_half;
            assign(first_half,prefix(chains[short_blocks[i].frag_links.D_chain],len));
            string new_link_1 = first_half;
            new_link_1.append(ch);
            table_entry* link_1 = new table_entry(new_link_1,fingerprint(first_half),fingerprint(ch));
            link_1->push_D_link(short_blocks[i].frag_links.D_chain);
            link_1->push_A_link(short_blocks[i].frag_links.A_chain);
            links.push_back(link_1);
            /*
            //Add the second link
            string second_half;
            assign(second_half,prefix(chains[short_blocks[i].frag_links.A_chain],len));
            string new_link_2 = ch;
            new_link_2.append(second_half);
            table_entry* link_2 = new table_entry(new_link_2,fingerprint(ch),fingerprint(second_half));
            link_2->push_D_link(short_blocks[i].frag_links.D_chain);
            link_2->push_A_link(short_blocks[i].frag_links.A_chain);
            links.push_back(link_2);
            */
            //::std::cout<<links[short_blocks[i].frag_links.D_chain]->get_short_read()->get_RNA_seq_sequence()<<::std::endl;
            //::std::cout<<links[short_blocks[i].frag_links.A_chain]->get_short_read()->get_RNA_seq_sequence()<<::std::endl;

            for(unsigned int j=0; j<short_blocks[i].other_links.size(); ++j){//Start_For_6
                string second_half;
                assign(first_half,prefix(chains[short_blocks[i].other_links[j].D_chain],len));
                string new_link_2 = second_half;
                new_link_2.append(ch);
                table_entry* link_2 = new table_entry(new_link_2,fingerprint(second_half),fingerprint(ch));
                link_2->push_D_link(short_blocks[i].other_links[j].D_chain);
                link_2->push_A_link(short_blocks[i].other_links[j].A_chain);
                links.push_back(link_1);
            }//End_For_6
            //}//End_If_5
    }
}
