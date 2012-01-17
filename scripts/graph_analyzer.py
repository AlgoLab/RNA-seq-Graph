import os,shutil
import re
import glob

def DFS (G, v):
	yield v
	visited = set ([v])
	S = G[v][0] | G[v][1]
	while S:
		w = S.pop()
		if w not in visited:
			yield w
			visited.add (w)
#			S.update (G[w][0] | G[w][1])
			S.update ([ x for x in G[w][0] if x not in visited ])
			S.update ([ x for x in G[w][1] if x not in visited ])

def graph_analyzer(graphml_file,threshold):
        nodes = {}
        sequences = {}
        (filename,extension) = os.path.splitext(graphml_file)
        seq_file = open(filename+'_cc',mode='w',encoding='utf-8')
        stat_file = open(filename+'_stats',mode='w',encoding='utf-8')

        with open(graphml_file,encoding='utf-8') as g_file:
                node_info = re.compile(r'\s+<node id="(\S+)">\s+$')
                length_info = re.compile(r'\s+<data key="key0">(\d+)</data>\s+$')
                seq_info = re.compile(r'\s+<data key="key1">(\S+)</data>\s+$')
                edge_info = re.compile(r'\s+<edge id="\S+" source="(\S+)" target="(\S+)">\s+$')
                end_node = re.compile(r'\s+</node>\s+$')
                for g_line in g_file:
                        edge = edge_info.match(g_line)
                        if(edge):
                                nodes[edge.group(1)][0].add(edge.group(2))
                                nodes[edge.group(2)][1].add(edge.group(1))
                                continue

                        n=node_info.match(g_line)
                        if(n):
                                node_id = n.group(1)
                                nodes[node_id]=[set(),set(),0]

                        l=length_info.match(g_line)
                        if(l):
                                length = int(l.group(1))
                                nodes[node_id][2]=length

                        s = seq_info.match(g_line)
                        if(s):
                                seq = s.group(1)
                                sequences[node_id] = seq

        isolated = set( [ vid for (vid,vinfo) in nodes.items() if not vinfo[0] and not vinfo[1] ] )
       
        stat_file.write("#Node_id\tLength\n")
        cc_num = 0
        for vid in isolated:
                stat_file.write("#{}\t{}\n".format(vid, nodes[vid][2]) )
                if nodes[vid][2] > int(threshold):
                        cc_num = cc_num + 1
                        seq_file.write("CC{}_NODE1:\t{}\n".format(cc_num,sequences[vid]) )
        stat_file.write("\n")
	
        vertexes = set(nodes.keys()) - isolated

        stat_file.write("@Node_id\tIn_deg\tOut_deg\n")			
        stat_file.write("\n".join(("@{}\t{}\t{}".format(vid,len(nodes[vid][1]),len(nodes[vid][0])) for vid in vertexes)))
        stat_file.write("\n")

        nodes_visited = set();
        stat_file.write("&CC_Nodes\tCC_Tot_Len\n")
        for vsrc in vertexes:
                if vsrc in nodes_visited:
                        continue
                cc_len = 0	
                cc_nodes = 0
                cc_set = set()
                for v in DFS (nodes, vsrc):
                        cc_len = cc_len + nodes[v][2]
                        cc_nodes = cc_nodes + 1
                        nodes_visited.add(v)
                        cc_set.add(v)
                stat_file.write("&"+str(cc_nodes)+"\t"+str(cc_len)+"\n")
                if cc_len > int(threshold):
                        cc_num = cc_num + 1
                        cc_num_node = 0
                        for vid in cc_set:
                                cc_num_node += 1
                                seq_file.write("CC{}_NODE{}:\t{}\n".format(cc_num,cc_num_node,sequences[vid]))
                        #seq_file.write("\n".join("CC{}_NODE:\t{}".format(cc_num,sequences[vid]) for vid in cc_set))
                        #seq_file.write("\n")
        seq_file.close()
        stat_file.close()

def main(arg1,arg2):
        graph_analyzer(arg1,arg2)
    
    
if __name__ == "__main__":
        import sys
        if len(sys.argv) == 3:
                main(sys.argv[1],sys.argv[2])
        else:
                print("Usage: python3 graph_analyzer.py <GraphML_file> <Length_Threshold>")
