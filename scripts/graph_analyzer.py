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

def graph_analyzer(graphml_file):
	nodes = {}
	with open(graphml_file,encoding='utf-8') as g_file:
		node_info = re.compile(r'\s+<node id="(\S+)">\s+$')
		length_info = re.compile(r'\s+<data key="key0">(\d+)</data>\s+$')
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

	isolated = set( [ vid for (vid,vinfo) in nodes.items() if not vinfo[0] and not vinfo[1] ] )

	print("#Node_id\tLength")
	for vid in isolated:		
		print("#{}\t{}".format(vid, nodes[vid][2]) )
	print()
	
	vertexes = set(nodes.keys()) - isolated

	print("@Node_id\tIn_deg\tOut_deg")			
	print("\n".join( ( "@{}\t{}\t{}".format(vid, len(nodes[vid][1]), len(nodes[vid][0])) for vid in vertexes ) ) )
	print()

	nodes_visited = set();
	print("&CC_Nodes\tCC_Tot_Len")
	for vsrc in vertexes:
		if vsrc in nodes_visited:
			continue
		cc_len = 0	
		cc_nodes = 0
		for v in DFS (nodes, vsrc):
			cc_len = cc_len + nodes[v][2]
			cc_nodes = cc_nodes + 1
			nodes_visited.add(v)
		print("&"+str(cc_nodes),"\t",cc_len)


def main(arg):
	graph_analyzer(arg)
    
    
if __name__ == "__main__":
	import sys
	main(sys.argv[1])
