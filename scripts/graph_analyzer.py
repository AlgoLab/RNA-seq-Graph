import os,shutil
import re
import glob

def graph_analyzer(graphml_file):
	nodes = {}
	with open(graphml_file,encoding='utf-8') as g_file:
		for g_line in g_file:
			node_info = re.compile(r'\s+<node id="(\S+)">\s+$')
			n=node_info.search(g_line)
			if(n):
				node_id = n.group(1)

			length_info = re.compile(r'\s+<data key="key0">(\d+)</data>\s+$')
			l=length_info.search(g_line)
			if(l):
				length = int(l.group(1))

			end_node = re.compile(r'\s+</node>\s+$')
			if(end_node.search(g_line)):
				nodes[node_id] = length
				print(node_id," ",nodes[node_id])

			edge_info = re.compile(r'\s+<edge id="(\S+)" source="(\S+)" target="(\S+)">\s+$')
			edge = edge_info.search(g_line)
			if(edge):
				if edge.group(2) in nodes:
					del nodes[edge.group(2)]
				if edge.group(3) in nodes:
					del nodes[edge.group(3)]
		
		for x in nodes.keys():
			print("Node: ", x, "\tLength:", nodes[x])
		print()

def main(arg):
	graph_analyzer(arg)
    
    
if __name__ == "__main__":
	import sys
	main(sys.argv[1])
