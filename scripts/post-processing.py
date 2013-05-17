#! /usr/bin/env python2.7

from optparse import OptionParser
import sys
import re

def main():
    parser = OptionParser(description = "Graph post-processing.", 
             usage = "%prog -G <GDL_graph_file> -R <PS_read_file> [ -T <TXT_graph_file> ]")
    parser.add_option("-G", metavar = "<GDL_graph_file>", help = "Graph file in GDL format.")
    parser.add_option("-R", metavar = "<PS_read_file>", help = "Perfectly spliced read file.")
    parser.add_option("-T", metavar = "<TXT_graph_file>", help = "Graph file in TXT format.")
    (options, args) = parser.parse_args()
    
    gdl_file = options.G
    read_file = options.R
    txt_file = options.T

    if not (gdl_file and read_file and txt_file):
        print "Error: missing argument(s)."
        return 

    if not (gdl_file.endswith(".gdl")):
        print "GDL file error. Expected .gdl"
        return

    with open(read_file, "r") as reads, open(txt_file, "r") as txt:
        out_txt_file = ".".join(txt_file.split(".")[0:-1]) + "_post." + "".join(txt_file.rsplit(".")[-1])
        out_gdl_file = ".".join(gdl_file.split(".")[0:-1]) + "_post." + "".join(gdl_file.rsplit(".")[-1])
        ss = []
        for r in reads:
            ss.append(r.rstrip().split(" ")[1])
            #print r.rstrip().split(" ")[1]
        i = len(ss)
        #print i
        snp1 = set()
        snp2 = set()
        snp = set()
        for j in range(0,i):
            for k in range(j+1,i):
                if ( sum( [1 for c1,c2 in zip(ss[j], ss[k]) if c1 != c2] ) == 1):
                    if ( ss[j][0:len(ss[j])/2] == ss[k][0:len(ss[k])/2]):
                        snp1.add((ss[j], ss[k]))
                        #print ss[j]
                        #print ss[k]
                        #print "1"
                    else:
                        snp2.add((ss[j], ss[k]))
                        #print ss[j]
                        #print ss[k]
                        #print "2"
        del ss
        for lseq1,lseq2 in snp1:
            ls1 = lseq1[1:]
            ls2 = lseq2[1:]
            for rseq1,rseq2 in snp2:
                rs1 = rseq1[:-1]
                rs2 = rseq2[:-1]
                if (ls1 == rs1 and ls2 == rs2):
                    #print lseq1
                    #print " {0}".format(rseq1)
                    #print lseq2
                    #print " {0}".format(rseq2)
                    snp.add((lseq1[0:len(lseq1)/2],rseq1[len(rseq1)/2:],lseq1[len(lseq1)/2],rseq2[len(rseq2)/2-1]))
                    snp.add((lseq2[0:len(lseq2)/2],rseq2[len(rseq2)/2:],rseq1[len(rseq1)/2-1],lseq2[len(lseq2)/2]))
                if (ls1 == rs2 and ls2 == rs1):
                    #print lseq1
                    #print " {0}".format(rseq2)
                    #print lseq2
                    #print " {0}".format(rseq1)
                    snp.add((lseq1[0:len(lseq1)/2],rseq2[len(rseq2)/2:],lseq1[len(lseq1)/2],rseq1[len(rseq1)/2-1]))
                    snp.add((lseq2[0:len(lseq2)/2],rseq1[len(rseq1)/2:],rseq2[len(rseq2)/2-1],lseq2[len(lseq2)/2]))
            #print "--------------"
        snp1.clear()
        snp2.clear()
        #print len(snp)
        
        graph_n = {}
        graph_e = {}
        node_re = re.compile(r'^node\#([0-9]+)')
        edge_re = re.compile(r'^edge\#([0-9]+)')
        for n in txt:
            nseq = n.rstrip().split(" ")
            node = re.search(node_re, nseq[0])
            edge = re.search(edge_re, nseq[0])
            seq = nseq[1]
            if node:
                graph_n[node.group(1)] = seq
                #print int(node.group(1))
            if edge:
                ed = nseq[1].split(';')
                graph_e[edge.group(1)] = (ed[0],ed[1])

        for lh,rh,c1,c2 in snp:
            #print "{0}[{2},{3}]{1}".format(lh,rh,c1,c2)
            nlh = 0
            nrh = 0
            for n in graph_n:
                #print lh
                #print s[-len(lh):]
                if (lh == graph_n[n][-len(lh):]):
                    nlh = n
                    #print "L: {0}".format(n)
                if (rh == graph_n[n][0:len(lh)]):
                    nrh = n
                    #print "R: {0}".format(n)
            if (nlh != 0 and nrh != 0):
                elh = 0
                erh = 0
                for e in graph_e:
                    if (graph_e[e][0] == nlh):
                        elh = 1
                    if (graph_e[e][1] == nrh):
                        erh = 1
                if (elh == 0 and erh == 0):
                    new_seq = graph_n[nlh] + "[{0}{1}]".format(c1,c2) + graph_n[nrh]
                    graph_n[nlh] = new_seq
                    for e in graph_e:
                        if (graph_e[e][0] == nrh):
                            graph_e[e] = (nlh, graph_e[e][1])
                    del graph_n[nrh]

        sort = {n:i for i,n in enumerate(graph_n)}
        with open(out_txt_file, "w") as out_txt, open(out_gdl_file, "w") as out_gdl:
            out_gdl.write("graph: {\n")
            out_gdl.write("\tnode.shape\t: circle\n")
            out_gdl.write("\tnode.height\t: 80\n")
            out_gdl.write("\tnode.width\t: 80\n")
            for n in graph_n:
                out_txt.write("node#{0} {1}\n".format(sort[n]+1, graph_n[n]))
                out_gdl.write("\tnode: {\n")
                out_gdl.write("\t\t title: \"{0}\"\n".format(sort[n]+1))
                out_gdl.write("\t\t label: \"{0} - {1}\"\n".format(sort[n]+1, len(graph_n[n])))
                out_gdl.write("\t\t //{0}\n".format(graph_n[n]))
                out_gdl.write("\t}\n")
            i = 1
            for e in graph_e:
                out_txt.write("edge#{0} {1};{2}\n".format(i, sort[graph_e[e][0]]+1, sort[graph_e[e][1]]+1))
                out_gdl.write("\t edge: {\n")
                out_gdl.write("\t\t source: \"{0}\"\n".format(sort[graph_e[e][0]]+1))
                out_gdl.write("\t\t target: \"{0}\"\n".format(sort[graph_e[e][1]]+1))
                out_gdl.write("\t}\n")
                i = i + 1

            out_gdl.write("}\n")

if __name__ == '__main__':
	main()
