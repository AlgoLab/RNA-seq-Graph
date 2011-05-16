#!/usr/bin/env python3

#
# Gianluca Della Vedova (http://gianluca.dellavedova.org)
# Released under Affero GPL license, version 3
# Example
# cluster_fingerprints.py > reads.clean
from pygraph.classes.graph import graph
from pygraph.algorithms.accessibility import connected_components
from itertools import *
import pprint
from operator import itemgetter
import os

def generate_reads ():
    with open('raw-reads', 'r') as f:
        with open('reads', 'w') as r:
            with open('fingerprints', 'w') as fw:
                for line in f:
                    line = line.rstrip()
                    l=len(line)
                    for i in range(l-63):
                        read=line[i:i+64]
                        r.write(read + "\n")
                        for p in get_fingerprints(read):
                            fw.write(p + "\n")

def get_fingerprints (s):
    return (s[:32], s[-32:])


def argmax(pairs):
    return max(pairs, key=itemgetter(1))[0]

def argmin(sequence, fn=None):
    """Two usage patterns:
    argmin([s0, s1, ...], fn)
    argmin([(fn(s0), s0), (fn(s1, s1), ...])
    Both return the si with lowest fn(si)"""
    if fn is None:
        return min(sequence)[1]
    else:
        return min((fn(e), e) for e in sequence)[1]

def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def max_hamming_distance(s1, strings):
    return max([hamming_distance(s1, x) for x in strings])

def consensus_string(strings):
    l=len(strings[0])
    count=[dict() for x in range(l)]
    maxv=[0  for x in range(l)]
    consensus=['' for x in range(l)]

    for s in strings:
        for pos in range(l):
            if s[pos] not in count[pos]:
                count[pos][s[pos]]=0
            count[pos][s[pos]] +=1
            if count[pos][s[pos]] > maxv[pos]:
                maxv[pos]+=1
                consensus[pos]=s[pos]
    return(''.join(consensus))

def center_string(strings):
    return(argmin(zip([max_hamming_distance(x,strings) for x in strings], strings)))



def get_strings ():
    with open('fingerprints.uniq', 'r') as f:
        fingerprints=[x.rstrip() for x in f.readlines()]

    gr = graph()
    gr.add_nodes(fingerprints)

    fs = dict()
    components = dict()
    for f in fingerprints:
        for h in (f[:15], f[16:]):
            if h not in fs:
                fs[h]=[]
            fs[h].append(f)

#    print(fs)
    edges=set()
    for shared in fs.values():
        for f1 in shared:
            for f2 in shared:
                if f1 < f2 and hamming_distance(f1, f2) < 2:
#                    print(f1.rstrip() + ":" + f2.rstrip())
                    # components[f1].extend(components[f2])
                    # components[f2]=components[f1]
                    edges.add((f1, f2))

    [gr.add_edge(e) for e in edges]
    components=connected_components(gr)

    real_components=dict()
    for (k,v) in connected_components(gr).items():
        if v in real_components:
            real_components[v].append(k)
        else:
            real_components[v]=[k]

    all_consensus={k:consensus_string(v) for (k,v) in real_components.items()}
    all_consensus_distances={k:max_hamming_distance(v, real_components[k]) for
                             (k,v) in all_consensus.items()}

    pp = pprint.PrettyPrinter(indent=4)
    for k in real_components.keys():
        (chosen, cost) = (consensus_string(real_components[k]), all_consensus_distances[k])
        if all_consensus_distances[k] > 1:
            c=center_string(real_components[k])
            if max_hamming_distance(c, real_components[k]) < cost:
                (all_consensus_distances[k], all_consensus_distances[k]) = (c, max_hamming_distance(c, real_components[k]))

    # from real_components we compute map, which is a dictionary mapping a fingerprint
                # to its representative
    map={}
    for (k,v) in real_components.items():
        map.update({(w, all_consensus[k]) for w in v})

    with open('clusters', 'w') as f:
        pp = pprint.PrettyPrinter(indent=4,stream=f)
        pp.pprint(real_components)
    return(map)


if __name__ == '__main__':
#    os.system('grep -v "^>" ARHGAP4_1_mismatch.fa | grep -v N')
    generate_reads()
    os.system("sort -u fingerprints > fingerprints.uniq")
    map=get_strings()
    with open('reads', 'r') as f:
        for line in f:
            line = line.rstrip()
            fc=get_fingerprints(line)
            fc2=[map.get(f, f) for f in fc]
            print('> ')
            print(''.join(fc2))
