def Overlap(a, b, min_overlap_length=3):
    ''' Finds the length of overlap greater than "min_overlap_length" between "a" and "b".'''
    start = 0
    while True:
        start = a.find(b[:min_overlap_length], start)

        #if no overlap of "b" prefix found in suffix of "a" then return 0
        if(start == -1):
            return 0
        
        #if overlap of "b" prefix is found in suffix of "a" then check if whole "b" matches with suffix in "a" starting at found overlap. 
        #if match found return the overlap length of "a"
        if(b.startswith(a[start:])):
            return len(a)-start
        
        #increment start
        start+=1

def Naive_Overlap_Map(reads, min_overlap_length):
    '''Finds the overlap greater than "min_overlap_length" between provided reads "reads".'''
    from itertools import permutations

    overlaps = {}
    #Find all permutations of 2 reads from the provided list of reads and loop through each permutation
    for a,b in permutations(reads, 2):
        # print(a,b)
        #Find if there is overlap of minimum length "min_overlap_length" in current pair of reads
        overlaplenght = Overlap(a, b, min_overlap_length)

        #If overlap is found, record it in the "overlaps" dictionary
        if(overlaplenght>0):
            overlaps[(a, b)] = overlaplenght
        
    return overlaps

from collections import defaultdict
def overlap_graph(reads, k):
    # Make index
    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i+k]].add(read)

    # print(index)
    # Make graph
    graph = defaultdict(set)
    for r in reads:
        for o in index[r[-k:]]:
            if r != o:
                if Overlap(r, o, k):
                    graph[r].add(o)

    # print(graph)
    edges = 0
    for read in graph:
        edges += len(graph[read])
    # print(edges)
    return(edges, len(graph))

import itertools
def ShortestCommonSuperstring(ss):
    shortest_sup = None
    all_scs = []
    for ssperm in itertools.permutations(ss):
        # print(ssperm)
        sup = ssperm[0]
        for i in range(len(ss)-1):
            olen = Overlap(ssperm[i], ssperm[i+1], min_overlap_length = 1)
            sup += ssperm[i+1][olen:]

        if shortest_sup is None or len(sup) <= len(shortest_sup):
            shortest_sup = sup
            all_scs.append(sup)
    return shortest_sup, all_scs

def PickMaximalOverlap(reads, k):
    reada, readb = None, None
    best_olen = 0
    for a, b in itertools.permutations(reads, 2):
        olen = Overlap(a, b, k)
        if olen > best_olen:
            reada = a
            readb = b
            best_olen = olen
    return reada, readb, best_olen

def greedy_ShortestCommonSuperstring(reads, k):
    reada, readb, olen = PickMaximalOverlap(reads, k)
    while olen>0:
        reads.remove(reada)
        reads.remove(readb)
        reads.append(reada + readb[olen:])
        reada, readb, olen = PickMaximalOverlap(reads, k)
    return ''.join(reads)

def de_bruijn_ise(st, k):
    edges = []
    node = set()
    for i in range(len(st) - k + 1):
        edges.append((st[i:i+k-1], st[i+1:i+k]))
        node.add(st[i:i+k-1])
        node.add(st[i+1:i+k])
    return node, edges

def ReadFastQ(filename):
    sequences = []
    qualities = []

    with open(filename, 'r') as f:
        while True:
            f.readline()
            seq = f.readline().rstrip()
            f.readline()
            qual = f.readline().rstrip()

            if len(seq) == 0:
                break;

            sequences.append(seq)
            qualities.append(qual)

    return sequences, qualities



# a = 'TGATCC'
# b = 'TCCGATTCC'
# print(Overlap(a, b))

# reads = ['TACTTGCA', 'GCATTGA', 'GGCTCA', 'TCGAAT']
# print(Naive_Overlap_Map(reads, 3))

# from Alignment_NaiveExactMatch import ReadFastQ
# reads, _ = ReadFastQ('/Users/rkandkurte/Documents/OneDrive/Learning/GenomicDataScience/Data/ERR266411_1.for_asm.fastq')
# edges, graphlen = overlap_graph(reads, 30)
# print(edges)
# print(graphlen)

# shortest_sup, all_scs = ShortestCommonSuperstring(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])
# print(shortest_sup)
# print(all_scs)

# node, edges = de_bruijn_ise('ACGCGTCG', 3)
# print(node)
# print(edges)

seqs, quals = ReadFastQ('/Users/rkandkurte/Documents/OneDrive/Learning/GenomicDataScience/Data/ads1_week4_reads.fq')
# print(seqs, quals)
from collections import defaultdict
for k in range(100, 1, -1):
        genome = greedy_ShortestCommonSuperstring(seqs, k)
        if len(genome) == 15894:
            print(genome.count('A'))
            print(genome.count('T'))
            # print(genome)
            break