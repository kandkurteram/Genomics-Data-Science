import bisect 

class Index(object):
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k
        self.index = []
        
        for i in range(len(t) - k + 1):
            self.index.append((t[i:i+k], i))
        
        self.index.sort()

    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]
        i = bisect.bisect_left(self.index, (kmer, -1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i+=1
        return hits

def approximate_match(p, t, n):
    segment_length = int(round(len(p)/(n+1)))
    kmer_index = Index(t, 8)
    all_matches = set()
    all_hits = 0
    for i in range(n+1):
        start = i * segment_length
        end = min(len(p), (i+1) * segment_length)

        matches = kmer_index.query(p[start:end])
        all_hits += len(matches)

        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0,start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches), all_hits

class SubseqIndex(object):
    """ Holds a subsequence index for a text T """

    def __init__(self, genome, subsequence_length, ival):
        """ Create index from all subsequences consisting of subsequence_length characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.subsequence_length = subsequence_length
        self.index = []
        self.ival = ival
        self.span = 1 + ival * (subsequence_length-1)
        print('span: ', self.span)
        
        for i in range(len(genome) - self.span + 1):
            self.index.append((genome[i:i+self.span:ival], i))
        
        self.index.sort()

    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        sebseq = p[:self.span:self.ival]
        i = bisect.bisect_left(self.index, (sebseq, -1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != sebseq:
                break
            hits.append(self.index[i][1])
            i+=1
        return hits

def approximate_match_subseq(pattern, genome, allowed_mismatch_count):
    segment_length = int(round(len(pattern)/(allowed_mismatch_count+1)))
    ival = 3
    kmer_index = SubseqIndex(genome, segment_length, ival)
    print(kmer_index.index[:10])

    all_matches = set()
    all_hits = 0
    for i in range(allowed_mismatch_count+1):
        # start = i * segment_length
        start = i
        # end = min(len(pattern), (i+1) * segment_length)
        
        matches = kmer_index.query(pattern[start:])

        for m in matches:
            all_hits += 1   
            if m < start or m-start+len(pattern) > len(genome):
                continue
            mismatches = 0
            for j in range(0,len(pattern)):
                if not pattern[j] == genome[m-start+j]:
                    mismatches += 1
                    if mismatches > allowed_mismatch_count:
                        break 

            if mismatches <= allowed_mismatch_count:
                all_matches.add(m - start)
    return list(all_matches), all_hits

def ReadGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0]=='>':
                genome += line.rstrip()
    return genome

# genome = ReadGenome('/Users/rkandkurte/Documents/OneDrive/Learning/GenomicDataScience/Data/chr1_GRCh38_excerpt.fasta')
# pattern = 'GGCGCGGTGGCTCACGCCTGTAAT'
# matchOffsets, all_hits = approximate_match(pattern, genome, 2)
# matchOffsets.sort()
# print(matchOffsets)
# print(len(matchOffsets))
# print(all_hits)


genome = ReadGenome('/Users/rkandkurte/Documents/OneDrive/Learning/GenomicDataScience/Data/chr1_GRCh38_excerpt.fasta')
pattern = 'GGCGCGGTGGCTCACGCCTGTAAT'
matchOffsets, all_hits = approximate_match_subseq(pattern, genome, 2)
matchOffsets.sort()
print(matchOffsets)
print(len(matchOffsets))
print(all_hits)
