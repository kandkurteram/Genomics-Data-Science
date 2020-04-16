from bm_preproc import BoyerMoore

def boyer_moore(p, t, p_bm):
    Occurences = []
    i=0
    while i < len(t) - len(p) + 1:
        shift = 1
        pattern_match=1
        for j in range(len(p)-1, -1, -1):
            if(p[j]!=t[i+j]):
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                pattern_match = 0;
                break;
        
        if pattern_match == 1:
            Occurences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        
        i += shift
            
    return Occurences

def boyer_moore_with_counts(p, t, p_bm):
    Occurences = []
    i=0
    nCharacterComparisons = 0
    nAlignmentsTried = 0
    while i < len(t) - len(p) + 1:
        shift = 1
        pattern_match=1
        nAlignmentsTried += 1
        for j in range(len(p)-1, -1, -1):
            nCharacterComparisons += 1

            if(p[j]!=t[i+j]):
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                pattern_match = 0;
                break;
        
        if pattern_match == 1:
            Occurences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        
        i += shift
            
    return Occurences, nAlignmentsTried, nCharacterComparisons

def ReadGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0]=='>':
                genome += line.rstrip()
    return genome

def approximate_match(p, t, n):
    segment_length = int(round(len(p)/(n+1)))
    all_matches = set()
    all_hits = 0
    for i in range(n+1):
        start = i * segment_length
        end = min(len(p), (i+1) * segment_length)
        p_bm = BoyerMoore(p[start:end], alphabet = 'ACGT')
        matches = boyer_moore(p[start:end], t, p_bm)
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

# t = 'CATGCGTT'
# p = 'TGCG'
# p_bm = BoyerMoore(p)

# Occurences = boyer_moore(p, t, p_bm)
# print(Occurences)

# genome = ReadGenome('/Users/rkandkurte/Documents/OneDrive/Learning/GenomicDataScience/Data/chr1_GRCh38_excerpt.fasta')
# pattern = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
# p_bm = BoyerMoore(pattern)
# Occurences, nAlignmentsTried, nCharacterComparisons = boyer_moore_with_counts(pattern, genome, p_bm)
# print(Occurences, nAlignmentsTried, nCharacterComparisons)

# t = 'CACTTAATTTG'
# p = 'AACTTG'
# print(approximate_match(p, t, 2))

genome = ReadGenome('/Users/rkandkurte/Documents/OneDrive/Learning/GenomicDataScience/Data/chr1_GRCh38_excerpt.fasta')
pattern = 'GGCGCGGTGGCTCACGCCTGTAAT'
matchOffsets, all_hits = approximate_match(pattern, genome, 2)
matchOffsets.sort()
print(matchOffsets)
print(len(matchOffsets))
print(all_hits)