def Naive(p, t):
    import re
    matchOffsets= [m.start() for m in re.finditer('(?=' + p + ')', t)]
    return matchOffsets

def Naive_With_Counts(p, t):
    matchOffsets = []
    nCharacterComparisons = 0
    nAlignmentsTried = 0

    for i in range(0, len(t)-len(p)+1):
        nAlignmentsTried += 1
        match = True

        for j in range(len(p)):
            nCharacterComparisons += 1
            if p[j]!=t[i+j]:
                match=False
                break
            
        if match:
            matchOffsets.append(i)

    return matchOffsets, nAlignmentsTried, nCharacterComparisons


def Naive_2mm(p, t):
    matchOffsets = []
    # p_compliment = Complement(p)
    # print(p_compliment)
    for i in range(0, len(t)):
        matchBaseCount = 0
        for j in range(0, len(p)):
            # print(i, j)
            if len(t) - i > j:
                if p[j]==t[i+j]:# or p_compliment[j]==t[i+j]:
                    matchBaseCount +=1

        if len(p) - matchBaseCount <= 2:
            matchOffsets.append(i)
        
    return matchOffsets


def Naive_With_RC(p, t):
    import re
    matchOffsets= [m.start() for m in re.finditer('(?=' + p + ')', t)]
    matchOffsets += [m.start() for m in re.finditer('(?=' + ReverseComplement(p) + ')', t)]
    return list(set(matchOffsets))

def ReverseComplement(g):
    basePairs = {'A': 'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    reverseComplement = ''
    for base in g:
        reverseComplement = basePairs[base] + reverseComplement
    return reverseComplement

def Complement(g):
    basePairs = {'A': 'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    complement = ''
    for base in g:
        complement = complement + basePairs[base]
    return complement

def ReadGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0]=='>':
                genome += line.rstrip()
    return genome

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



# p = 'TA'
# t = 'TTAACTTCCAGGG'
# matchOffsets = Naive_With_RC(p, t)
# print(matchOffsets)

# g = 'ATGCTACGGGTTAAAAAAAACC'
# reverse = ReverseComplement(g)
# print(reverse)

# g = 'ATT'
# complement = Complement(g)
# print(complement)

# genome = ReadGenome('Data/lambda_virus.fa')
# pattern = 'AGGAGGTT'
# matchOffsets = Naive_2mm(pattern, genome)
# print(min(matchOffsets))

# sequences, qualities = ReadFastQ('Data/ERR037900_1.first1000.fastq')
# print(sequences[:1])
# print(qualities[:1])

# p = 'CTGT'
# t = 'AAAAAAAAAACTGTAAAAAAAAAACTTTAAAAAAAAAACGGGAAAAAAAAAA'
# matchOffsets = Naive_2mm(p, t)
# print(matchOffsets)

# p = 'word'
# t = 'there would have been a time for such a word'
# print(len(t))
# matchOffsets, nAlignmentsTried, nCharacterComparisons = Naive_With_Counts(p, t)
# print(matchOffsets, nAlignmentsTried, nCharacterComparisons)


# genome = ReadGenome('/Users/rkandkurte/Documents/OneDrive/Learning/GenomicDataScience/Data/chr1_GRCh38_excerpt.fasta')
# pattern = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
# matchOffsets, nAlignmentsTried, nCharacterComparisons = Naive_With_Counts(pattern, genome)
# print(matchOffsets, nAlignmentsTried, nCharacterComparisons)

# genome = ReadGenome('/Users/rkandkurte/Documents/OneDrive/Learning/GenomicDataScience/Data/chr1_GRCh38_excerpt.fasta')
# pattern = 'GGCGCGGTGGCTCACGCCTGTAAT'
# matchOffsets = Naive_2mm(pattern, genome)
# matchOffsets.sort()
# print(matchOffsets)