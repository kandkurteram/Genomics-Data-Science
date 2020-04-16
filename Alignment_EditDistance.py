
def editDistanceRecursive(x, y):
    ''' Recursive version (underperforming)
        Trying to find number of changes needed in "x" to convert into into "y".'''

    if len(x)==0:
        return len(y)
    elif len(y)==0:
        return len(x)
    else:
        delta = 0 if x[-1]==y[-1] else 1
        return min(
                    editDistanceRecursive(x,y[:-1]) + 1,
                    editDistanceRecursive(x[:-1],y) + 1,
                    editDistanceRecursive(x[:-1],y[:-1]) + delta,
                )

def editDistance(x, y):
    '''Dynamic Programming Version (optimized to give good performance).
       Trying to find number of changes needed in "x" to convert into into "y".'''
    matrix =[]
    matrix = [[0]*(len(y)+1) for _ in range(len(x)+1)]
    
    for i in range(len(x)+1):
        matrix[i][0] = i

    for i in range(len(y)+1):
        matrix[0][i] = i

    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHorizontal = matrix[i][j-1] + 1
            distVertical = matrix[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiagonal = matrix[i-1][j-1]
            else:
                distDiagonal = matrix[i-1][j-1] + 1

            matrix[i][j] = min(distHorizontal, distVertical, distDiagonal)
    
    return matrix[-1][-1]


def editDistanceApproximateMatch(x, y):
    '''To find minimum edit distance with approximate match. 
       i.e. We are trying to find the number of changes needed in "x" or "y" to find at least a substring match for "x" in "y"'''
    matrix =[]
    matrix = [[0]*(len(y)+1) for _ in range(len(x)+1)]
    
    for i in range(len(x)+1):
        matrix[i][0] = i

    for i in range(len(y)+1):
        matrix[0][i] = 0

    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHorizontal = matrix[i][j-1] + 1
            distVertical = matrix[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiagonal = matrix[i-1][j-1]
            else:
                distDiagonal = matrix[i-1][j-1] + 1

            matrix[i][j] = min(distHorizontal, distVertical, distDiagonal)
    
    # for row in matrix:
    #     for col in row:
    #         print(col, end=" ")
    #     print()
    return min(matrix[-1][:])


def GlobalAlignment(x, y):
    alphabet = ['A', 'C', 'G', 'T']
    score = [[0, 4, 2, 4, 8],
             [4, 0, 4, 2, 8],
             [2, 4, 0, 4, 8],
             [4, 2, 4, 0, 8],
             [8, 8, 8, 8, 8]]
    matrix = []
    matrix = [[0]*(len(y)+1) for _ in range(len(x)+1)]
    
    for i in range(1, len(x)+1):
        matrix[i][0] = matrix[i-1][0] + score[alphabet.index(x[i-1])][-1]

    for i in range(1, len(y)+1):
        matrix[0][i] = matrix[0][i-1] + score[-1][alphabet.index(y[i-1])]

    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHorizontal = matrix[i][j-1] + score[-1][alphabet.index(y[j-1])]
            distVertical = matrix[i-1][j] +  + score[alphabet.index(x[i-1])][-1]
            if x[i-1] == y[j-1]:
                distDiagonal = matrix[i-1][j-1]
            else:
                distDiagonal = matrix[i-1][j-1] +  + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]

            # print('i: ', i)
            # print('j: ', j)
            # print('x[i-1]: ', x[i-1])
            # print('y[j-1]: ', y[j-1])
            # print('distVertical: ', distVertical)
            # print('distHorizontal: ', distHorizontal)
            # print('distDiagonal: ', distDiagonal)
            matrix[i][j] = min(distHorizontal, distVertical, distDiagonal)
    
    # for row in matrix:
    #     for col in row:
    #         print(col, end=" ")
    #     print()

    return matrix[-1][-1]

# x = 'TATGCTA'
# y = 'TATGCTA'
# import time
# startTime = time.time()
# print(editDistanceRecursive(x, y))
# print("------%s seconds-------" % (time.time() - startTime))

# startTime = time.time()
# print(editDistanceOptimized(x,y))
# print("------%s seconds-------" % (time.time() - startTime))


from Alignment_NaiveExactMatch import ReadGenome
y = ReadGenome('/Users/rkandkurte/Documents/OneDrive/Learning/GenomicDataScience/Data/chr1_GRCh38_excerpt.fasta')
x = 'GATTTACCAGATTGAG'
print(editDistanceApproximateMatch(x,y))