
# Global alignment

from turtle import up
import numpy as np
import itertools

match_score = 1
gap_cost = -2
seq = 'AGCT'
ref = 'ATGCT'


def score_matrix_v0(seq, ref, match_score, gap_cost):
    # initialize matrix
    matrix = np.zeros((5,6))
    matrix[0,1:] = [i for i in range(-2,(len(ref)*gap_cost),-2)]
    matrix[1:,0] = [i for i in range(-2,(len(seq)*gap_cost),-2)]

    for i in range(1, len(seq)+1):
        possible_scores = 0
        for j in range(1, len(ref)+1):
            gap_j = matrix[i-1,j] + gap_cost
            gap_i = matrix[i,j-1] + gap_cost
            possible_scores = max(gap_j, gap_i)
            if seq[i-1] == ref[j-1]:
                diag = matrix[i-1,j-1] + match_score
            else:
                diag = matrix[i-1,j-1] - match_score
            if diag > possible_scores:
                possible_scores = diag
            matrix[i,j] = possible_scores

    return matrix





def score_matrix_v1(seq, ref, match_score, gap_cost):
    '''more concise using itertools'''
    matrix = np.zeros((len(seq) + 1, len(ref) + 1))
    matrix[0,1:] = [i for i in range(-2,(len(ref)*gap_cost),-2)]
    matrix[1:,0] = [i for i in range(-2,(len(seq)*gap_cost),-2)]        
    for i, j in itertools.product(range(1, matrix.shape[0]), range(1, matrix.shape[1])):
        match = matrix[i-1,j-1] + (match_score if seq[i-1] == ref[j-1] else - match_score)
        deletion = matrix[i-1,j] + gap_cost
        insertion = matrix[i, j-1] + gap_cost
        matrix[i,j] = max(match, deletion, insertion)
    return matrix

matrix = score_matrix_v1(seq, ref, match_score, gap_cost)



def func():
    '''returns starting index of backtraced value'''

np.unravel_index(29, (5,6))
matrix[4,5]

m,n = np.unravel_index(np.argmax(matrix), (5,6))
if seq[m-1] == ref[n-1]:
    new_index = (m-1,n-1)
    
else: 
    # find max
    max(matrix[m-1, n-1], matrix[m-1,n], matrix[m,n-1])
    np.argmax(matrix[m-1:m+1,n-1:n+1])
m=4
n=5
np.unravel_index(np.argmax(matrix[m-1:m+1,n-1:n+1]), (2,2))


np.flip(matrix,0)
m_flip = np.flip(np.flip(matrix, 0), 1)
# index of max value in matrix
i, j = np.unravel_index(m_flip.argmax(), m_flip.shape)


# find max value 
# check string match
# 
np.subtract(matrix.shape, (1,1))


