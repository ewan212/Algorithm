
# Global alignment

from nis import match
from turtle import up
import numpy as np
import itertools

for i in range(-2,((len(ref)+1)*gap_cost),-2):
    print(i)



def score_matrix_v0(seq, ref, match_score, gap_cost):
    # initialize matrix
    matrix = np.zeros((5,6))
    matrix[0,1:] = [i for i in range(-2,((len(ref)+1)*gap_cost),-2)]
    matrix[1:,0] = [i for i in range(-2,((len(seq)+1)*gap_cost),-2)]

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



# matrix = np.zeros((len(b) + 1, len(a) + 1))
# matrix[0,1:] = [i for i in range(0, (len(a)*gap_cost),-2)]
# matrix[0,1:] = [i for i in range(0,(len(b)*gap_cost),-2)]

def score_matrix_v1(seq, ref, match_score, gap_cost):
    '''more concise using itertools'''
    matrix = np.zeros((len(seq) + 1, len(ref) + 1))
    matrix[0,1:] = [i for i in range(-2,((len(ref)+1)*gap_cost),-2)]
    matrix[1:,0] = [i for i in range(-2,((len(seq)+1)*gap_cost),-2)]        
    for i, j in itertools.product(range(1, matrix.shape[0]), range(1, matrix.shape[1])):
        match = matrix[i-1,j-1] + (match_score if seq[i-1] == ref[j-1] else - match_score)
        deletion = matrix[i-1,j] + gap_cost
        insertion = matrix[i, j-1] + gap_cost
        matrix[i,j] = max(match, deletion, insertion)
    return matrix


match_score = 1
gap_cost = -2
seq = 'AGCT'
ref = 'ATGCT'
matrix = score_matrix_v1(seq, ref, match_score, gap_cost)





# find max value in matrix
# check where current value is calculated from
# append letter to seq and ref
# subtract index



(m,n) = np.subtract(matrix.shape, (1,1))

align_seq = ''
align_ref = ''


diagonal = (1,1)
left = (0,1)
up = (1,0)

while m>0 and n>0:
    m_temp = matrix[m-1:m+1,n-1:n+1]
    if seq[m-1] == ref[n-1]: 
        align_seq += seq[m-1]
        align_ref += ref[n-1] 
        (m, n) = np.subtract((m,n), diagonal)


    else: 
        # max value neighbor 
        if np.argsort(m_temp.flatten()[:-1])[-1] == 1:
            align_ref += '-'
            align_seq += ref[n-1]
            (m,n) = np.subtract((m,n), up)

        elif np.argsort(m_temp.flatten()[:-1])[-1] == 2:
            
            align_seq += '_'
            align_ref += ref[n-1]
            (m,n) = np.subtract((m,n), left)

print(align_ref[::-1])
print(align_seq[::-1])




    



        

