
from nis import match
from tempfile import tempdir
from cv2 import ml_TrainData
import numpy as np
import itertools



'''
ON ALIGNMENT:

Global Seqnce alignment: 
Assumes two sequences are homologous (share common evolutionary ancestor) 
across their whole length. Tries to align the whole genome from tail to nose. 
Assumption for two genes in the same gene family 

Local Alignment: 
Assumes that  there are regions of homology within a bigger sequence, 
but also regions that do not share common ancestry. 
Reasonable assumption if checking how many genes two species share. 
If the species aren't extremely closely related, we expect shared genes and 
large portions of genome that's not shared. 

'''

'''
LOGIC: 

1) fill in matrix:
horizontal movements: gap penalty
diagonal movements: match / mismatch
store all 3 possible values (dynamic programming) and take the maximum as final value

2) Backtracing:
From bottom of matrix where highest score is, continue up to the upper left conner
If the bses match, go up diagonally, if not,  go to the highest value neighbor
'''



def score_matrix_v0(seq, ref, match_score, gap_cost):
    '''
    (v0): this funciton builds a score matrix using dynamic programming

    :param seq: user defined sequence
    :param ref: user defined reference sequence
    :param match_socre: user defined match score
    :param gap_cost: user defined gap penalty
    '''
    # initialize matrix
    matrix = np.zeros((5,6))
    matrix[0,1:] = [i for i in range(gap_cost,((len(ref)+1)*gap_cost),gap_cost)]
    matrix[1:,0] = [i for i in range(gap_cost,((len(seq)+1)*gap_cost),gap_cost)]

    for i in range(1, len(seq)+1): # +1 here bc matrix has one additional cell (0)
        possible_scores = 0
        for j in range(1, len(ref)+1):
            gap_j = matrix[i-1,j] + gap_cost # (i-1,j): left
            gap_i = matrix[i,j-1] + gap_cost # (i,j-1):up
            possible_scores = max(gap_j, gap_i)
            if seq[i-1] == ref[j-1]: # if match
                diag = matrix[i-1,j-1] + match_score
            else: # if mismatch
                diag = matrix[i-1,j-1] - match_score
            if diag > possible_scores:
                possible_scores = diag
            matrix[i,j] = possible_scores

    return matrix



def score_matrix_v1(seq, ref, match_score, gap_cost):
    '''
    (v1) this function builds a score matrix using itertools and numpy (more space efficient)
    :param seq:
    :param ref:
    :param match_socre: 
    :param gap_cost: 
    '''
    matrix = np.zeros((len(seq) + 1, len(ref) + 1))
    matrix[0,:] = np.linspace(0, (len(ref)*gap_cost),len(ref)+1) 
    matrix[:,0] = np.linspace(0, (len(seq)*gap_cost),len(seq)+1) 
    for i, j in itertools.product(range(1, matrix.shape[0]), range(1, matrix.shape[1])):
        match = matrix[i-1,j-1] + (match_score if seq[i-1] == ref[j-1] else - match_score)
        deletion = matrix[i-1,j] + gap_cost
        insertion = matrix[i, j-1] + gap_cost
        matrix[i,j] = max(match, deletion, insertion)
    return matrix




def needleman_wunsch(matrix):
    ''' 
    This function backtraces the score matrix and returns alignment result.
    LOGIC: 
    1) go to the last filled value of the matrix 
    2) breaks down socre matrix into 2x2 sub-matrices, find largest neighboring value or move diagonally
    3) append base to sequence to reconstruct alignment
    4) returns alignment (flips strings since backtracing)

    :param matrix: score matrix constructed 
    '''
    # for indexing: since seq and ref are one short than matrix dim
    (m,n) = np.subtract(matrix.shape, (1,1))

    # store alignment info
    align_seq = ''
    align_ref = ''

    # scores and directions
    diagonal = (1,1)
    left = (0,1)
    up = (1,0)


    while m>0 and n>0:
        # look at 2x2 matrix
        m_temp = matrix[m-1:m+1,n-1:n+1]
        # if sequences match, go back diagonally
        if seq[m-1] == ref[n-1]: 
            align_seq += seq[m-1]
            align_ref += ref[n-1] 
            (m, n) = np.subtract((m,n), diagonal)

        else: 
            # find max value neighbors (3 neighbor values)
            if np.argsort(m_temp.flatten()[:-1])[-1] == 1:
                align_ref += '-'
                align_seq += ref[n-1]
                (m,n) = np.subtract((m,n), up)

            elif np.argsort(m_temp.flatten()[:-1])[-1] == 2:
                align_seq += '-'
                align_ref += ref[n-1]
                (m,n) = np.subtract((m,n), left)

    return (align_ref[::-1], align_seq[::-1])
    



if __name__ == "__main__":
    match_score = 1
    gap_cost = -2
    seq = 'AGCT'
    ref = 'ATGCT'
    matrix = score_matrix_v1(seq, ref, match_score, gap_cost)
    align_result = needleman_wunsch(matrix)
    print("alignment result:")
    print("ref:", align_result[0])
    print("seq:", align_result[1])







