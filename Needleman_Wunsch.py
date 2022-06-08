
# Global alignment

from nis import match
from tempfile import tempdir
from turtle import setundobuffer, up
from cv2 import ml_TrainData
import numpy as np
import itertools


def score_matrix_v0(seq, ref, match_score, gap_cost):
    '''first version'''
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




def score_matrix_v1(seq, ref, match_score, gap_cost):
    '''using itertools and numpy'''
    matrix = np.zeros((len(seq) + 1, len(ref) + 1))
    matrix[0,:] = np.linspace(0, (len(ref)*gap_cost),len(ref)+1) 
    matrix[:,0] = np.linspace(0, (len(seq)*gap_cost),len(seq)+1) 
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

h = matrix
h_flip = np.flip(np.flip(matrix,0), 1)
i_, j_ = np.unravel_index(h_flip.argmax(), h_flip.shape)
i, j = np.subtract(h.shape, (i_ + 1, j_ + 1))




# find max value in matrix
# check where current value is calculated from
# append letter to seq and ref
# subtract index

def needleman_wunsch(matrix):
    # for indexing matrix
    (m,n) = np.subtract(matrix.shape, (1,1))

    # store alignment info
    align_seq = ''
    align_ref = ''

    # scores and directions
    diagonal = (1,1)
    left = (0,1)
    up = (1,0)


    while m>0 and n>0:
        m_temp = matrix[m-1:m+1,n-1:n+1]
        # if sequences match, go back diagonally
        if seq[m-1] == ref[n-1]: 
            align_seq += seq[m-1]
            align_ref += ref[n-1] 
            (m, n) = np.subtract((m,n), diagonal)

        else: 
            # find max value neighbor 
            if np.argsort(m_temp.flatten()[:-1])[-1] == 1:
                align_ref += '-'
                align_seq += ref[n-1]
                (m,n) = np.subtract((m,n), up)

            elif np.argsort(m_temp.flatten()[:-1])[-1] == 2:
                align_seq += '_'
                align_ref += ref[n-1]
                (m,n) = np.subtract((m,n), left)

    return align_ref[::-1], align_seq[::-1]
    

# base case: m,n==0
max(m_temp.flatten()[:-1])
back_max = np.where(m_temp == max(m_temp.flatten()[:-1]))
(a,b) = (back_max[0][0], back_max[1][0])

def traceback(m,n):
    if m==0 or n==0:
        return align_seq, align_ref 
    else: 
        if seq[m-1] == ref[n-1]: 
            align_seq += seq[m-1]
            align_ref += ref[n-1] 
            (m, n) = np.subtract((m,n), diagonal)

        else: 
            m_temp = matrix[m-1:m+1,n-1:n+1]
            # find max value neighbor 
            back_max = np.where(m_temp == max(m_temp.flatten()[:-1]))
            (a,b) = (back_max[0][0], back_max[1][0])

    return traceback(m,n)








