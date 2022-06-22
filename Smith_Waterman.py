
'''
Score matrix: exceptg all negativ values become 0
(similar problem: finidng the longest substring)

'''
# from nis import match
from tempfile import tempdir
# from cv2 import ml_TrainData
import numpy as np
import itertools

def score_matrix(seq, ref, match_score, gap_cost):
    '''
    This function builds a score matrix
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
        matrix[i,j] = max(match, deletion, insertion, 0)
    return matrix



seq = 'AGCT'
ref = 'ATGCT'
match_score = 1
gap_cost = -2
matrix = score_matrix(seq, ref, match_score, gap_cost)



i,j  = np.unravel_index(matrix.argmax(), matrix.shape) # max value index


def trace_back(matrix):
    seq_align = ''
    fin = []
    i,j  = np.unravel_index(matrix.argmax(), matrix.shape)
    if matrix[i-1,j-1] == 0:
        return seq_align, fin.append((i,j))
     
    seq_align += seq[i-1]
    matrix = matrix[:i, :j]
    i, j = i-1, j-1
    return trace_back(matrix)





def traceback(matrix, b, b_='', old_i=0):
    '''
    this function uses recursion to trace back the scoring matrix
    '''
    # flip H to get index of **last** occurrence of H.max() with np.argmax()
    H_flip = np.flip(np.flip(matrix, 0), 1)
    i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape) #0,0
    i, j = np.subtract(matrix.shape, (i_ + 1, j_ + 1))  # (i, j) are **last** indexes of H.max()
    if matrix[i, j] == 0:
        return b_, j # seq string and index
    b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_
    return traceback(matrix[0:i, 0:j], b, b_, i)





