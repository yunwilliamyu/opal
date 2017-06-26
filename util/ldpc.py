# -*- coding: utf-8 -*-
"""
Generates LDPC hash function locations for k-mers
"""

import numpy as np
import argparse
import sys

def ldpc(k, t, _m):
    '''Generates a low density code matrix.

    k (int):        width of matrix / length of k-mer
    t (int):        number of 1's per row / row_weight
    _m (int):       suggestion for height of matrix / number of hashes
                    (may not be exactly what is returned as the actual
                    height of the matrix is computed below, but the result
                    will be at least _m in height
    
    returns 0/1 matrix
    '''
    if k % t != 0:
        raise ValueError('k should be multiple of t!')

    m = (int(np.ceil(_m*1.0/(k/t)) + 1)) * (k//t)
    w = m * t // k

    H_basic = np.zeros((m//w, k), dtype=np.bool)
    for i in range(m//w):
        for j in range(i * t, (i + 1) * t):
            H_basic[i, j] = 1

    H = H_basic.copy()
    for p in range(w - 1):
        perm_idx = np.random.permutation(k)
        H = np.vstack((H, H_basic[:, perm_idx]))
    return H

def hierarchical_ldpc(k, t1, t2, _m):
    '''Generates a low density code matrix using a hierarchical approach

    k (int):        width of matrix / length of k-mer
    t1 (int):       first level partitioning
    t2 (int):       number of 1's per row in the second level
    _m (int):       suggestion for height of matrix / number of hashes
                    (may not be exactly what is returned as the actual
                    height of the matrix is computed below, but the result
                    will be at least _m in height
    
    returns 0/1 matrix
    '''
    H = np.ndarray(0)
    H1 = ldpc(k, t1, _m)
    H2 = ldpc(t1, t2, _m)
    H_rows = []
    for row in range(_m):
        new_row_indices = H1[row].nonzero()[0][H2[row]]
        new_row = np.zeros(k, dtype=bool)
        new_row[new_row_indices] = True
        H_rows.append(new_row)
    H = np.vstack(H_rows)
    return H



def write_out(H, d, _m):
    '''Writes out LDPC matrix in text file for use with modified
    fasta2skm

    H (2D ndarray): 0/1 matrix, output of LDPC above
    d (string):     filename to write out too
    _m (int):       suggestion for height of matrix / number of hashes
                    (may not be exactly what is returned as the actual
                    height of the matrix is computed in LDPC above)

    '''
    k = len(H[0])
    t = H[0].sum()
    m = (int(np.ceil(_m*1.0/(k/t)) + 1)) * (k//t)
    w = m * t // k
    with open(d, 'w') as fout:    
        fout.write('%d %d\n'%(_m + 1, t))
        for j in range(t):
            fout.write('%d '%(j))
        fout.write('\n')
        st = m//w
        for i in range(_m):
            #sys.stdout.write('%d: '%(i))
            for j in range(k):
                if H[st + i, j] == 1:
                    fout.write('%d '%(j))
            fout.write('\n')

def ldpc_write(k, t, _m, d):
    '''Generates and writes out LDPC matrix'''
    H = ldpc(k, t, _m)
    write_out(H, d, _m)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', nargs=1)
    parser.add_argument('-t', nargs=1)
    parser.add_argument('-m', nargs=1)
    parser.add_argument('-d', nargs=1)
    args = parser.parse_args()
    #print args.k, args.t, args.m
    k = int(args.k[0])
    t = int(args.t[0])
    _m = int(args.m[0])
    d = args.d[0]
    ldpc_write(k, t, _m, d)

