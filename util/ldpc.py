# -*- coding: utf-8 -*-
"""
Generates LDPC hash function locations for k-mers
"""

import numpy as np
import argparse
import sys

def ldpc(k, t, _m, d):
    if k % t != 0:
        raise ValueError('k should be multiple of t!')

    m = (int(np.ceil(_m*1.0/(k/t)) + 1)) * (k/t)
    w = m * t / k

    H_basic = np.zeros((m/w, k), dtype=np.bool)
    for i in range(m/w):
        for j in range(i * t, (i + 1) * t):
            H_basic[i, j] = 1

    H = H_basic.copy()
    for p in range(w - 1):
        perm_idx = np.random.permutation(k)
        H = np.vstack((H, H_basic[:, perm_idx]))

    with open(d, 'w') as fout:    
        fout.write('%d %d\n'%(_m + 1, t))
        for j in range(t):
            fout.write('%d '%(j))
        fout.write('\n')
        st = m/w
        for i in range(_m):
            #sys.stdout.write('%d: '%(i))
            for j in range(k):
                if H[st + i, j] == 1:
                    fout.write('%d '%(j))
            fout.write('\n')


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
    ldpc(k, t, _m, d)

