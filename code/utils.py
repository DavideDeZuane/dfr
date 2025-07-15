import random
from math import floor

import numpy as np

def sample_supp(m, w):
    '''
    Sample support of vector with length m and weight w
    '''
    
    supp = []
    while len(supp) < w:
        
        #new position
        pos = random.randrange(0,m)

        #see if position was already in supp
        duplicated = 0
        i = 0
        while (duplicated == 0) & (i<len(supp)):
            if supp[i] == pos:
                duplicated = 1
            else:
                i += 1
                
        if duplicated == 0:
            supp.append(pos)
    
    #sorting
    supp.sort()
    
    return supp

def sample_near_codeword(h1_supp, h2_supp, p, v, t, u):
    '''
    sample error vector that has intersection of size u with a column of H
    '''    
    e_supp = []
    
    #select column
    pos = random.randrange(0,2*p)
    
    #set target col
    if pos < p:
        target_col = [(i+pos)%p for i in h1_supp]
    else:
        target_col = [p+(i+pos)%p for i in h2_supp]
                
    #generate random permutation 
    perm = np.array(range(v))
    random.shuffle(perm)

    for i in range(u):    
        e_supp.append(target_col[i])

    #sample remaining positions
    while len(e_supp)<t:
        
        #new position
        pos = random.randrange(0,2*p)

        if pos not in target_col:
            
            #see if position was already in supp
            duplicated = 0
            i = 0
            while (duplicated == 0) & (i<len(e_supp)):
                if e_supp[i] == pos:
                    duplicated = 1
                else:
                    i += 1
                    
            if duplicated == 0:
                e_supp.append(pos)
    
    #sorting
    e_supp.sort()
    
    return e_supp


def sample_H(p, v):
    '''
    Sample double circulant matrix with circulant size p, circulant weight v.
    h1 and h2 are the support of the first columns of each circulant
    '''
    h1_supp = sample_supp(p, v)
    h2_supp = sample_supp(p,v)
    
    #create first row for first circulant
    first_row_supp = []
    for i in h1_supp:
        if i == 0:
            first_row_supp.append(i)
        else:
            first_row_supp.append(p-i)
    first_row_supp.sort()

    
    #compute first row for second circulant    
    for i in h2_supp:
        if i == 0:
            first_row_supp.append(p+i)
        else:
            first_row_supp.append(2*p-i)
    first_row_supp.sort()
    
    
    #create row supports
    row_supps = []
    for i in range(p):
        this_row_supp = [p*floor(j/p) + (i+j)%p for j in first_row_supp]
        this_row_supp.sort()
        row_supps.append(this_row_supp)

    return h1_supp, h2_supp, row_supps



def compute_syndrome(h1_supp, h2_supp, e_supp, p):
    '''
    Compute syndrome of vector with support e_supp
    '''
    
    s = np.zeros(p, dtype = np.uint8)

    for i in e_supp:
        if i < p:
            for j in h1_supp:
                s[(j + i) % p] ^= 1
        else:
            for j in h2_supp:
                s[(j + i) % p] ^= 1
    return s

def hamming_wt(s):
    '''
    return hamming of input vector
    '''

    return np.sum(s)