from utils import compute_syndrome, hamming_wt
import numpy as np

def create_lookup(h1_supp, h2_supp, p, v):
    '''
    Creates a look-up table with pairs of values (s, i), where s is the sydnrome corresponding to trapping sey i
    '''
    
    L = {} # Dictionary

    #get entry for table, coming from fist circulant
    for i in range(p):
        
        e_supp = [(i+j)%p for j in h1_supp]
        
        #get lookup table entry
        s = compute_syndrome(h1_supp, h2_supp, e_supp, p)
        supp_s = []
        for j in range(p):
            if s[j]:
                supp_s.append(j)
                
        L[''.join(str(supp_s))] = e_supp

    # print(f"supp_e = {supp_e}")

    
    for i in range(p):
        
        e_supp = [p+(i+j)%p for j in h2_supp]
                
        #get lookup table entry
        s = compute_syndrome(h1_supp, h2_supp, e_supp, p)
        supp_s = []
        for j in range(p):
            if s[j]:
                supp_s.append(j)
                
        L[''.join(str(supp_s))] = e_supp

    return L



def bf_out_of_place(h1_supp, h2_supp, s, L, p, thresholds):
    '''
    Out of place bf
    '''
    
    num_iter = len(thresholds)
    iter = 0
    flag_lookup = 0
    
    # Error estimate: initialize the error vector as a list of zeros
    dec_e  = np.zeros(2*p, dtype = np.uint8)  

    while (iter < num_iter) & (hamming_wt(s) != 0):
        
        delta_s = np.zeros(p, dtype = np.uint8)
        
        #flips for first circulant
        for i in range(p):
            # qui per la posizione i si calcola quante componenti della sindrome sono a 1, tra quelle influenzate da r_i
            sigma_i = sum([s[(j+i)%p] for j in h1_supp])
            # se supera la soglia allora si flippa il bit
            if sigma_i >= thresholds[iter]:
                
                #flip error position
                dec_e[i] = (1+dec_e[i])%2

                #flip syndrome
                for j in h1_supp:
                    pos = (i+j)%p
                    delta_s[pos] = (delta_s[pos]+1)%2
        
        #flips for second circulant
        for i in range(p):
            sigma_i = sum([s[(j+i)%p] for j in h2_supp])
            
            if sigma_i >= thresholds[iter]:
                
                #flip error position
                dec_e[p+i] = (1+dec_e[p+i])%2

                #flip syndrome
                for j in h2_supp:
                    pos = (i+j)%p
                    delta_s[pos] = (delta_s[pos]+1)%2                    
                    
        
        #flip syndrome
        s = (s+delta_s)%2
        
        #see if syndrome is in look-up table        
        if flag_lookup == 0:
            supp_s = []
            for i in range(p):
                if s[i]:
                    supp_s.append(i)
            target_key = ''.join(str(supp_s))
            
            e_supp = L.get(target_key)

            #if a position is found, use look-up
            if e_supp is not None:
                
               # print("######## ooooooook, pos = ",e_supp)
                
                flag_lookup = 1
                
                #update error vector
                dec_e_2 = dec_e.copy()
                
                for i in e_supp:
                    dec_e_2[i] = (dec_e_2[i]+1)%2

        iter += 1

    # At the end of the loop, if any decoding is found.
    if flag_lookup == 0:
        dec_e_2 = dec_e.copy()

    # Return the two decoded vectors        
    return dec_e, dec_e_2


def my_BF_max_improved(supp_col_tot, supp_row_tot, s, n, num_iterations, L):
    '''
    BF-Max algorithm, with the improvement of the look-up table.

    Parameters:
        supp_col_tot (list): the list of the whole columns supports,
        supp_row_tot (list): the list of the whole rows supports,
        s (list): The binary syndrome vector,
        n (int): The number of bits (length of the syndrome vector),
        t (int): The number of iterations to perform for the bit-flipping algorithm
        L (list): list of trapping sets.
    
    Returns:
        dec_e (np.array): The support of the estimated error vector.
    '''
    
    # Compute counters (the count of '1's in the syndrome for each column index)
    counters = [compute_counter(s, supp_col_tot[i]) for i in range(n)]
    
    # Error estimate: initialize the error vector as a list of zeros
    dec_e  = np.zeros(n, dtype = np.uint8)  

    # Flag for look-up table usage
    flag_lookup = 0
    Iter = 0
        
    while (Iter < num_iterations) & (hamming_weight(s) != 0):
        
        # Find the index to flip based on the counter values
        i = index_to_flip(counters, n)
        
        # Flip the bit in the error vector (toggle between 0 and 1)
        dec_e[i] = (dec_e[i] + 1) % 2
        
        # Update syndrome and counters
        for j in supp_col_tot[i]:
            if s[j] == 1:
                incr = -1
            else:
                incr = 1
                
            # Update counters
            for ell in supp_row_tot[j]:
                counters[ell] += incr
            
            # Update the syndrome vector (flip the bit at position j)
            s[j] = (s[j] + 1) %2
            

    # Search in the look-up table, given the syndrome. 
        if flag_lookup == 0:
            target_s = str(compute_support(s)) 
            delta_e = L.get(target_s)

            if delta_e is not None:
                print("\n   ok, sono entrato")
                dec_e_2 = dec_e.copy() 
                for i in delta_e:
                    dec_e_2[i] = (dec_e_2[i] + 1) %2
                    flag_lookup = 1

    # end while

    # At the end of the loop, if any decoding is found.
    if flag_lookup == 0:
        dec_e_2 = dec_e.copy()

    # Return the two decod vectors        
    return dec_e, dec_e_2