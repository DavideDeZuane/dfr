from utils import sample_supp, sample_H, compute_syndrome, sample_near_codeword, hamming_wt
from decoders import create_lookup, bf_out_of_place
import time

import math

p = 1021 #dimensione del blocco della matrice
v = 13
t = 30 


max_failures = 10; max_tx = 10000; num_backup = 100

# thresholds (abbiamo una soglia per ogni iterazione) questo fanno da regola decisionale, 11 iterazioni con soglia decrescente
# per ogni bit della parola ricevuta si vede quante equazioni di parità viola, se questo numero è più grande della threshold si flippa il bit
# ci dicono quanto deve essere colpevole un bit per essere flippato 
thresholds = [12, 11, 10, 10, 9, 9, 8, 8, 7, 7, 7]
h1_supp, h2_supp, row_supps = sample_H(p, v)


in_time = time.time()    
L = create_lookup(h1_supp, h2_supp, p, v)
out_time = time.time()

print("Table computed, time = ",out_time - in_time)

e_supp = sample_supp(2*p, t)
s = compute_syndrome(h1_supp, h2_supp, e_supp, p)

in_time = time.time()
dec_e, dec_e_2 = bf_out_of_place(h1_supp, h2_supp, s, L, p, thresholds)
out_time = time.time()

print("Decoding end, time = ",out_time - in_time)


floor_std = 0
floor_look_up = 0

for u in range(math.floor(v/2), 3, -1):
    
    pr_u = 2*p*math.comb(v,u)*math.comb(2*p-v, t-u)/math.comb(2*p,t)
    
    #start simulation
    dfr_std = 0; dfr_look_up = 0; num_tx = 0

    while (min(dfr_std, dfr_look_up)< max_failures) & (num_tx < max_tx):

        num_tx += 1
                
        #sample error    
        e_supp = sample_near_codeword(h1_supp, h2_supp, p, v, t, u)
   #     print(e_supp)
      #  e_supp = sample_supp(2*p, t)

        #compute syndrome
        s = compute_syndrome(h1_supp, h2_supp, e_supp, p)

        dec_e, dec_e_2 = bf_out_of_place(h1_supp, h2_supp, s, L, p, thresholds)

 #       print(e_supp)
        #check if decoding is successful for standard decoder
        for i in e_supp:
            dec_e[i] = (dec_e[i]+1)%2

        if hamming_wt(dec_e) > 0:
            dfr_std += 1

        #check if decoding is successful for improved decoder
        for i in e_supp:
            dec_e_2[i] = (dec_e_2[i]+1)%2
            
        if hamming_wt(dec_e_2) > 0:
      #      print("------ ERROR")
            dfr_look_up += 1
            
        #print results
        if (num_tx%num_backup) == 0:
            print("--> Num tx = ",num_tx,", u = ",u,", DFR(STD) = ",dfr_std/num_tx,", DFR(LU) = ", dfr_look_up/num_tx)
#            print("--> Floor contribution, STD : ",math.log(pr_u) + math.log(dfr_std/num_tx,2),", LU : ", math.log(pr_u) + math.log(dfr_look_up/num_tx,2))                 
            print("----> FLOOR(std) = ",floor_std + pr_u * dfr_std/num_tx,", FLOOR(look_up) = ",floor_look_up + pr_u * dfr_look_up/num_tx)

    
    floor_std += (pr_u * dfr_std/num_tx)
    floor_look_up += (pr_u * dfr_look_up/num_tx)
        
    print("u = ",u,", Num tx = ",num_tx,", DFR(STD) = ",dfr_std/num_tx,", DFR(LU) = ", dfr_look_up/num_tx)
    print("FLOOR(STD) >= ",floor_std,", FLOOR(LU) >= ",floor_look_up)
    print("-------------------------------------------")
