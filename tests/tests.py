import numpy as np
import pandas as pd
import spkmenas
import create_random_tests
import os
import time
import subprocess
import random

def load():

    #Load Everything in C
    W_c = pd.read_csv(r"C:\Users\Tomer\CLionProjects\spectral_k_means\tests\out_w.txt", header=None).to_numpy()
    D_c = pd.read_csv(r"C:\Users\Tomer\CLionProjects\spectral_k_means\tests\out_D.txt", header=None).to_numpy()
    D_half_c = pd.read_csv(r"C:\Users\Tomer\CLionProjects\spectral_k_means\tests\out_D_half.txt", header=None).to_numpy()
    lap_c = pd.read_csv(r"C:\Users\Tomer\CLionProjects\spectral_k_means\tests\out_lap.txt", header=None).to_numpy()

    #Load data to Python
    data = pd.read_csv(r"C:\Users\Tomer\CLionProjects\spectral_k_means\tests\input.txt",header=None).to_numpy()

    W_py = spkmenas.weighted_adj_mat(data)
    D_py = spkmenas.diagonal_degree_mat(W_py)
    D_half_py = D_py.copy()
    spkmenas.pow_diag(D_half_py,-0.5)
    lap_py = spkmenas.l_norm_mat(W_py,D_half_py)

    return (W_c, W_py,D_c,D_py,D_half_c,D_half_py,lap_c,lap_py)





def is_same(A, B):
    n1, m1 =  A.shape
    n2, m2 = B.shape

    if (n1!=n2 or m1!= m2):
        print("Dimentions aren't equal problem")
        print (f"C dimentins are {A.shape}")
        print (f"Python dimentins are {B.shape}")
        return False
    bol = True
    max_error =0
    c=0
    for i in range(n1):
        for j in range (m2):
            if (abs(A[i, j]- B[i, j])>10**-5): #Error Should be very small
                #print(f"The [{i},{j}] element isn't the same. ")
                #print(f"in C: {A[i, j]} , in Python {B[i, j]} diffrent is : {abs((A[i, j])- (B[i, j]))}")
                bol = False
                val = abs (A[i,j]-B[i,j])
                if (val>max_error):
                    max_error = max(max_error,val)
                    max_i =i
                    max_j = j
                    c+=1

    if (bol):
        print ("Passed Successfully")
    else: #Repory after an Error
        print ("*****************************************")
        print (f"Max error is {max_error}")
        print (f"largest error is in indices {(max_i,max_j)}")
        print (f"total errors number is {c} out of {A.shape[0]*A.shape[1]}")
        print ("*****************************************")
    return bol

def test():
    W_c, W_py,D_c,D_py,D_half_c,D_half_py,lap_c,lap_py = load()
    print ("*"*20+ " W testing "+"*"*20)
    if not is_same(W_c,W_py):
        return False
    print ("*"*20+ " D testing "+"*"*20)
    if not is_same(D_c,D_py):
        return False
    print ("*"*20+ " D_half testing "+"*"*20)
    if not is_same(D_half_c,D_half_py):
        return False
    print ("*"*20+ " Laplacian testing "+"*"*20)
    if not is_same(lap_c,lap_py):
        return False
    return True

def create_random_test(n,m):

    mat = np.random.uniform(-10,10,(n,m))
    mat = (np.round(mat*10**5))/10**5
    df = pd.DataFrame(mat)
    df.to_csv("input.txt", header = None,index = False )




def repeated_test(iter):
    for i in range (iter):
        n= random.randint(2,1001)
        m = random.randint(1,11)
        create_random_test(n,m)
        print (f"Iteration {i+1} out of {iter}. shape is {(n,m)}")
        subprocess.call([r"C:\Users\Tomer\CLionProjects\spectral_k_means\cmake-build-debug\SpectralKMeans.exe"])
        if not test():
            return None


if __name__== '__main__':
    repeated_test(50)


