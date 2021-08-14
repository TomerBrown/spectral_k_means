import pandas as pd
import numpy as np
import random;


def create_random_test(n= random.randint(2,1001),m = random.randint(1,11)):
    mat = (np.random.rand(n,m)-1)*10
    mat = np.round(mat*10**5)/10**5
    df = pd.DataFrame(mat)
    df.to_csv("input.txt", header = None,index = False )
