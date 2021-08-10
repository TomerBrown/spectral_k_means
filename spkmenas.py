import pandas as pd
import numpy as np
import math


def l2_norm(point1, point2):
    """Computes the l2 norm of two given points"""
    result = np.sqrt(np.sum(((point1 - point2) ** 2)))
    return result


def weighted_adj_mat(matrix):
    """ Computes the matrix W given it's definition as exponent of half of l2 norm """
    n, m = matrix.shape
    W = np.zeros((n, n))
    for i in range(n):
        # It is  written to go only on upper triangle of matrix and copy the rest for efficiency
        for j in range(i, n):
            if (i != j):
                point1 = matrix[i, :]
                point2 = matrix[j, :]
                W[i, j] = math.exp(-l2_norm(point1, point2) / 2)
                # Just copying the value for symmetric
                W[j, i] = W[i, j]
    return W


def diagonal_degree_mat(W):
    """Given W - The Wighted Adjacency matrix, computes it's diagonal degree matrix - D"""
    n, m = W.shape
    D = np.zeros((n, n))
    for i in range(n):
        D[i, i] = np.sum(W[i, :])
    return D


def pow_diag(D, pow):
    """ given a matrix and a real num - pow, Changes D so with D_new[i,i] = D[i,i]**pow only on diagonal entries"""
    n, m = D.shape
    for i in range(n):
        D[i, i] = D[i, i] ** pow


def l_norm_mat(W, D_half):
    """Computes the Lnorm Matrix given it's definition"""
    n, m = W.shape
    return np.identity(n) - D_half @ W @ D_half


def max_off_diagonal(A):
    """ :returns tuple (x,y) so A[x,y] is the biggest (by absolute value) of all non- diagonal entries of A"""
    n, m = A.shape
    max_val = abs(A[1, 0])
    x, y = 1, 0
    for i in range(n):
        for j in range(m):
            num = A[i, j]
            if (i != j and abs(num) > abs(max_val)):
                max_val = abs(num)
                x = i
                y = j
    return (x, y)


def sign(num):
    """:returns -1 if num<0 else 1"""
    if (num < 0):
        return -1
    else:
        return 1


def create_p(A):
    """:returns the matrix P corresponding to the given matrix A"""
    n, m = A.shape
    i, j = max_off_diagonal(A)
    theta = (A[i, j] - A[i, i]) / (2 * A[i, j])
    t = sign(theta) / (abs(theta) + math.sqrt(theta ** 2 + 1))
    c = 1 / (math.sqrt(t ** 2 + 1))
    s = t * c
    p = np.identity(n)
    p[i, i] = c
    p[j, j] = c
    p[i, j] = s
    p[j, i] = -s
    return p


def is_diag(A):
    """:returns whether the given matrix A is diagonal """
    n,m = A.shape
    for i in range(n):
        for j in range(m):
            if (i!=j and A[i,j]!=0):
                return False
    return True


def off(A):
    """Implementing the off of matrix as described in project """
    n,m = A.shape
    frob = np.sum(np.sum(np.power(A,2)))
    diag = 0
    for i in range (n):
        diag += A[i,i]**2
    return frob-diag


def jacobi_algo(A , epsilon = 0.001):
    n,m = A.shape
    p  = create_p(A)
    V= p
    A_prime = (p.T@A)@p
    print(off(A))
    print(off(A_prime))
    while (off(A)-off(A_prime)>epsilon):
        c+=1
        print(c)
        A = A_prime
        p  = create_p(A)
        V = V@p
        A_prime = p.T@A@p
    eigval = []
    for i in range (n):
        eigval.append(A[i,i])
    return (eigval,V)




def test_max_d():
    arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    print(max_off_diagonal(arr))
    print((is_diag(arr)))
    arr2 = np.array([[1,2],[3,4]])
    print(off(arr2))


path = "tests/input_1.txt"

# Load File and convert it to numpy array
df = pd.read_csv(path, header=None)
matrix = df.to_numpy()

# Calculate W
W = weighted_adj_mat(matrix)

# Calculate D
D = diagonal_degree_mat(W)

# Calculate D_half
pow_diag(D, -0.5)
D_half = D

l_norm = l_norm_mat(W, D_half)
eigval , eigvectors = jacobi_algo(l_norm)
eigval = sorted(eigval)
print(eigval)

