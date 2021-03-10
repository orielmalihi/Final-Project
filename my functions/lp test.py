from sympy import *
import random


# Pivot operation for simplex method
# http://web.tecnico.ulisboa.pt/mcasquilho/acad/or/ftp/FergusonUCLA_LP.pdf pp.16
def pivot(M, i, j):
    p = M[i, j]
    MM = Matrix.zeros(M.rows, M.cols)
    for ii in range(M.rows):
        for jj in range(M.cols):
            if ii == i and jj == j:
                MM[ii, jj] = 1 / M[i, j]
            elif ii == i and jj != j:
                MM[ii, jj] = M[ii, jj] / M[i, j]
            elif ii != i and jj == j:
                MM[ii, jj] = -M[ii, jj] / M[i, j]
            else:
                MM[ii, jj] = M[ii, jj] - M[ii, j] * M[i, jj] / M[i, j]
    return MM


# Simplex method with randomized pivoting
# http://web.tecnico.ulisboa.pt/mcasquilho/acad/or/ftp/FergusonUCLA_LP.pdf
def simplex(M, R, S):
    while True:
        B = M[:-1, -1]
        A = M[:-1, :-1]
        C = M[-1, :-1]
        if all(B[i] >= 0 for i in range(B.rows)):
            piv_cols = []
            for j in range(C.cols):
                if C[j] < 0:
                    piv_cols.append(j)

            if not piv_cols:
                return M, R, S
            random.shuffle(piv_cols)
            j0 = piv_cols[0]

            piv_rows = []
            for i in range(A.rows):
                if A[i, j0] > 0:
                    ratio = B[i] / A[i, j0]
                    piv_rows.append((ratio, i))

            if not piv_rows:
                assert False, 'Unbounded'
            piv_rows = sorted(piv_rows, key=lambda x: (x[0], x[1]))
            piv_rows = [(ratio, i) for ratio, i in piv_rows if ratio == piv_rows[0][0]]
            random.shuffle(piv_rows)
            _, i0 = piv_rows[0]

            M = pivot(M, i0, j0)
            R[j0], S[i0] = S[i0], R[j0]
        else:
            for k in range(B.rows):
                if B[k] < 0:
                    break

            piv_cols = []
            for j in range(A.cols):
                if A[k, j] < 0:
                    piv_cols.append(j)

            if not piv_cols:
                assert False, 'Infeasible'
            random.shuffle(piv_cols)
            j0 = piv_cols[0]

            ratio = B[k] / A[k, j0]
            piv_rows = [(ratio, k)]
            for i in range(A.rows):
                if A[i, j0] > 0 and B[i] > 0:
                    ratio = B[i] / A[i, j0]
                    piv_rows.append((ratio, i))

            piv_rows = sorted(piv_rows, key=lambda x: (x[0], x[1]))
            piv_rows = [(ratio, i) for ratio, i in piv_rows if ratio == piv_rows[0][0]]
            random.shuffle(piv_rows)
            _, i0 = piv_rows[0]

            M = pivot(M, i0, j0)
            R[j0], S[i0] = S[i0], R[j0]


# Maximize Cx + D constrained with Ax <= B and x >= 0
# Minimize y^{T}B constrained with y^{T}A >= C^{T} and y >= 0
def linear_programming(A, B, C, D):
    M = Matrix([[A, B], [-C, D]])
    r_orig = ['x_{}'.format(j) for j in range(M.cols - 1)]
    s_orig = ['y_{}'.format(i) for i in range(M.rows - 1)]

    r = r_orig.copy()
    s = s_orig.copy()

    M, r, s = simplex(M, r, s)

    argmax = []
    argmin_dual = []

    for x in r_orig:
        for i, xx in enumerate(s):
            if x == xx:
                argmax.append(M[i, -1])
                break
        else:
            argmax.append(S.Zero)

    for x in s_orig:
        for i, xx in enumerate(r):
            if x == xx:
                argmin_dual.append(M[-1, i])
                break
        else:
            argmin_dual.append(S.Zero)

    return M[-1, -1], argmax, argmin_dual
A = Matrix([[0, 1, 2], [-1, 0, -3], [2, 1, 7]])
B = Matrix([3, -2, 5])
C = Matrix([[1, 1, 5]])
D = Matrix([0])

# (11/3, [0, 1/3, 2/3], [0, 2/3, 1])
a,b,c = linear_programming(A, B, C, D)
print(a)
print(b)
print(c)


################################

obj = [-1,-1,-5]
#      ─┬  ─┬
#       │   └┤ Coefficient for y
#       └────┤ Coefficient for x

lhs_ineq = [[ 0,1,2],
            [-1,0,-3],
            [2,1,7]]

rhs_ineq = [3,
            -2,
            5]


bnd = [(0, float("inf")),  # Bounds of x
       (0, float("inf")),
       (0, float("inf"))]  # Bounds of y

from scipy.optimize import linprog

opt = linprog(c=obj, A_ub=lhs_ineq, b_ub=rhs_ineq, bounds=bnd,
              method="revised simplex")
print(opt)

s = sqrt(2)
print(s.evalf())


# arr = [1,2,3,-4]
# arr = Matrix(arr)
# print(arr*-1)
