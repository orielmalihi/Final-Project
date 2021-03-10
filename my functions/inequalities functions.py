from sympy import *
import random
from typing import List

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


# A = Matrix([[0, 1, 2], [-1, 0, -3], [2, 1, 7]])
# B = Matrix([3, -2, 5])
# C = Matrix([[1, 1, 5]])
# D = Matrix([0])

# (11/3, [0, 1/3, 2/3], [0, 2/3, 1])
# a,b,c = linear_programming(A, B, C, D)
# print(a)
# print(b)
# print(c)

x,y,z = symbols("x y z")


def rearrange_ineq(ineq : Rel):
    """

    Parameters:
    inequality to rearrange.

    Returns
    rearranged inequality so that the symbols are on the left side
    and the constants are on the right side.

    >>> e = 2 - x  >= y + 4
    >>> rearrange_ineq(e)
    x + y <= -2
    >>> e = (1/2)*x - sqrt(2) - 4 <= 0
    >>> rearrange_ineq(e)
    0.5*x <= sqrt(2) + 4
    >>> e = 3*x + 2*y - 8 - sqrt(2) >= x + 4
    >>> rearrange_ineq(e)
    -2*x - 2*y <= -12 - sqrt(2)

    """
    nl = ineq.lhs - ineq.rhs
    cons = 0
    op = ineq.rel_op
    if op != ">=" and op != "<=":
        assert False, 'mismatch relation'
    if len(nl.free_symbols) == 1 and isinstance(nl, Mul) and op == "<=":
        return Le(nl, 0)
    if len(nl.free_symbols) == 1 and isinstance(nl, Mul) and op == ">=":
        return Le(-1*nl, 0)
    for v in nl.args:
        if len(v.free_symbols) >= 1:
            continue
        else:
            cons += v
    if op == ">=":
        return Le(-1*(nl-cons), cons)
    elif op == "<=":
        return Le(nl - cons, -1*cons)


def get_coef(expr: Expr):
    temp_dic = {}
    if len(expr.free_symbols) == 1 and isinstance(expr, Symbol):
        temp_dic[expr] = 1
        return temp_dic
    if len(expr.free_symbols) == 1 and isinstance(expr, Mul):
        if isinstance(expr.args[0], Symbol):
            temp_dic[expr.args[0]] = expr.args[1]
            return temp_dic
        elif isinstance(expr.args[1], Symbol):
            temp_dic[expr.args[1]] = expr.args[0]
            return temp_dic
    for v in expr.args:
        if isinstance(v, Symbol):
            temp_dic[v] = temp_dic.get(v, 0) + 1
        if isinstance(v, Mul):
            if isinstance(v.args[0], Symbol):
                temp_dic[v.args[0]] = temp_dic.get(v.args[0], 0) + v.args[1]
            elif isinstance(v.args[1], Symbol):
                temp_dic[v.args[1]] = temp_dic.get(v.args[1], 0) + v.args[0]
    return temp_dic



def find_valus_interval(iset: List[Rel], expr: Expr):
    _symbols = set()
    for r in iset:
        _symbols = _symbols.union(r.free_symbols)
    al = []
    bl = []
    cl = []
    dl = []
    for r in iset:
        r = rearrange_ineq(r)
        tl = []
        temp_dic = get_coef(r.lhs)
        for s in _symbols:
            tl.append(temp_dic.get(s, 0))
        al.append(tl)
        bl.append(r.rhs)
        print(temp_dic)
    print("al: ", al)
    print("bl: ", bl)
    temp_dic = get_coef(expr)
    for s in _symbols:
        cl.append(temp_dic.get(s, 0))
    print("cl: ", cl)
    if len(expr.free_symbols) == 1 and (isinstance(expr, Mul) or isinstance(expr, Symbol)):
        dl.append(0)
    else:
        dl.append(0)
        cons = 0
        for v in expr.args:
            if len(v.free_symbols) >= 1:
                continue
            else:
                cons += v
        dl[0] = cons
    print("dl: ", dl)
    A = Matrix(al)
    B = Matrix(bl)
    C = Matrix([cl])
    D = Matrix(dl)
    max,b,c = linear_programming(A, B, C, D)
    min, b,c = linear_programming(A, B, C*-1, D*-1)
    return [-1*min, max]
    # print(a)
    # print(b)
    # print(c)





eq1 = 2*x + y <= 20
eq2 = 4*x + 5*y <= 10
eq3 = -x + 2*y >= -2
# print(type(eq1))

ans = find_valus_interval([eq1, eq2, eq3], x + 2*y)
print(ans)

# e = 3*x + 2*y - 8 - sqrt(2) >= x + 4
# e = rearrange_ineq(e)
# print(e)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
