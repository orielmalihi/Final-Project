from sympy import *
import random
from typing import List


def _pivot(M, i, j):
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


def _simplex(M, R, S):
    """
    Simplex method with randomized pivoting
    http://web.tecnico.ulisboa.pt/mcasquilho/acad/or/ftp/FergusonUCLA_LP.pdf
    """
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

            M = _pivot(M, i0, j0)
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

            M = _pivot(M, i0, j0)
            R[j0], S[i0] = S[i0], R[j0]



def linear_programming(A, B, C, D):
    """
    When x is a column vector of variables, y is a column vector of dual variables,
    and when the objective is either:

    *) Maximizing Cx constrained to Ax <= B and x >= 0.
    *) Minimizing y^{T}B constrained to y^{T}A >= C^{T} and y >= 0.

    Thw method s eturns a triplet of solutions optimum, argmax, argmax_dual where
    optimum is the optimum solution, argmax is x and argmax_dual is y.

    Examples
    ========

    >>> A = Matrix([[0, 1, 2], [-1, 0, -3], [2, 1, 7]])
    >>> B = Matrix([3, -2, 5])
    >>> C = Matrix([[1, 1, 5]])
    >>> D = Matrix([0])
    >>> linear_programming(A, B, C, D)
    (11/3, [0, 1/3, 2/3], [0, 2/3, 1])

    >>> A = Matrix([[2, 1, 3], [-1, -2, 3], [1, 1, 5]])
    >>> linear_programming(A, B, C, D)
    (35/9, [0, 5/3, 4/9], [13/9, 2/9, 0])
    """

    M = Matrix([[A, B], [-C, D]])
    r_orig = ['x_{}'.format(j) for j in range(M.cols - 1)]
    s_orig = ['y_{}'.format(i) for i in range(M.rows - 1)]

    r = r_orig.copy()
    s = s_orig.copy()

    M, r, s = _simplex(M, r, s)

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
        assert False, 'Error: expected <= or >= realtion, got ' + str(op)
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
    for s in expr.free_symbols:
        temp_dic[s] = 0
    if len(expr.free_symbols) == 1 and isinstance(expr, Symbol):
        temp_dic[expr] = 1
        return temp_dic
    if len(expr.free_symbols) == 1 and isinstance(expr, Mul):
        smb = None
        mul = 1
        for v in expr.args:
            if isinstance(v, Symbol):
                smb = v
            else:
                mul *= nsimplify(v)
        temp_dic[smb] = mul
        return temp_dic
    for v in expr.args:
        if isinstance(v, Symbol):
            temp_dic[v] = temp_dic.get(v, 0) + 1
        if isinstance(v, Mul):
            smb = None
            mul = 1
            for val in v.args:
                if isinstance(val, Symbol):
                    smb = val
                else:
                    mul *= nsimplify(val)
            temp_dic[smb] += mul
    return temp_dic


# Maximize and Minimize Cx + D (expr) constrained with Ax <= B and x >= 0 (iset), returning [min, max] interval
def find_valus_interval(iset: List[Rel], expr: Expr):
    """
    Parameters
    ----------
    iset - linear inequalities set of constraints
    expr - target expression

    Returns the the possible values of expr as an interval of [min, max].
    -------
    >>> eq1 = x <= 2
    >>> eq2 = y <= 1
    >>> find_valus_interval([eq1, eq2], x)
    [0, 2]
    >>> find_valus_interval([eq1, eq2], y)
    [0, 1]
    >>> find_valus_interval([eq1, eq2], y + 2*x)
    [0, 5]

    >>> eq1 = z >= 2
    >>> eq2 = z <= sqrt(5)
    >>> eq3 = x + z <= 3
    >>> find_valus_interval([eq1, eq2, eq3], z)
    [2, sqrt(5)]
    >>> find_valus_interval([eq1, eq2, eq3], x)
    [0, 1]
    >>> find_valus_interval([eq1, eq2, eq3], 2*x + z) # we want more value to x and less to z, so z = 2 and x = 1.
    [2, 4]
    >>> find_valus_interval([eq1, eq2, eq3], x + 2*z) # we want more value to z and less to x, so z = sqrt(5) and x = 3 - sqrt(5)
    [4, sqrt(5) + 3]

    >>> eq1 = sqrt(3)*z >= 2
    >>> eq2 = 1.0*z <= sqrt(5)
    >>> eq3 = 1/2*x + 0.5*x + z <= 3
    >>> find_valus_interval([eq1, eq2, eq3], z )
    [2*sqrt(3)/3, sqrt(5)]
    >>> find_valus_interval([eq1, eq2, eq3], x)
    [0, 3 - 2*sqrt(3)/3]
    >>> find_valus_interval([eq1, eq2, eq3], z - x)
    [-3 + 4*sqrt(3)/3, sqrt(5)]

    >>> eq1 = x + y <= 2
    >>> eq2 = 2*x - y <= 6
    >>> eq3 = -y >= -sqrt(2)
    >>> find_valus_interval([eq1, eq2, eq3], y)
    [0, sqrt(2)]
    >>> find_valus_interval([eq1, eq2, eq3], x)
    [0, 2]
    >>> find_valus_interval([eq1, eq2, eq3], x + 2*y)
    [0, sqrt(2) + 2]
    >>> eq1 = x <= 6
    >>> find_valus_interval([eq1], x)
    [0, 6]
    >>> eq1 = x + y <= 3
    >>> eq2 = x + y >= 2
    >>> find_valus_interval([eq1, eq2], x + y)
    [2, 3]
    >>> find_valus_interval([eq1, eq2], -x - y)
    [-3, -2]
    """
    _symbols = set()
    for r in iset:
        _symbols = _symbols.union(r.free_symbols)
    al = []
    bl = []
    cl = []
    dl = []
    for r in iset:
        r = rearrange_ineq(r)
        # print("rearenged: ", r)
        tl = []
        temp_dic = get_coef(r.lhs)
        # print(temp_dic)
        for s in _symbols:
            tl.append(temp_dic.get(s, 0))
        al.append(tl)
        bl.append(r.rhs)
    # print("al: ", al)
    # print("bl: ", bl)
    temp_dic = get_coef(expr)
    for s in temp_dic:
        if s not in _symbols:
            assert False, 'Error: there is no data about Symbol ' + str(s)
    for s in _symbols:
        cl.append(temp_dic.get(s, 0))
    # print("cl: ", cl)
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
    # print("dl: ", dl)
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

def is_implied_by(iset: List[Rel], target: Rel):
    """

    :param iset: linear inequalities set
    :param target: target linear inequality
    :return: true if the traget is implied by iset or false
    if it is not.

    >>> eq1 = 2*x + y <= 20
    >>> eq2 = -4*x + 5*y <= 10
    >>> eq3 = -x + 2*y >= -2
    >>> eq4 = x >= 0
    >>> eq5 = y >= 0
    >>> is_implied_by([eq1, eq2, eq3, eq4, eq5], x + 2*y >= 10)  # the optional values of x+2y are between [0, 145/7]
    False
    >>> is_implied_by([eq1, eq2, eq3, eq4, eq5], x + 2*y >= -0.00001)
    True
    >>> is_implied_by([eq1, eq2, eq3, eq4, eq5], x + y >= 0.00001 - y)
    False
    >>> is_implied_by([eq1, eq2, eq3, eq4, eq5], x + 2*y <= 145/7 + 0.0001)
    True
    >>> is_implied_by([eq1, eq2, eq3, eq4, eq5], x + 2*y <= 145/7 - 0.0001)
    False

    >>> eq1 = z >= 2
    >>> eq2 = z <= sqrt(5)
    >>> eq3 = x + z <= 3
    >>> is_implied_by([eq1, eq2, eq3], 2*x + z >= 1.9999) # the optinal values of 2*x + z are between [2, 4]
    True
    >>> is_implied_by([eq1, eq2, eq3], 2*x + z >= 2.00001)
    False
    >>> is_implied_by([eq1, eq2, eq3], 2*x + z <= 4.00001)
    True
    >>> is_implied_by([eq1, eq2, eq3], 2*x + z <= 3.9999)
    False
    """
    target = rearrange_ineq(target)
    min, max = find_valus_interval(iset, target.lhs)
    rel = target.rel_op
    if rel == ">=" and min >= target.rhs:
        return True
    if rel == ">=" and min < target.rhs:
        return False
    if rel == "<=" and max <= target.rhs:
        return True
    if rel == "<=" and max >= target.rhs:
        return False
    return None


def simplify_linear_inequalites(iset: List[Rel]):
    """
    :param iset: linear inequalities set to be simplified ( redaundent inequalities
    would be removed).
    :return: simplified set of inequalities.
    >>> eq1 = x + y >= 10
    >>> eq2 = x + y >= 20
    >>> eq3 = x + y <= 30
    >>> simplify_linear_inequalites([eq1, eq2, eq3])
    [x + y >= 20, x + y <= 30]
    >>> eq4 = x + y >= 5
    >>> eq5 = x + y <= 40
    >>> simplify_linear_inequalites([eq1, eq2, eq3, eq4, eq5])
    [x + y >= 20, x + y <= 30]
    >>> eq1 = x + y >= 10
    >>> eq2 = x - y >= 20
    >>> eq2 = 2*x >= 30
    >>> eq4 = x + y <= 30
    >>> simplify_linear_inequalites([eq1, eq2, eq3, eq4])
    [2*x >= 30, x + y <= 30]
    """
    # reduced = []
    # for target in iset:
    #     temp = []
    #     for ineq in iset:
    #         if ineq is not target:
    #             temp.append(ineq)
    #     try:
    #         if not is_implied_by(temp, target):
    #             reduced.append(target)
    #     except:
    #         reduced.append(target)
    reduced = []
    tiset = iset
    while len(tiset) != 0:
        target = tiset[0]
        temp = reduced + tiset[1:]
        try:
            if not is_implied_by(temp, target):
                reduced.append(target)
        except:
            reduced.append(target)
        tiset = tiset[1:]
    return reduced




# eq1 = 2*x + y <= 20
# eq2 = 4*x + 5*y <= 10
# eq3 = -x + 2*y >= -2
# print(type(eq1))

# ans = find_valus_interval([eq1, eq2, eq3], x + 2*y)
# print(ans)

# e = 3*x + 2*y - 8 - sqrt(2) >= x + 4
# e = rearrange_ineq(e)
# print(e)

# eq1 = x  >= 2
# eq2 = x <= 5
# ans = find_valus_interval([eq1, eq2], -x)
# print("ans: ", ans)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
