from sympy import *
import numpy as np
x,y,z = symbols("x y z")
init_printing(use_unicode=True)

expr = (z-5) <= 4*z
print(expr)
print(type(expr))

print("check z=1 --> "+str(expr.subs(z, 1)))

print("expr left side --> "+str(expr.lhs))
print("expr right side --> "+str(expr.rhs))

expr = -x + 4 < 0

p = Poly(expr.lhs, x)
print("p: "+str(p))
print(solve_poly_inequality(p, '<'))

expr = x**2 >= 4
print("x^2 >= 4  ?")
print(solve_univariate_inequality(expr, x, relational=True))

print("expr: "+str(expr))
expr2 = expr.rhs - expr.lhs <= 0
print("expr2: "+str(expr2))

rel1 = Rel(1 , x, ">=")
rel1 = Ge(1 + 5*y + sqrt(2), -x)
print("rel1 print:")
print(type(rel1))
print(rel1)
print("arg[0] print:")
print(type(rel1.args[0]))
# print(sympify(rel1.args[0]).atoms().pop())
print(rel1.args[0])
print("rel op print:")
print(type(rel1.rel_op))
print(rel1.rel_op)
print("lhs print:")
print(type(rel1.lhs))
print(rel1.lhs)
print("reverse op:")
print(rel1.negated)
# print("x type = "+str(type(x)))
gen = rel1.free_symbols
coef_dic = {}
for s in gen:
    for e in preorder_traversal(rel1):
        if isinstance(e, Mul) and s in e.args:
            coef = 0;
            for val in e.args:
                if not isinstance(val, Symbol):
                    coef += val;
                    print("val type = " + str(type(val)))
            print("s = "+str(s)+", coef = "+str(coef)+" , e = "+str(e) + " , type of e = "+str(type(e)))
            coef_dic.update({s : coef})
            break
        elif isinstance(e, Symbol) and e is s:
            print("s = "+str(s)+", coef = 1" + " , e = "+str(e) + " , type of e = "+str(type(e)))
            coef_dic.update({s : 1})
            break
print(coef_dic)
# testi = 1*x
# print(type(x))

# x = sqrt(16)
# print(" x now is type = "+str(type(x)))
# print(x.evalf())


# obj = [-1, -2]
# #      ─┬  ─┬
# #       │   └┤ Coefficient for y
# #       └────┤ Coefficient for x
#
# lhs_ineq = [[ 2,  1],  # Red constraint left side
#             [-4,  5],  # Blue constraint left side
#             [ 1, -2]]  # Yellow constraint left side
#
# rhs_ineq = [20,  # Red constraint right side
#             10,  # Blue constraint right side
#              2]  # Yellow constraint right side
#
# lhs_eq = [[-1, 5]]  # Green constraint left side
# rhs_eq = [15]       # Green constraint right side
#
# bnd = [(0, float("inf")),  # Bounds of x
#        (0, float("inf"))]  # Bounds of

# opt = linprog(c=obj, A_ub=lhs_ineq, b_ub=rhs_ineq,
#               A_eq=lhs_eq, b_eq=rhs_eq, bounds=bnd,
#               method="revised simplex")

# obj = [-1, -2]
# #      ─┬  ─┬
# #       │   └┤ Coefficient for y
# #       └────┤ Coefficient for x
#
# lhs_ineq = [[ 2,  1],  # Red constraint left side
#             [-4,  5],  # Blue constraint left side
#             [ 1, -2]]  # Yellow constraint left side
#
# rhs_ineq = [20,  # Red constraint right side
#             10,  # Blue constraint right side
#              2]  # Yellow constraint right side
#
# # lhs_eq = [[-1, 5]]  # Green constraint left side
# # rhs_eq = [15]       # Green constraint right side
# #
# bnd = [(0, float("inf")),  # Bounds of x
#        (0, float("inf"))]  # Bounds of y

obj = [-1, -1]
#      ─┬  ─┬
#       │   └┤ Coefficient for y
#       └────┤ Coefficient for x

lhs_ineq = [[ 1,  1],
            [0, 1],
            [1, 0]]

rhs_ineq = [10,
            5,
            3]


bnd = [(0, float("inf")),  # Bounds of x
       (1, float("inf"))]  # Bounds of y

from scipy.optimize import linprog

opt = linprog(c=obj, A_ub=lhs_ineq, b_ub=rhs_ineq, bounds=bnd,
              method="revised simplex")
print(opt)
# x = Integer(0)
# sol = x + opt.fun
# print(sol*-1)
# c = [-1, 4]
# A = [[-3, 1], [1, 2]]
# b = [6, 4]
# x0_bounds = (None, None)
# x1_bounds = (-3, None)
#
# res = linprog(c, A_ub=A, b_ub=b, bounds=[x0_bounds, x1_bounds])

A = Matrix([[0, 1, 2], [-1, 0, sqrt(2)], [2, 1, 7]])
print(str(type(A[1,2])))

rel1 = 2*x + sqrt(2) - 3 + y >= 5
print(rel1)
print(srepr(rel1.lhs))
print(rel1.lhs.args)
cons = 0
for v in rel1.lhs.args:
    if isinstance(v, Symbol) or isinstance(v, Mul):
        continue
    else:
        cons += v
print(cons)


print("****")

arr = np.array([[1,2,3], [4,5,6]])
arr2 = arr
arr = None
print(arr)
