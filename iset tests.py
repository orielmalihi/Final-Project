from sympy import *

def is_implied_by(iset, target):
    """

    :param iset: linear inequalities set
    :param target: target linear inequality
    :return: true if the trager is implied by iset or false
    if it is not.

>>> eq1 = 2*x + y <= 20
>>> eq2 = -4*x + 5*y <= 10
>>> eq3 = -x + 2*y >= -2
>>> eq4 = x >= 0
>>> eq5 = y >= 0
>>> is_implied_by([eq1, eq2, eq3, eq4, eq5], x + 2*y > 10)  # the optinal values of x+2y are between [0, 20.714285714285715]
False
>>> is_implied_by([eq1, eq2, eq3, eq4, eq5], x + 2*y > -1)
True
    """

def simplify_linear_inequalites(iset):
    """
    :param iset: linear inequalities set to be simplified ( redaundent inequalities
    would be removed.
    :return: simplified set of inequalities.
>>> eq1 = x + y > 20
>>> eq2 = x + y > 10
>>> simplify_linear_inequalites([eq1, eq2])
[x + y > 10 ]
    """

def find_values_interval(iset, expr):
    """
    :param iset: linear inequalities set that would be the constraints.
    :param expr: target expression
    :return: the optional values for the expression given the constraints.
>>> eq1 = 2*x + y <= 20
>>> eq2 = -4*x + 5*y <= 10
>>> eq3 = -x + 2*y >= -2
>>> eq4 = x >= 0
>>> eq5 = y >= 0
>>> find_valus_interval([eq1, eq2, eq3, eq4, eq5], x + 2*y)
[0, 20.714285714285715]
>>> eq1 = x + y < 10
>>> eq2 = y <= 8
>>> eq3 = x >= 0
>>> eq4 = y >= 2
>>>> find_valus_interval([eq1, eq2, eq3, eq4], x + y)
[2, 10]
>>> eq1 = x + y < 10
>>> eq2 = y <= 5
>>> eq3 = y >= 1
>>> eq4 = x <= 3
>>> eq5 = x >= 0
>>>> find_valus_interval([eq1, eq2, eq3, eq4, eq5], x + y)
[1, 8]

    """
