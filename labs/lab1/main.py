# from intervalpy import Interval
# from intervalpy import Interval
from typing import List
import math as m
from Interval import *


def det(A: List[List[Interval]]):
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


def create_int_matrix(A: List[List[int]], delta: int | float):
    A_int = [[Interval(A[0][0] - delta, A[0][0] + delta), Interval(A[0][1], A[0][1])],
             [Interval(A[1][0] - delta, A[1][0] + delta), Interval(A[1][1], A[1][1])]]
    return A_int


def min_delta_search(A: List[List[int]], delta: int | float) -> bool:
    def f(A: List[List[int]], delta: int | float):
        A_int = create_int_matrix(A, delta)
        det_A = det(A_int)
        if 0.0 in det_A:
            return True
        else:
            return False

    min_delta = delta
    a = 0
    b = delta
    while not m.isclose(a, b, rel_tol=1e-14):
        min_delta = (a + b) / 2
        if f(A, min_delta):
            b = min_delta
        else:
            a = min_delta
    return (a + b) / 2


if __name__ == "__main__":
    print("Enter A_11, A_12, A_21, A_22: ", end='')
    a, b, c, d = [float(n) for n in input().split()]
    if (a < 0 or b < 0 or c < 0 or d < 0):
        raise ValueError("Values must be non-negative.")
    delta = min(a, b, c, d)
    print("delta = ", delta)
    A_orig = [[a, b], [c, d]]
    A = [[Interval(a - delta, a + delta), Interval(b, b)],
         [Interval(c - delta, c + delta), Interval(d, d)]]

    det_A = det(A)

    print("det_A = ", det_A)
    if m.isclose(det_A.mid, 0.0):
        print("delta = 0")
    elif 0.0 in Interval(det_A.left, det_A.right):
        min_delta = min_delta_search(A_orig, delta)
        print("min delta = ", min_delta)
    print(det(create_int_matrix(A_orig, min_delta)))