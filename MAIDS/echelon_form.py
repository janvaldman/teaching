import numpy as np
import sympy as sym

A1 = np.matrix([[1, 1, -1, -1], [2, 5, -7, -5], [2, -1, 1, 3], [5, 2, -4, -2]])
b1 = np.matrix([1, -2, 4, 6]).T

A2 = np.matrix([[1, -1, 0, 0, 1], [1, 1, 0, -3, 0], [2, -1, 0, 1, -1],
                [-1, 2, 0, -2, -1]])
b2 = np.matrix([3, 6, 5, -1]).T

A = A2
b = b2

Ab = np.hstack((A, b))
Abech = sym.Matrix(Ab).rref()

print("A:\n", A, "\n")
print("b (transpose):\n", b.T, "\n")
print("Ab:\n", Ab, "\n")
print("Ab - reduced echelon form:\n", np.matrix(Abech[0]), "\n")
print("basic variables: \n", Abech[1], "\n")

xp = np.linalg.pinv(A).dot(b)
print("xp - particular solution (transpose): \n", xp.T, "\n")
print("check if A*xp-b=0 (transpose): \n", (A * xp - b).T)