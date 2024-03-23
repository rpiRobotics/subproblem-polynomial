# FANUC CRX-10iA/L 6-DOF IK

Three pairs of intersecting axes

(This code works as long as there are two pairs of intersecting axes)

`crx.m`: Generate the 4 multivariate polynomials and export them. Given the zeros of the polynomial from Maple, find the remaining joint angles.

'crx.mpl': Maple file to find and solve resultant univariate polynomial.

'crx_out.txt': Output from Maple