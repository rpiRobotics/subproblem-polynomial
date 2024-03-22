# Sawyer 7-DOF IK with SEW angle

`sawyer.m`: Generate the 4 multivariate polynomials and export them. Given the zeros of the polynomial from Maple, find the remaining joint angles.

`sawyer_polynomials.txt`: System of 4 multivariate polynomials.

`sawyer.mw`: Maple worksheet to find and solve resultant univariate polynomial.

`sawyer.mpl`: Same as above but in plaintext.

`sawyer_maple_output.txt`: Results from Maple.

`IK_R_2R_2R_2R_q7.m`: Function - IK using a 1D search over q7. Polynomial method should match this.

