# Robotics Research Corporation (RRC) K-1207i (with fixed q_6) 6-DOF IK

Two intersecting axes

(This robot also has two non-consecutive intersecting axes)

`RRC.m`: Generate the 3 multivariate polynomials and export them. Given the zeros of the polynomial from Maple, find the remaining joint angles.

`RRC_fixed_q6_eqns.txt: System of 3 multivariate polynomials.

`automated_RRC_fixed_q6_rational.mpl`: Maple file to find and solve resultant univariate polynomial.

`automated_RRC_fixed_q6_rational_out.txt`: Output from Maple.