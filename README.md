# subproblem-polynomial

Whenever search is used with subproblems to solve kinematics problems, an alternative method can be used that guarantees finding all solutions. Each subproblem can be written as a multivariate polynomial in $x_i = \tan(q_i/2)$. The error function can also be written as a similar polynomial. Then, the resultant of this system of polynomials can be found to eliminate all but one variable. The zeros of this polynomial can be easily obtained, and these zeros correspond to the zeros of the error function. Subproblems can then be used to find the remaining joint angles.

For closed-form solutions for robots with three intersecting or parallel axes, see [ik-geo](https://github.com/rpiRobotics/ik-geo). Code in this repo depends on that code.


## Folder Breakdown

`+polynomial_IK`: Functions to generate system of polynomial equations, and functions to find remaining joint angles given 1 or 2 joint solutions

`examples`: Code to generate systems of polynomials for specific robots

`maple`: Code to solve systems of polynoimals in Maple (Reference code)

`mathematica`: Code to solve systems of polynomials in Mathematica (Not as good as Maple code)

`polynomial_subproblems`: Subproblem solutions and rotation matrix in terms of tangent half angle

`robot_definitions`: Kinematic parameters for example robots

`unit_tests`: Code to test certain MATLAB files