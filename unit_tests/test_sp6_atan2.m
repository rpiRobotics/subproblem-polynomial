syms x1 x2 real

[P, S] = subproblem_setups.sp_6.setup;

[p1, p2] = half_tan_sp6(P.H, P.K, P.P, P.d1, P.d2, [x1 x2]);
double(subs(p1, [x1 x2], tan(1/2*[S.theta1 S.theta2])))
double(subs(p2, [x1 x2], tan(1/2*[S.theta1 S.theta2])))