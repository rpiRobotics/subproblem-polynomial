syms x1 x2 x3 real

[P, S] = subproblem_setups.sp_5.setup

[p1, p2, p3] = half_tan_sp5(P.p0,P.p1,P.p2,P.p3,P.k1,P.k2,P.k3, [x1 x2 x3])
%%
double(subs(p1, [x1 x2 x3], tan(1/2*[S.theta1 S.theta2 S.theta3])))
double(subs(p2, [x1 x2 x3], tan(1/2*[S.theta1 S.theta2 S.theta3])))
double(subs(p3, [x1 x2 x3], tan(1/2*[S.theta1 S.theta2 S.theta3])))

%%
[p1, p2, R] = half_tan_sp5_R(P.p0,P.p1,P.p2,P.p3,P.k1,P.k2,P.k3, [x1 x3]);
double(subs(p1, [x1 x3], tan(1/2*[S.theta1 S.theta3])))
double(subs(p2, [x1 x3], tan(1/2*[S.theta1 S.theta3])))
double(subs(rot(P.k2, S.theta2)'*R, [x1 x3], tan(1/2*[S.theta1 S.theta3])))