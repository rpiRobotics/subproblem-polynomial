function eqns = three_pairs_intersecting_indep(kin, R_06, p_0T)

syms x1 x2 x3 x4 real
R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_23 = half_tan_rot(kin.H(:,3), x3);
R_34 = half_tan_rot(kin.H(:,4), x4);

% Error eqn
R_04 = R_01 * R_12 * R_23 * R_34;
eqns.err = kin.H(:,5)' * R_04.' * R_06* kin.H(:,6) - kin.H(:,5)'*kin.H(:,6);

% x3 eqn
p_16 = p_0T - kin.P(:,1) - R_06*kin.P(:,7);
R_34 = half_tan_rot(kin.H(:,4), x4);
eqns.x3 = half_tan_sp3(R_34*kin.P(:,5), -kin.P(:,3), kin.H(:,3), norm(p_16));

% x1 and x2 eqn
[eqns.x1, eqns.x2] = half_tan_sp2(p_16, kin.P(:,3)+R_23*R_34*kin.P(:,5), -kin.H(:,1), kin.H(:,2));

eqns.err = simplify(eqns.err);
eqns.x3 = simplify(eqns.x3);
eqns.x1 = simplify(eqns.x1);
eqns.x2 = simplify(eqns.x2);

end