function eqns_num = two_intersecting(kin, R_06, p_06, filename)
syms x1  x3 x4 real
R_01 = half_tan_rot(kin.H(:,1), x1);
% R_12 = half_tan_rot(kin.H(:,2), x2);
R_23 = half_tan_rot(kin.H(:,3), x3);
R_34 = half_tan_rot(kin.H(:,4), x4);

p_35_3 = kin.P(:,4) + R_34*kin.P(:,5);

[lhs_1, lhs_2, R_12] = half_tan_sp5_R( ...
               -kin.P(:,2), p_06, kin.P(:,3), p_35_3, ...
               -kin.H(:,1), kin.H(:,2), kin.H(:,3), [x1 x3]);

% Error eqn
R_04 = R_01 * R_12 * R_23 * R_34;
lhs_err = kin.H(:,5)' * R_04' * R_06* kin.H(:,6) - kin.H(:,5)'*kin.H(:,6);


% lhs_err = expand(lhs_err)
eqns_lhs = [lhs_err; lhs_1; lhs_2];
eqns_lhs = simplify(eqns_lhs);

eqns_frac = simplifyFraction(eqns_lhs);
eqns_num = numden(eqns_frac);

fileID = fopen(filename,'w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    fprintf(fileID, string(eqns_num(i)) + '\n');
end
end