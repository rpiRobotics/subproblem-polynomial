
syms x1 x2 x3 x6 x7 real
R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_23 = half_tan_rot(kin.H(:,3), x3);
R_56 = half_tan_rot(kin.H(:,6), x6);
R_67 = half_tan_rot(kin.H(:,7), x7);

p_07 = p_0T - R_07 * kin.P(:,8);
[~, n_SEW] = SEW.inv_kin([0;0;0], p_07, psi);

% q6 (SP4 based on q7 search)
[~, a, b, c] = half_tan_sp4(R_67*R_07'*n_SEW, kin.P(:,6), -kin.H(:,6), n_SEW'*p_07);
lhs_x6 = a*x6^2 + b*x6 + c;

p_WE = R_07 * R_67' * R_56' * kin.P(:,6);

% q1
[~, a, b, c] = half_tan_sp3(kin.P(:,2), p_07-p_WE, kin.H(:,1), norm(kin.P(:,4)));
lhs_x1 = a*x1^2 + b*x1 + c;

% Find (q_2, q_3) with Subproblem 2
[~, ~, a3, b3, c3, a2, b2, c2] = half_tan_sp2(kin.P(:,4), R_01'*p_07-R_01'*p_WE-kin.P(:,2), kin.H(:,3), -kin.H(:,2));
lhs_x3 = a3*x3^2 + b3*x3 + c3;
lhs_x2 = a2*x2^2 + b2*x2 + c2;

% Error eqn
lhs_err = kin.H(:,4)' * (R_01 * R_12 * R_23)' * R_07 * (R_56 * R_67)'* kin.H(:,5) - kin.H(:,4)' * kin.H(:,5);

%%
eqns_lhs = [lhs_x1; lhs_x2; lhs_x3; lhs_x6; lhs_err];
eqns_lhs = simplify(eqns_lhs)

eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)

filename = "sawyer_polynomials_alt.txt";

fileID = fopen(filename,'w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i := ", i);
    % fprintf(fileID, string(vpa(eqns_num(i))) + '\n');
    fprintf(fileID, string(eqns_num(i)) + ':\n\n');
end
