function [eqns_num, plot] = three_pairs_intersecting(kin, R_06, p_0T, filename)

syms x1 x2 x3 x4 real
R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_23 = half_tan_rot(kin.H(:,3), x3);
R_34 = half_tan_rot(kin.H(:,4), x4);

% Error eqn
R_04 = R_01 * R_12 * R_23 * R_34;
lhs_err = kin.H(:,5)' * R_04.' * R_06* kin.H(:,6) - kin.H(:,5)'*kin.H(:,6);

% x3 eqn
p_16 = p_0T - kin.P(:,1) - R_06*kin.P(:,7);
R_34 = half_tan_rot(kin.H(:,4), x4);
[plot.x3, a, b, c] = half_tan_sp3(R_34*kin.P(:,5), -kin.P(:,3), kin.H(:,3), norm(p_16));
lhs_x3 = a*x3^2 + b*x3 + c;

% x1 and x2 eqn
[plot.x1, plot.x2, a1, b1, c1, a2, b2, c2] = half_tan_sp2(p_16, kin.P(:,3)+R_23*R_34*kin.P(:,5), -kin.H(:,1), kin.H(:,2));

lhs_x1 = a1*x1^2 + b1*x1 + c1;
lhs_x2 = a2*x2^2 + b2*x2 + c2;

eqns_lhs = [lhs_err; lhs_x1; lhs_x2; lhs_x3];
eqns_lhs = simplify(eqns_lhs)


eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)

fileID = fopen(filename,'w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    % fprintf(fileID, string(vpa(eqns_num(i))) + '\n');
    fprintf(fileID, string(eqns_num(i)) + '\n');
end


% These equtions can be used for plotting error function
plot.x3 = simplify(plot.x3);
plot.x1 = simplify(plot.x1);
plot.x2 = simplify(plot.x2);
plot.err = eqns_lhs(1);

end