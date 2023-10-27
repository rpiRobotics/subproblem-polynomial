function eqns_num = general_6R(kin, R_06, p_06, filename)
syms x1 x2 x3 x5 real
R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_23 = half_tan_rot(kin.H(:,3), x3);
% R_34 = half_tan_rot(kin.H(:,4), x4);
R_45 = half_tan_rot(kin.H(:,5), x5);

[lhs_1, lhs_2, R_34] = half_tan_sp5_R( ...
    -kin.P(:,4), ...
    R_12'*(R_01'*p_06-kin.P(:,2))-kin.P(:,3), ...
    kin.P(:,5), ...
    kin.P(:,6), ...
    -kin.H(:,3), kin.H(:,4), kin.H(:,5), [x3 x5]);

% Error eqns
% R_05 = R_01 * R_12 * R_23 * R_34 * R_45;
%lhs_err = R_05 * kin.H(:,6) - R_06 * kin.H(:,6); % 3 equations
%lhs_err = (R_06 * kin.H(:,6))'*  R_05 * kin.H(:,6) -1
lhs_err_1 = kin.H(:,3)'*R_34*R_45*kin.H(:,6) - kin.H(:,3)'*(R_01*R_12)'*R_06*kin.H(:,6); % No x3
lhs_err_2 = kin.H(:,4)'*R_45*kin.H(:,6) - kin.H(:,4)'*(R_01*R_12*R_23)'*R_06*kin.H(:,6); % No x4

% lhs_err_3 = kin.H(:,1)'*R_12 * R_23 * R_34 * R_45 * kin.H(:,6) - kin.H(:,1)'*R_06*kin.H(:,6); % no x1
% lhs_err_4 = kin.H(:,2)'*R_23 * R_34 * R_45 * kin.H(:,6) - kin.H(:,2)'* R_01' * R_06*kin.H(:,6); % no x2
% lhs_err_5 = kin.H(:,5)'*kin.H(:,6) - kin.H(:,5)'*(R_01*R_12*R_23*R_34)'*R_06*kin.H(:,6); % No x5

eqns_lhs = [lhs_err_1; lhs_err_2; lhs_1; lhs_2];

eqns_lhs = simplify(eqns_lhs);

eqns_frac = simplifyFraction(eqns_lhs);
eqns_num = numden(eqns_frac);

fileID = fopen(filename,'w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    fprintf(fileID, string(eqns_num(i)) + '\n');
end
end