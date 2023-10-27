kin = define_easy_1D_unit();

R_06 = eye(3);
p_06 = [1;1;1]*1.25;

[Q, is_LS_vec] = IK.IK_2_intersecting(R_06, p_06, kin);

%%
kin = define_easy_1D_varying();

R_06 = eye(3);
p_06 = [1;1;1]*2;

[Q, is_LS_vec] = IK.IK_2_intersecting(R_06, p_06, kin);

%%
[R_t, P_t] = fwdkin(kin, Q(:,1))


%% Use R for sp5

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
eqns_lhs = simplify(eqns_lhs)

%%
eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)
%%
fileID = fopen('easy_1D_eqns.txt','w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    %fprintf(fileID, string(vpa(eqns_num(i))) + '\n');
    fprintf(fileID, string(eqns_num(i)) + '\n');
end

%%
[Q, is_LS_vec] = IK.IK_2_intersecting(R_06, p_06, kin);
s = [-2.37742, -0.773904, -0.343369, 0.0691691, 0.198726, 1.74226]
xline(2*atan(s), 'r', LineWidth=2);