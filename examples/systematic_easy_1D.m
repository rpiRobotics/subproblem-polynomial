ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

kin.joint_type = [0 0 0 0 0 0];
kin.P = [zv ex ex ex ex zv zv];
kin.H = [ez ey ez ey ez ey];

R_06 = eye(3);
p_06 = [1;1;1]*1.25;

[Q, is_LS_vec] = IK.IK_2_intersecting(R_06, p_06, kin);

%%
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

kin.joint_type = [0 0 0 0 0 0];
kin.P = [zv 4*ex 3*ex 2*ex ex zv zv];
kin.H = [ez ey ez ey ez ey];

R_06 = eye(3);
p_06 = [1;1;1]*2;

[Q, is_LS_vec] = IK.IK_2_intersecting(R_06, p_06, kin);

%%
[R_t, P_t] = fwdkin(kin, Q(:,1))

%%
syms x1 x2 x3 x4 real
R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_23 = half_tan_rot(kin.H(:,3), x3);
R_34 = half_tan_rot(kin.H(:,4), x4);

% Error eqn
R_04 = R_01 * R_12 * R_23 * R_34;
lhs_err = kin.H(:,5)' * R_04' * R_06* kin.H(:,6) - kin.H(:,5)'*kin.H(:,6);

% [lhs_1, lhs_2, lhs_3] = half_tan_sp5(-kin.P(:,2), p_06, kin.P(:,3), kin.P(:,4)+R_34*kin.P(:,5), -kin.H(:,1), kin.H(:,2), kin.H(:,3), [x1 x2 x3]);

p_35_3 = kin.P(:,4) + R_34*kin.P(:,5);

[lhs_1, lhs_2, lhs_3] = half_tan_sp5( ...
               -kin.P(:,2), p_06, kin.P(:,3), p_35_3, ...
               -kin.H(:,1), kin.H(:,2), kin.H(:,3), [x1 x2 x3]);

eqns_lhs = [lhs_err; lhs_1; lhs_2; lhs_3];
eqns_lhs = simplify(eqns_lhs)

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


%% Use R for sp5 without tan(theta/2)

syms x1  x3 x4 real
R_01 = half_tan_rot(kin.H(:,1), x1);
% R_12 = half_tan_rot(kin.H(:,2), x2);
R_23 = half_tan_rot(kin.H(:,3), x3);
R_34 = half_tan_rot(kin.H(:,4), x4);

p_35_3 = kin.P(:,4) + R_34*kin.P(:,5);

[lhs_1, lhs_2, alpha_sq, R_12_bar_1, R_12_bar_2] = half_tan_sp5_R_numden( ...
               -kin.P(:,2), p_06, kin.P(:,3), p_35_3, ...
               -kin.H(:,1), kin.H(:,2), kin.H(:,3), [x1 x3]);

% Error eqn
a = kin.H(:,5)' * R_34' * R_23';
b = R_01' * R_06* kin.H(:,6);
c = kin.H(:,5)'*kin.H(:,6);
lhs_err = alpha_sq*(a*R_12_bar_1*b-c)^2-(a*R_12_bar_2*b)^2;


%R_04 = R_01 * R_12 * R_23 * R_34;
%lhs_err = kin.H(:,5)' * R_04' * R_06* kin.H(:,6) - kin.H(:,5)'*kin.H(:,6);


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
%% Test eqns

q = Q(:,1);
x = tan(q/2)

double(subs(eqns_num(3), [x1 x3 x4], x([1 3 4])'))

%%
[Q, is_LS_vec] = IK.IK_2_intersecting(R_06, p_06, kin);
s = [-2.37742, -0.773904, -0.343369, 0.0691691, 0.198726, 1.74226]
xline(2*atan(s), 'r', LineWidth=2);