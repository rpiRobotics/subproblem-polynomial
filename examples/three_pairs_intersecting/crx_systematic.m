kin = define_CRX;

%% Easier example pose
R_06 = eye(3);
p_0T = [0.25; 0.25; 0.25];

%% Harder example pose
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

% H_6 is e_y, so make sure right-most rotation is about e_y

% R_06 = rot(ey, deg2rad(30))*rot(ex, deg2rad(20))*rot(ey, deg2rad(10));
% R_06 = rot(ey, deg2rad(3))*rot(ex, deg2rad(2))*rot(ey, deg2rad(1));
R_06 = rot(ey, deg2rad(3.12))*rot(ex, deg2rad(2.25))*rot(ey, deg2rad(1));

% R_06 = rot(ey, deg2rad(sym(30)))*rot(ex, deg2rad(sym(20)))*rot(ey, deg2rad(sym(10)));
% p_0T = [1415; 9265; 3589]/10000;
p_0T = [3141; 5926; 5358]/10000;


% [R_06, p_0T] = fwdkin(kin, rand_angle([7 1]))
%% Symbolic pose
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

syms a b g
alpha = a
beta= b
gamma = g
R_06 = rot(ey, gamma)*rot(ex, beta)*rot(ey, alpha);

% syms R11 R12 R13 R21 R22 R23 R31 R32 R33  real
% R_06 = [R11 R12 R13; R21 R22 R23; R31 R32 R33]

syms t1 t2 t3 real
p_0T = [t1 t2 t3]' + kin.P(:,1) - R_06*kin.P(:,7);

%% Solve using search-based method
[Q, is_LS_vec] = IK_3_pairs_intersecting(double(R_06),p_0T,kin, true)

%%
[R_test, T_test] = fwdkin(kin, Q(:,1))

%% Alternate approach with Subproblem 5
syms x1 x2 x4 real
R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_34 = half_tan_rot(kin.H(:,4), x4);

p_16 = p_0T - kin.P(:,1) - R_06*kin.P(:,7);
[lhs_1, lhs_2, R_23] = half_tan_sp5_R( ...
    -kin.P(:,3), ...
    R_01'*p_16, ...
    zv, ...
    kin.P(:,5), ...
    -kin.H(:,2), kin.H(:,3), kin.H(:,4), [x2 x4]);
R_04 = R_01 * R_12 * R_23 * R_34;
lhs_err = kin.H(:,5)' * R_04.' * R_06* kin.H(:,6) - kin.H(:,5)'*kin.H(:,6);

eqns_lhs = [lhs_err; lhs_1; lhs_2];
eqns_lhs = simplify(eqns_lhs)


eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)

fileID = fopen('crx_eqns.txt','w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    % fprintf(fileID, string(vpa(eqns_num(i))) + '\n');
    fprintf(fileID, string(eqns_num(i)) + '\n');
end

%%
polynomial_IK.three_pairs_intersecting(kin, R_06, p_0T, 'crx_eqns.txt')

% dr = ResourceFunction["DixonResultant"][{p1, p2, p3, p4}, {x1, x2, x3}]
% Factor[dr]
% pstar = FactorList[dr][[12]][[1]]
% x4 /. NSolve[pstar, x4, Reals]
%%
x_roots = [-3.07428487489949, -1.7429090251912485, -0.6410722219542363, ...
-0.21689616352170946, 0.3252788992213006, 0.5737534119947945, ...
1.5598866488889704, 4.610501097682664]
theta_roots = 2*atan(x_roots);

real_theta = theta_roots(imag(theta_roots)==0)

%%
[Q, is_LS_vec] = IK_3_pairs_intersecting(R_06,p_0T,kin, true)
% xline(real(theta_roots), 'r')
xline(real_theta, 'r')

%% Find remaining joint angles

Q = polynomial_IK.three_pairs_intersecting_given_q4(kin, R_06, p_0T, real_theta)

[R_06_t, p_0T_t] = fwdkin(kin, Q(:,1));
R_06 - R_06_t
p_0T - p_0T_t



