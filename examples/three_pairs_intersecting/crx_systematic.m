kin = define_CRX;

%% Easier example pose
R_06 = eye(3);
p_0T = [0.25; 0.25; 0.25];


%% Example pose with endpoints
R_06 = eye(3);
p_0T = [   0.1; 0.1; 0.3];

%% Harder example pose
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];


[Nb, Db] = rat(tan(deg2rad(225/100)/2), 1e-6)
[Ng, Dg] = rat(tan(deg2rad(312/100)/2), 1e-6)

R_06 = rot(ey, 2*atan(Ng/Dg))*rot(ex, 2*atan(Nb/Db))*rot(ey, deg2rad(1));

p_0T = [3141; 5926; 5358]/10000;

%% harder pose, rational



%% Symbolic pose
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

syms a b g real
alpha = a
beta= b
gamma = g
R_06 = rot(ey, gamma)*rot(ex, beta)*rot(ey, alpha);

% syms R11 R12 R13 R21 R22 R23 R31 R32 R33  real
% R_06 = [R11 R12 R13; R21 R22 R23; R31 R32 R33]

syms t1 t2 t3 real
p_0T = [t1 t2 t3]' + kin.P(:,1) - R_06*kin.P(:,7);


%% Symbolic rational pose

ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

syms xa xb xg
R_06 = half_tan_rot(ey, xg)*half_tan_rot(ex, xb)*half_tan_rot(ey, xa);

syms t1 t2 t3 real
p_0T = [t1 t2 t3]';
%% Solve using search-based method
[Q, is_LS_vec] = IK_3_pairs_intersecting(double(R_06),p_0T,kin, true)
tan(Q(4,:)/2)
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
[eqns_num, plot] = polynomial_IK.three_pairs_intersecting(kin, R_06, p_0T, 'crx_eqns.txt')

% dr = ResourceFunction["DixonResultant"][{p1, p2, p3, p4}, {x1, x2, x3}]
% Factor[dr]
% pstar = FactorList[dr][[12]][[1]]
% x4 /. NSolve[pstar, x4, Reals]

%%
eqns = polynomial_IK.three_pairs_intersecting_indep(kin, R_06, p_0T)
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



