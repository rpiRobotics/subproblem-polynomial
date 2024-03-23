kin = define_CRX;

R_06 = eye(3);
p_0T = [0.25; 0.25; 0.25];

%% Solve using search-based method
[Q, is_LS_vec] = IK_3_pairs_intersecting(double(R_06),p_0T,kin, true)
tan(Q(4,:)/2)

[R_test, T_test] = fwdkin(kin, Q(:,1))

%% Symbolic rational pose

ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

syms xa xb xg
R_06 = half_tan_rot(ey, xg)*half_tan_rot(ex, xb)*half_tan_rot(ey, xa);

syms t1 t2 t3 real
p_0T = [t1 t2 t3]';

%% Export equations in terms of symbolic pose
[eqns_num, plot] = polynomial_IK.three_pairs_intersecting(kin, R_06, p_0T, 'crx_eqns.txt')


%% Find q4 solutions given roots of polynomial
x_roots = [-3.07428487489949, -1.7429090251912485, -0.6410722219542363, ...
-0.21689616352170946, 0.3252788992213006, 0.5737534119947945, ...
1.5598866488889704, 4.610501097682664]
q4 = 2*atan(x_roots)


%% Find remaining joint angles
Q = polynomial_IK.three_pairs_intersecting_given_q4(kin, R_06, p_0T, q4)

%% Double check forward kinematics
[R_06_t, p_0T_t] = fwdkin(kin, Q(:,1));
R_06 - R_06_t
p_0T - p_0T_t



