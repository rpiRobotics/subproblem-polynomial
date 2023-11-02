%% zero offsets 
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];
kin.H = [ez ey ez ey ez ey];
kin.P = [zv zv ex zv ex zv zv];
kin.joint_type = zeros(1,6);

R_06 = eye(3);
p_0T = [0.25; -1.25; 0.25];
%%
[Q, is_LS_vec] = IK_3_pairs_intersecting(R_06,p_0T,kin, true)
%%
syms R [3 3] real
syms X Y Z real
R_06 = R
p_0T = [X; Y; Z]

%% general 2R-2R-2R robot
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];
kin.H = rand_normal_vec(6);
kin.P = [zv zv rand_vec zv rand_vec zv zv];
kin.joint_type = zeros(1,6);
% 
% q = rand_angle([6 1]);
% [R_06,p_0T] = fwdkin(kin, q)
R_06 = eye(3);
p_0T = [0; -1.5; 0.5];

[Q, is_LS_vec] = IK_3_pairs_intersecting(R_06,p_0T,kin, true)

%% Easy numbers

ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];
kin.H = [ez ey ez ey ez ey];
kin.P = [zv zv rand_vec zv rand_vec zv zv];
kin.P = round(kin.P, 2);
kin.joint_type = zeros(1,6);

R_06 = eye(3);
p_0T = [0; -1.5; 0.5];

[Q, is_LS_vec] = IK_3_pairs_intersecting(R_06,p_0T,kin, true)

%% CRX

kin = define_CRX;
R_06 = eye(3);
p_0T = [0.25; 0.25; 0.25];
% [R_06, p_0T] = fwdkin(kin, rand_angle([7 1]))
%%
[Q, is_LS_vec] = IK_3_pairs_intersecting(R_06,p_0T,kin, true)
%%
syms x1 x2 x3*x4 real
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
[~, a, b, c] = half_tan_sp3(R_34*kin.P(:,5), -kin.P(:,3), kin.H(:,3), norm(p_16));
lhs_x3 = a*x3^2 + b*x3 + c;

% x1 and x2 eqn
[~, ~, a1, b1, c1, a2, b2, c2] = half_tan_sp2(p_16, kin.P(:,3)+R_23*R_34*kin.P(:,5), -kin.H(:,1), kin.H(:,2));

lhs_x1 = a1*x1^2 + b1*x1 + c1;
lhs_x2 = a2*x2^2 + b2*x2 + c2;

eqns_lhs = [lhs_err; lhs_x1; lhs_x2; lhs_x3];
eqns_lhs = simplify(eqns_lhs)
%%
eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)

%%
tic
res = resultant(eqns_num(1), eqns_num(2), x1);
% res = simplify(sqrt(-res), "IgnoreAnalyticConstraints",true);
res = resultant(res, eqns_num(3), x2);
% res = simplify(res/(16*(x3^2 + 1)^4*(x4^2 + 1)^4), "IgnoreAnalyticConstraints",true);
[rem, res] = polynomialReduce(res, 256*(x3^2+1)^4*(x4^2+1)^4);
res = resultant(res, eqns_num(4), x3);
res = simplify(sqrt(res), "IgnoreAnalyticConstraints",true);


c = sym2poly(res);
x_roots = roots(c);
theta_roots = 2*atan(x_roots);

real_theta = theta_roots(imag(theta_roots)==0);
toc
%%
res = resultant(eqns_num(1), eqns_num(3), x2)
[rem, res] = polynomialReduce(res, 4*(x3^2+1)^2);
res = resultant(res, eqns_num(4), x3)
res = resultant(res, eqns_num(2), x1)
[rem, res] = polynomialReduce(res, (x4^2+1)^8);
%%
[r, q] = polynomialReduce(res, (x3^2+1)^4)

%%
factors = factor(res)';

c = sym2poly(factors(3));
x_roots = roots(c);
theta_roots = 2*atan(x_roots);

real_theta = theta_roots(imag(theta_roots)==0);

%%
[Q, is_LS_vec] = IK_3_pairs_intersecting(R_06,p_0T,kin, true)
% xline(real(theta_roots), 'r')
xline(real_theta, 'r')

%% CRX polynomial from Mathematica
P = 592628990976 + 493649640000*x4 - 7228671199758*x4^2 + ...
 8539013655000*x4^3 - 29945708282222*x4^4 - 2896653453750*x4^5 + ...
 45929194539844*x4^6 - 39808558518750*x4^7 + 177776490417945*x4^8 + ...
 39808558518750*x4^9 + 45929194539844*x4^10 + 2896653453750*x4^11 - ...
 29945708282222*x4^12 - 8539013655000*x4^13 - 7228671199758*x4^14 - ...
 493649640000*x4^15 + 592628990976*x4^16;

c = sym2poly(P);
x_roots = roots(c);
theta_roots = 2*atan(x_roots);

real_theta = theta_roots(imag(theta_roots)==0);